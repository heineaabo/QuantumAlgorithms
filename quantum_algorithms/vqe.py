import numpy as np
import qiskit as qk
from qiskit.extensions.standard import *
from tools import print_state,get_state_count
from ansatz import UnitaryCoupledCluster,RYRZ
from algorithm import QuantumAlgorithm
from attributes import QuantumComputer

class VQE(QuantumAlgorithm):
    def __init__(self,
                 hamiltonian,
                 optimizer,
                 ansatz = 'UCCSD',
                 ansatz_depth = 1,
                 options = {}):
        """
        Input:
            hamiltonian  (class) - Hamiltonian with system information
                                   and circuit_list.
            optimizer    (class) - Optimization method.
            ansatz       (class) - Ansatz.
            ansatz_depth (class) - Ansatz depth of Trotter expansion.

        """
        if isinstance(hamiltonian,dict):
            self.n_qubits = hamiltonian['l']
            self.n_fermi = hamiltonian['n']
            self.circuit_list = hamiltonian['circuit']
            self.conv = 1 if hamiltonian.get('conv') == None else hamiltonian['conv'] # Occupation convention
        else:
            self.n_qubits = hamiltonian.l
            self.n_fermi = hamiltonian.n
            self.circuit_list = hamiltonian.circuit_list('vqe')
            self.conv = hamiltonian.conv # Occupation convention

        super().__init__(self.n_qubits,options)
        # Unitary Coupled Cluster ansatz
        if isinstance(ansatz,str):
            if ansatz[:3].upper() == 'UCC':
                self.ansatz = UnitaryCoupledCluster(self.n_fermi,
                                                    self.n_qubits,
                                                    ansatz[3:].upper(),
                                                    depth=ansatz_depth)
                if isinstance(hamiltonian,dict):
                    self.theta = hamiltonian.theta
                else:
                    self.theta = self.ansatz.new_parameters(hamiltonian.h,
                                                            hamiltonian.v)
            elif ansatz.upper() == 'RYRZ':
                self.ansatz = RYRZ(self.n_fermi,self.n_qubits,depth=ansatz_depth)
                self.theta = self.ansatz.new_parameters()

        # Custom ansatz
        else:
            self.ansatz = ansatz
            self.theta = 0

        # For optimization
        self.optimizer = optimizer
        self.optimizer.set_loss_function(self.expval)

        # For counting states
        self.legal = 0
        self.illegal = 0
        self.evals = 0

        # For plotting progression from optimization
        self.energies = []

    def expval(self,theta=None):
        if theta is None:
            theta = self.theta
        E = 0
        # Prepare qiskit circuit and registers.
        qb = qk.QuantumRegister(self.n_qubits)
        cb = qk.ClassicalRegister(self.n_qubits)
        qc = qk.QuantumCircuit(qb,cb)
        # Prepare ansatz to be reused.
        qc_ansatz = self.ansatz(theta,qc,qb)
        for i,pauli_string in enumerate(self.circuit_list):
            # Skip constant term
            if len(pauli_string) == 0:
                E += pauli_string.factor
                continue
            # New circuit
            qc = qk.QuantumCircuit(qb,cb)
            # Apply measurement transformation
            qc = pauli_string.prepare(qc,qb) 
            # Combine circuit with ansatz circuit
            qc = qc_ansatz + qc
            # Measure circuit and add expectation value
            measurement = self.measure(qc,qb,cb)
            if self.meas_fitter != None:
                measurement = self.meas_fitter.filter.apply(measurement)
            E += pauli_string.expectation(measurement,self.shots)
        if self.prnt:
            #print('⟨E⟩ = {}, θ = {}'.format(E,theta))
            print('⟨E⟩ = {}'.format(E))
        self.energies.append(E)
        self.evals += 1
        return E

    def optimize(self,theta=None):
        if theta == None:
            theta = self.theta
        return self.optimizer(theta)

    def optimize_gradient(self,
                          theta,
                          max_iters=200,
                          max_evals=200,
                          step=1e-02, # Step length
                          tol=1e-08):
        self.energies = []
        qc,qb,qa,cb,ca = None,None,None,None,None
        new_theta = theta
        for i in range(max_iters):
            E_inter = self.expval(new_theta)
            print('<E> =',E_inter)
            for d_j in range(len(theta)):
                im = 0
                for pauli_string in self.circuit_list:
                    factor = pauli_string[0].real
                    qubit_list = []
                    qb = qk.QuantumRegister(self.n_qubits)
                    cb = qk.ClassicalRegister(self.n_qubits)
                    qc = qk.QuantumCircuit(qb,cb)
                    # Prepare ancilla qubit with Hadamard
                    qa = qk.QuantumRegister(1)
                    ca = qk.ClassicalRegister(1)
                    qc.add_register(qa,ca)
                    qc.h(qa[0])
                    qc = self.ansatz(new_theta,qc,qb,qa,d_j)
                    for qubit,gate in pauli_string[1:]:
                        qc = self.add_gate(qubit,gate,qc,qb,qa)
                        qubit_list.append(qubit)
                    qc.h(qa[0])
                    qc.rx(np.pi/2,qa[0])
                    im += self.ancilla_measure(factor,qc,qa,ca)
                print('Updating theta {} -> {} - {} = {}'.format(d_j,new_theta[d_j],im,new_theta[d_j]-im*step))
                new_theta[d_j] -= im*step
            print('Iteration {} finished!'.format(i),new_theta)
        return theta

    def ancilla_measure(self,factor,qc,qa,ca):
        qc.measure(qa,ca)
        job = qk.execute(qc, 
                        backend = self.backend, 
                        shots=self.shots,
                        seed_transpiler=self.seed,
                        seed_simulator=self.seed)
        result = job.result().get_counts(qc)
        im = 0
        for key, value in result.items():
            if key[0] == '0':
                im += 1 - 2*(value/self.shots)
            elif key[0] == '1':
                im += 2*(value/self.shots) - 1
        return factor*im

