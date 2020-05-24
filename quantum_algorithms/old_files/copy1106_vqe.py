import numpy as np
import random

import qiskit as qk
from qiskit.extensions.standard import *

from ansatz import UnitaryCoupledCluster
from attributes import QuantumComputer
from tools import print_state,get_state_count
from algorithm import QuantumAlgorithm

class VQE(QuantumAlgorithm):
    def __init__(self,
                 hamiltonian,
                 #optimizer,
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
        super().__init__(options)
        if isinstance(hamiltonian,dict):
            self.n_qubits = hamiltonian['l']
            self.n_fermi = hamiltonian['n']
            self.circuit_list = hamiltonian['circuit']
            self.conv = 1 if hamiltonian.get('conv') == None else hamiltonian['conv'] # Occupation convention
        else:
            self.n_qubits = hamiltonian.l
            self.n_fermi = hamiltonian.n
            self.circuit_list = hamiltonian.to_circuit_list(ptype='vqe')
            self.conv = hamiltonian.conv # Occupation convention

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
        # Custom ansatz
        else:
            self.ansatz = ansatz
            self.theta = ansatz.theta

        # For optimization
        if options.get('optimizer') != None:
            self.optimizer = options.get('optimizer')
            self.optmizer.set_loss_function(self.expval)

        # For counting states
        self.legal = 0
        self.illegal = 0
        self.evals = 0

        # For plotting progression from optimization
        self.energies = []

    def initialize_circuit(self):
        n_qubits = self.n_qubits
        qb = qk.QuantumRegister(n_qubits)
        cb = qk.ClassicalRegister(n_qubits)
        qc = qk.QuantumCircuit(qb,cb)
        return(qb,qc,cb)

    def add_gate(self,qubit,gate,qc,qb,qa=None):
        """
        Add gate with measurement transformation.
        If qa register -> controlled operations.
        """
        if qa == None:
            if gate == 'x':
                qc.h(qb[qubit])
                #qc.x(qb[qubit])
            elif gate == 'y':
                #qc.rx(-np.pi/2,qb[qubit])
                qc.sdg(qb[qubit])
                qc.h(qb[qubit])
                #qc.y(qb[qubit])
            #elif gate == 'z':
            #    qc.z(qb[qubit])
        else:
            if gate == 'x':
                #qc.cx(qa[0],qb[qubit])
                qc.ch(qa[0],qb[qubit])
            elif gate == 'y':
                #qc.cy(qa[0],qb[qubit])
                qc.crx(-np.pi/2,qa[0],qb[qubit])
            #elif gate == 'z':
            #    qc.cz(qa[0],qb[qubit])
        return(qc)

    def prepare_pauli(self,pauli_string,qc,qb,ancilla=True):
        """
        Add gate with measurement transformation.
        If qa register -> controlled operations.
        """
        qubit_list = []
        for qubit,gate in pauli_string[1:]:
            if gate == 'x':
                qc.h(qb[qubit])
            elif gate == 'y':
                qc.sdg(qb[qubit])
                qc.h(qb[qubit])
            qubit_list.append(qubit)
        target = qubit_list[0]
        if ancilla:
            for qbit in qubit_list[1:]:
                qc.cx(qb[qbit],qb[target])
            return qc,target
        else:
            return qc,qubit_list

    def measure(self,target,factor,qc,qb,cb):
        qc.measure(qb,cb)
        job = qk.execute(qc, 
                        backend = self.backend, 
                        shots=self.shots,
                        noise_model=self.noise_model,
                        coupling_map=self.coupling_map,
                        basis_gates=self.basis_gates,
                        seed_transpiler=self.seed,
                        seed_simulator=self.seed)
        result = job.result().get_counts(qc)
        E = 0
        if self.options.get('count_states'):
            legal,illegal = get_state_count(result,self.ansatz.n,self.ansatz.l)
            self.legal += legal
            self.illegal += illegal
        for state,num_measure in result.items():
            state = state[::-1]
            if isinstance(target,int):
                eigenval = 1 if state[target] == '0' else -1
            elif isinstance(target,(list,tuple)):
                eigenval = 1
                for i in target:
                    eigenval *= 1 if state[i] == '0' else -1
            E += eigenval*num_measure
        E /= self.shots
        return factor*E

    def expval(self,theta=None):
        if theta is None:
            theta = self.theta
        E = 0
        qb,qc,cb = self.initialize_circuit()
        qc_ansatz = self.ansatz(theta,qc,qb)
        Z_list = []
        for i,pauli_string in enumerate(self.circuit_list):
            factor = pauli_string[0].real
            if len(pauli_string) == 1:
                E += factor
                continue
            only_Z = True
            for action in pauli_string[1:]:
                if action[1] in ['x','y']:
                    only_Z = False
            if only_Z:
                Z_list.append(pauli_string)
                continue
            qc = qk.QuantumCircuit(qb,cb)
            qc,target = self.prepare_pauli(pauli_string,qc,qb,ancilla=self.ancilla_measure)
            qc = qc_ansatz + qc
            E += self.measure(target,factor,qc,qb,cb)
        E += self.Z_expval(Z_list,qc_ansatz,qb,cb)
        if self.prnt:
            print('<E> = ', E,', theta =',theta)
        if self.options.get('max_energy'):
            E = -E
        self.energies.append(E)
        self.evals += 1
        return E

    def Z_expval(self,Z_list,qc,qb,cb):
        if len(Z_list) == 0:
            return 0
        result = self.Z_measure(qc,qb,cb)
        all_E = 0
        for i,pauli_string in enumerate(Z_list):
            factor = pauli_string[0]
            Zs = [action[0] for action in pauli_string[1:]]
            E = 0
            for key,val in result.items():
                state = key[::-1]
                e = 1
                for qbit in Zs:
                    e *= 1 if state[qbit] == '0' else -1
                E += e*val
            E /= self.shots
            all_E += factor*E
        return all_E

    def Z_measure(self,qc,qb,cb):
        qc.measure(qb,cb)
        job = qk.execute(qc, 
                        backend = self.backend, 
                        shots=self.shots,
                        noise_model=self.noise_model,
                        coupling_map=self.coupling_map,
                        basis_gates=self.basis_gates,
                        seed_transpiler=self.seed,
                        seed_simulator=self.seed)
        result = job.result().get_counts(qc)
        return result

    def optimize_classical(self,
                           method='L-BFGS-B', # Minimization method.
                           max_iters = 200,   # Minimizer iterations.
                           max_eval = 200,    # Funtion evaluations.
                           tol=1e-08,
                           adaptive=False):   # Nelder mead to be adaptive
        from scipy.optimize import minimize
        theta = self.theta
        self.energies = []
        max_eval_str = 'maxfev' # Different string notation for different methods
        bounds = None
        if method == 'L-BFGS-B':
            bounds = [(0,2*np.pi) for i in theta]
            max_eval_str = 'maxfun'
        # Set scipy.minimize options
        options= None
        #options={'disp':True,
        #         'maxiter':max_iters,
        #         max_eval_str:max_eval}
        if method == 'Nelder-Mead':
            options['adaptive'] = adaptive
        result = minimize(self.expval,
                            theta,
                            bounds=bounds,
                            method=method,
                            options=options,
                            tol=tol)
        params = result.x
        return params

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
                    qb,qc,cb =self.initialize_circuit()
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

