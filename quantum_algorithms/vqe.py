import numpy as np
import random

import qiskit as qk
from qiskit.extensions.standard import *

from ansatz import UnitaryCoupledCluster
from attributes import QuantumComputer
from tools import print_state,get_state_count

class VQE:
    def __init__(self,
                 hamiltonian,
                 ansatz = 'UCCSD',
                 ansatz_depth = 1,
                 options = {},
                 callback = {}):
        """
        Input:
            hamiltonian  (int)   - Hamiltonian with system information
                                   and circuit_list.
            ansatz       (class) - Ansatz.
            ansatz_depth (class) - Ansatz depth of Trotter expansion.

            options:
                * seed         (int)  - Seed for simulator and transpiler.
                * shots        (int)  - Number of execution shots.
                * max_energy   (bool) - If maximum or minimum energy is of 
                                        interest.
                * backend      (str)  - Name of backend.
                * device       (str)  - Name of IBMQ quantum computer.
                * noise_model  (bool) - Noise model of device.
                * coupling_map (bool) - Coupling map of device.

                * print        (bool) - Print for every function evaluation.
                * count_states (bool) - Count legal and illegal states for
                                        all function evaluations.
        """
        self.n_qubits = hamiltonian.l
        self.n_fermi = hamiltonian.n
        self.circuit_list = hamiltonian.to_circuit_list(ptype='vqe')
        if ansatz[:3].upper() == 'UCC':
            self.ansatz = UnitaryCoupledCluster(self.n_fermi,
                                                self.n_qubits,
                                                ansatz[3:].upper(),
                                                depth=ansatz_depth)
            self.theta = self.ansatz.new_parameters(hamiltonian.h,
                                                    hamiltonian.v)

        #### Setup options
        self.options = options
        # For execution
        self.shots = 1000 if options.get('shots') == None\
                          else options.get('shots')
        self.seed = options.get('seed')
        if self.seed != None:
            from qiskit.aqua import aqua_globals
            aqua_globals.random_seed = self.seed
        self.max_energy = options.get('max_energy')
        # For optimization
        if options.get('optimizer') != None:
            self.optimizer = options.get('optimizer')
            self.optmizer.set_loss_function(self.expval)
        # For Backend
        if options.get('backend') == None:
            self.options['backend'] = 'qasm_simulator' 
        self.backend = qk.Aer.get_backend(options['backend'])
        # For noise model and coupling map
        self.noise_model, self.coupling_map  = QuantumComputer(options.get('device'),options.get('noise_model'),options.get('coupling_map'))
        # GPU accelerated
        if options.get('gpu'):
            from qiskit_qcgpu_provider import QCGPUProvider
            Provider = QCGPUProvider()
            self.backend = Provider.get_backend(options['backend'])

        self.legal = 0
        self.illegal = 0
        self.evals = 0

        # For plotting after calculations
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
                qc.x(qb[qubit])
                qc.h(qb[qubit])
            elif gate == 'y':
                qc.y(qb[qubit])
                qc.sdg(qb[qubit])
                qc.h(qb[qubit])
                #qc.rx(np.pi/2,qb[qubit])
            elif gate == 'z':
                #qc.z(qb[qubit])
                qc.append(ZGate(),[qubit],[])
        else:
            if gate == 'x':
                qc.cx(qa[0],qb[qubit])
                qc.ch(qa[0],qb[qubit])
            elif gate == 'y':
                qc.cy(qa[0],qb[qubit])
                qc.crx(np.pi/2,qa[0],qb[qubit])
            elif gate == 'z':
                qc.cz(qa[0],qb[qubit])
        return(qc)

    def measure(self,qubit_list,factor,qc,qb,cb):
        if len(qubit_list) == 0:
            return(factor)
        qc.measure(qb,cb)
        job = qk.execute(qc, 
                        backend = self.backend, 
                        shots=self.shots,
                        seed_transpiler=self.seed,
                        seed_simulator=self.seed)
        result = job.result().get_counts(qc)
        E = 0
        #print_state(result,self.ansatz.n,self.ansatz.l)
        if self.options.get('count_states'):
            legal,illegal = get_state_count(result,self.ansatz.n,self.ansatz.l)
            self.legal += legal
            self.illegal += illegal
        for key, value in result.items():
            key1 = key[::-1]
            eigenval = 1
            for idx in qubit_list:
                e =  1 if key1[idx] == '0' else -1
                eigenval *= e
            E += eigenval*value
        E /= self.shots
        return(factor*E)

    def expval(self,theta=None):
        if theta == None:
            theta = self.theta
        E = 0
        qc,qb,cb = None,None,None
        for i,pauli_string in enumerate(self.circuit_list):
            factor = pauli_string[0].real
            qubit_list = []
            qb,qc,cb =self.initialize_circuit()
            qc = self.ansatz(theta,qc,qb)
            for qubit,gate in pauli_string[1:]:
                qc = self.add_gate(qubit,gate,qc,qb)
                qubit_list.append(qubit)
            #print('Term {}: {}'.format(i,pauli_string))
            E += self.measure(qubit_list,factor,qc,qb,cb)
        if self.options.get('print'):
            print('<E> = ', E)
        if self.options.get('max_energy'):
            E = -E
        self.energies.append(E)
        self.evals += 1
        return(E)

    def optimize_classical(self,
                           theta,
                           method='L-BFGS-B', # Minimization method.
                           max_iters = 200,   # Minimizer iterations.
                           max_eval = 200,    # Funtion evaluations.
                           tol=1e-08,
                           adaptive=False):   # Nelder mead to be adaptive
        from scipy.optimize import minimize
        self.energies = []
        max_eval_str = 'maxfev' # Different string notation for different methods
        bounds = None
        if method == 'L-BFGS-B':
            bounds = [(0,2*np.pi) for i in theta]
            max_eval_str = 'maxfun'
        # Set scipy.minimize options
        options={'disp':True,
                 'maxiter':max_iters,
                 max_eval_str:max_eval}
        if method == 'Nelder-Mead':
            options['adaptive'] = adaptive
        result = minimize(self.expval,
                            theta,
                            bounds = bounds,
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
                    # Prepare ancilla qubit with hadamard
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

