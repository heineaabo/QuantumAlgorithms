import numpy as np
import qiskit as qk
from qiskit.extensions.standard import *
from scipy.optimize import minimize
import random

from tools import print_state

class VQE:
    def __init__(self,
                n_qubits = None,
                n_fermi = None,
                ansatz = None,
                circuit_list = None,
                optimizer = None,
                shots = 1000,
                seed=None,
                ancilla=0,
                max_energy=False,
                prnt=False):
        self.circuit_list = circuit_list
        self.ansatz = ansatz
        self.n_qubits = n_qubits
        self.n_fermi = n_fermi
        self.shots = shots
        self.seed = seed
        self.ancilla = ancilla
        self.max_energy = max_energy
        self.prnt=prnt

        # Arbitrary optimization class
        if optimizer != None:
            self.optimizer = optimizer
            self.optmizer.set_loss_function(self.expval)

        # For plotting after calculations
        self.energies = []

    def initialize_circuit(self):
        n_qubits = self.n_qubits
        qb = qk.QuantumRegister(n_qubits+self.ancilla)
        cb = qk.ClassicalRegister(n_qubits+self.ancilla)
        qc = qk.QuantumCircuit(qb,cb)
        return(qb,qc,cb)

    def add_gate(self,qubit,gate,qc):
        if gate == 'x':
            qc.append(XGate(),[qubit],[])
            qc.append(HGate(),[qubit],[])
        elif gate == 'y':
            qc.append(YGate(),[qubit],[])
            qc.append(SdgGate(),[qubit],[])
            qc.append(HGate(),[qubit],[])
        elif gate == 'z':
            qc.append(ZGate(),[qubit],[])
        return(qc)

    def measure(self,qubit_list,factor,qc,qb,cb):
        if len(qubit_list) == 0:
            return(factor)
        qc.measure(qb,cb)
        job = qk.execute(qc, 
                        backend = qk.Aer.get_backend('qasm_simulator',seed=self.seed), 
                        shots=self.shots)
        result = job.result().get_counts(qc)
        E = 0
        #print_state(result,self.ansatz.n,self.ansatz.l)
        for key, value in result.items():
            key1 = key[::-1]
            eigenval = 1
            for idx in qubit_list:
                e =  1 if key1[idx] == '0' else -1
                eigenval *= e
            E += eigenval*value
        E /= self.shots
        return(factor*E)

    def expval(self,theta):
        E = 0
        qc,qb,cb = None,None,None
        for pauli_string in self.circuit_list:
            factor = pauli_string[0].real
            qubit_list = []
            qb,qc,cb =self.initialize_circuit()
            qc = self.ansatz(theta,qc,qb)
            for qubit,gate in pauli_string[1:]:
                qc = self.add_gate(qubit,gate,qc)
                qubit_list.append(qubit)
            E += self.measure(qubit_list,factor,qc,qb,cb)
        if self.prnt:
            print('<E> = ', E)
        if self.max_energy:
            E = -E
        self.energies.append(E)
        return(E)

    def optimize_classical(self,theta,method='L-BFGS-B',max_iters = 1000):
        self.energies = []
        if method == 'L-BFGS-B':
            bounds = [(0,2*np.pi) for i in theta]
        else:
            bounds = None
        result = minimize(self.expval,
                            theta,
                            bounds = bounds,
                            method=method,options={'disp':True,'maxiter':max_iters})
        params = result.x
        return params

    def gradient_descent(self,theta,max_iters=100,dt=1e-2):
        # Iteration
        # Calculate for original thetas
        og_exp = expcal(theta,prnt=False)
        expvals = np.zeros_like(theta)
        for i,t in enumerate(theta):
            theta[i] = t+dt
            expvals[i] = expcal(theta,prnt=False)
            theta[i] = t
        grad = np.ones_like(theta)*og_exp - expvals




