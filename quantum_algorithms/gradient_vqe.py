from vqe import VQE
import numpy as np
import qiskit as qk

def gradient(self,theta):
    qb = qk.QuantumRegister(self.n_qubits)
    cb = qk.ClassicalRegister(self.n_qubits)
    grad = np.zeros_like(theta)
    for d_j in range(len(theta)):
        im = 0
        for pauli_string in self.circuit_list:
            qubit_list = []
            # New circuit
            qc = qk.QuantumCircuit(qb,cb)
            # Prepare ancilla qubit with Hadamard
            qa = qk.QuantumRegister(1)
            ca = qk.ClassicalRegister(1)
            qc.add_register(qa,ca)
            qc.h(qa[0])
            # Apply measurement transformation
            qc = pauli_string.prepare(qc,qb,qa=qa) 
            # Finish Hadamard test
            qc.h(qa[0])
            qc.rx(np.pi/2,qa[0])
            # Measure circuit and add expectation value
            measurement = self.measure(qc,qa,ca)
            im += self.ancilla_expectation(measurement)
        #print('Updating theta {} -> {} - {} = {}'.format(d_j,new_theta[d_j],im,new_theta[d_j]-im*step))
        grad[d_j] += im*2
    return grad
VQE.gradient = gradient

def simple_gradient(self,theta):
    grad = np.zeros_like(theta)
    for k in range(len(theta)):
        e_k = np.zeros_like(theta)
        e_k[k] = 0.5*np.pi
        left = self.expval(theta+e_k)
        right = self.expval(theta-e_k)
        grad[k] = 0.5*(right - left)
    return grad
VQE.simple_gradient = simple_gradient

def ancilla_expectation(self,measurement):
    im = 0
    for key, value in measurement.items():
        if key[0] == '0':
            im += 1 - 2*(value/self.shots)
        elif key[0] == '1':
            im += 2*(value/self.shots) - 1
    return im
VQE.ancilla_expectation = ancilla_expectation
