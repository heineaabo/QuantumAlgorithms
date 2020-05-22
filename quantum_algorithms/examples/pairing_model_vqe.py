import numpy as np
import sys
sys.path.append('..')
sys.path.append('../../../QuantumCircuitOptimizer')
from vqe import VQE
from quantum_circuit import SecondQuantizedHamiltonian


l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals
delta = 1 # Level spacing
g = 1     # Interaction strength

# Matrix elements
h_pq = np.identity(l)
for p in range(l):
	h_pq[p,p] *= delta*(p - (p%2))/2
	
h_pqrs = np.zeros((l,l,l,l))
for p in range(0,l-1,2):
	for r in range(0,l-1,2):
		h_pqrs[p,p+1,r,r+1] = -0.5*g


options = {'seed':1}
pairing_model = SecondQuantizedHamiltonian(n,l,h_pq,h_pqrs)
vqe = VQE(pairing_model,ansatz='UCCD',options=options)
print('Energy expectation value with MP2 parameters:',vqe.expval())
