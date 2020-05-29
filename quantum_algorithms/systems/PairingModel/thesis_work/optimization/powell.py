import numpy as np

import sys
sys.path.append('../../../..')
from vqe import VQE
from fci import FCI

sys.path.append('../../../../optimizers')
from optimizer import Minimizer

sys.path.append('../../../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian

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


Efci = FCI(n,l,h_pq,h_pqrs)
print('FCI energy :',Efci)


pairing =  PairingHamiltonian(n,l,h_pq,h_pqrs)

options = {'shots':1000,'print':True}
#options = {'shots':1000,
#           'optimization_level':1,
#           'device':'ibmq_london',
#           'layout':[0,3,2,1],  
#           #'noise_model':True,
#           'basis_gates':True,
#           'coupling_map':True,
#           'print':True}
model = VQE(pairing,Minimizer('Powell',tol=1e-03),'RYPAIRING',options=options)
param = model.optimize()
param = [param]
num = 10
Es = np.zeros(num)
for i in range(num):
    Es[i] = model.expval(param)
print('E =',np.mean(Es),'+-',np.std(Es))

