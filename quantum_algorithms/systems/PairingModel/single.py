import numpy as np

import sys
sys.path.append('../..')
sys.path.append('../../optimizers')
sys.path.append('../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian
from fci import FCI

from vqe import VQE
from optimizer import Minimizer

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


import time
t1 = time.time()
pairing =  PairingHamiltonian(n,l,h_pq,h_pqrs)
t2 = time.time()
print('Time:',t2-t1)

pairing.group_paulis(qwc=False,gc=True)



theta = [5.829889373194686] # Hardcode good parameter


#options = {'count_states':False,'shots':1000,'print':True}
#options = {'shots':10000,'print':True}
options = {'shots':1000,
           'seed':1,
           #'optimization_level':1,
           #'device':'ibmq_london',
           #'meas_fit':True,
           #'layout':[0,3,2,1],  
           #'noise_model':True,
           #'basis_gates':True,
           #'coupling_map':True,
           'print':True}
model = VQE(pairing,Minimizer('Powell'),'RYPAIRING',options=options)

#param = model.optimize()

#print(model.theta)
print(model.expval(theta))
