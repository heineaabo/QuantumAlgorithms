import numpy as np

import sys
sys.path.append('../..')
sys.path.append('../../optimizers')
sys.path.append('../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,FermionHamiltonian
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

#theta = [5.829889373194686] # Hardcode good parameter

import time
t1 = time.time()
pairing =  SecondQuantizedHamiltonian(n,l,h_pq,h_pqrs,anti_symmetric=True)
t2 = time.time()
print('Time:',t2-t1)

for i in pairing.circuit_list('vqe'): print(i)

#options = {'count_states':False,'shots':100000,'print':False,'seed':5,'ancilla':False}
#options = {'shots':10000,'print':True}
#options = {'shots':1000,
#           'print':True,
#           'device':'ibmq_16_melbourne',
#           'noise_model':True} #,
#           'coupling_map':True}
#model = VQE(pairing,Minimizer('Powell'),'UCCD',options=options)
#for i in model.circuit_list: print(i)

#param = model.optimize()
#print(param)

#print(model.theta)
#print(model.expval(theta))
