import numpy as np

import sys
sys.path.append('../..')
sys.path.append('../../optimizers')
sys.path.append('../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,FermionHamiltonian

from optimizer import Minimizer

from vqe import VQE

l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals

h_pq = np.load('matrix/h_n{}_l{}.npy'.format(n,l))
h_pqrs = np.load('matrix/v_n{}_l{}.npy'.format(n,l))
E = np.load('matrix/E_n{}_l{}.npy'.format(n,l))

from time import time
#t1 = time()
#qdot = FermionHamiltonian(n,l,h_pq,h_pqrs,anti_symmetric=True)
#t2 = time()
t3 = time()
qdot2 = SecondQuantizedHamiltonian(n,l,h_pq,h_pqrs,anti_symmetric=True)
t4 = time()

#for i in qdot.circuit_list: print(i)
#print('NEXT')
print('Time:',t4-t3)
for i in qdot2.circuit_list('vqe'): print(i)


options = {'shots':10000,'print':True}
#options = {'shots':1000,
#           'print':True,
#           'device':'ibmq_16_melbourne',
#           'noise_model':True} #,
#           'coupling_map':True}

#model = VQE(qdot,Minimizer('powell'),'UCCSD',options=options)

#param = model.optimize()
#print(param)

#print('FCI:',E)
#print(model.expval())

