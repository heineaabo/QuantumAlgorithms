import numpy as np

import sys
sys.path.append('../..')
sys.path.append('../../optimizers')
sys.path.append('../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian
from fci import FCI

from vqe import VQE
from optimizers import QuantumGradientDescent,Minimizer

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
#for i,circ in enumerate(pairing.circuit_list('vqe')): print(i,'\n',circ)


#options = {'count_states':False,'shots':1000,'print':False}
options = {'shots':1000,
           'optimization_level':1,
           'device':'ibmq_london',
           'layout':[1,0,2,3],
           'noise_model':True,
           'basis_gates':True,
           'coupling_map':True,
           'print':False}
#methods = ['','adam','adagrad','rmsprop']
methods = ['Powell','Nelder-Mead']
maxit = 100
#data = np.zeros((len(methods),maxit))
#theta = 2*np.pi*np.random.randn(1)
theta = np.array([-9.31816914])
for i,method in enumerate(methods):
    vqe = VQE(pairing,Minimizer(method,max_evals=200),'RY',options=options)
    #vqe = VQE(pairing,QuantumGradientDescent(method,max_iter=maxit),'RY',options=options)
    vqe.optimize(theta.copy())
    #data[i] = vqe.optimizer.energies
    np.save('data/optimization/{}_n.npy'.format(method),np.array(vqe.energies))
#np.save('data/optimization/QGD.npy',data)
    


