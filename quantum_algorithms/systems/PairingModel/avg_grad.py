import numpy as np
from matrix import pairing_mat_elems

import sys
sys.path.append('../..')
from vqe import VQE
sys.path.append('../../optimizers')
from optimizers import QuantumGradientDescent
sys.path.append('../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian

n = 2
l = 4
delta = 1
g = 1

h,v = pairing_mat_elems(n,l,delta,g)


pairing =  PairingHamiltonian(n,l,h,v)


options = {'shots':1000,
           'optimization_level':1,
           'device':'ibmq_london',
           'layout':[1,0,2,3],
           'noise_model':True,
           'basis_gates':True,
           'coupling_map':True,
           #'seed':1,
           'print':False}
methods = ['','adam','adagrad','rmsprop']
maxit = 60
avg_nums = [0,2,5,10]

theta = np.array([-9.31816914]) # Bad theta value
for method in methods:
    data = np.zeros((len(avg_nums),maxit))
    for j,avg_num in enumerate(avg_nums):
        vqe = VQE(pairing,
                  QuantumGradientDescent(method,
                                         max_iter=maxit,
                                         avg_num=avg_num),
                  'RY',
                  options=options)
        vqe.optimize(theta.copy())
        data[j] = vqe.optimizer.energies
    if method == '':
        method = 'GD'
    np.save('data/optimization/avg_grad/{}.npy'.format(method),data)
