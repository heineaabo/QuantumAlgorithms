import numpy as np
from matrix import pairing_mat_elems

import sys
sys.path.append('../..')
from vqe import VQE
sys.path.append('../../optimizers')
from optimizers import QuantumGradientDescent
sys.path.append('../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian

import matplotlib.pyplot as plt

n = 2
l = 4
delta = 1
g = 1

h,v = pairing_mat_elems(n,l,delta,g)


pairing =  PairingHamiltonian(n,l,h,v)


options = {'shots':1000,
           'optimization_level':1,
           #'device':'ibmq_london',
           #'layout':[1,0,2,3],
           #'noise_model':True,
           #'basis_gates':True,
           #'coupling_map':True,
           #'seed':1,
           'print':False}
method = 'rmsprop'
maxit = 100
tols = [1e-04,1e-05]

avgs = [2,5,10,20]

theta = np.array([-9.31816914]) # Bad theta value
for j,tol in enumerate(tols):
    print('For tolerance',tol)
    vqe = VQE(pairing,
              QuantumGradientDescent(method,
                                     max_iter=maxit,
                                     tol=tol),
              'RY',
              options=options)
    vqe.optimize(theta.copy())
    Es = vqe.optimizer.energies
    mean = [np.mean(Es[-i:]) for i in avgs]
    var  = [np.var(Es[-i:]) for i in avgs]
    print('Num    Mean        Var')
    for i,a in enumerate(avgs):
        print('{:<8}{:<12.8f}{:.8f}'.format(a,mean[i],var[i]))
    #plt.figure(j)
    plt.plot(vqe.optimizer.thetas,label='tol {}'.format(tol))
plt.legend()
plt.show()
