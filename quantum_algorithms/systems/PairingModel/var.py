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
           'device':'ibmq_london',
           'layout':[1,0,2,3],
           'noise_model':True,
           'basis_gates':True,
           'coupling_map':True,
           #'seed':1,
           'print':False}
method = 'rmsprop'
maxit = 100
tols = [1e-04,1e-05]

avgs = [2,5,10,20]

theta = np.array([-9.31816914]) # Bad theta value
vqe = VQE(pairing,
          QuantumGradientDescent(method,
                                 max_iter=maxit,
                                 step_length=1,
                                 tol=1e-05), #tol),
          'RY',
          options=options)
vqe.optimize(theta)
plt.plot(vqe.optimizer.thetas)
plt.show()
#Es = []
#for i in range(20):
#    Es.append(vqe.expval(theta.copy()))
#print('For 2:')
#print('Mean:',np.mean(Es[:2]))
#print('Variance:',np.var(Es[:2]))
#print('')
#print('For 5:')
#print('Mean:',np.mean(Es[:5]))
#print('Variance:',np.var(Es[:5]))
#print('')
#print('For 10:')
#print('Mean:',np.mean(Es[:10]))
#print('Variance:',np.var(Es[:10]))
#print('')
#print('For 20:')
#print('Mean:',np.mean(Es[:20]))
#print('Variance:',np.var(Es[:20]))
#print('')
