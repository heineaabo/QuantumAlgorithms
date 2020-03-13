import numpy as np
import matplotlib.pyplot as plt 
from tools import get_H2

import sys
sys.path.append('../..')
sys.path.append('../../../../QuantumCircuitOptimizer')
from vqe import *
#from spsa import SPSA
from qiskit.aqua.components.optimizers import COBYLA,SPSA,POWELL,NELDER_MEAD

from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian

def get_avg(L, num):
    L = np.asarray(L[len(L)-num:])
    return np.mean(L),np.std(L)


l = 4
n = 2

max_iter = 1 #200
max_evals = 1 # 200

bonds = np.linspace(0.5,1.5,11)

result = []
results_evals = []

for i,R in enumerate(bonds):
    print('For bond length {}'.format(R))
    # Get precalculated matrix elements
    h_pq,h_pqrs,Enuc = get_H2(R)

    # Prepare circuit list
    H2 = SecondQuantizedHamiltonian(n,l,h_pq,h_pqrs,
                                    nuclear_repulsion=Enuc,add_spin=True)

    options = {'shots':1000}

    for method in methods:
        print('    - now running {}'.format(method),end='')
        vqe = VQE(H2,options=options)
        theta = vqe.theta
        optimizer = COBYLA(maxiter=max_evals)
        theta,E,_ = optimizer.optimize(len(theta),vqe.expval,initial_point=theta) 
        result.append(E)
        results_evals.append(vqe.energies)
        print(', {} function evaluations'.format(vqe.evals))
        
np.save('data/optimizers/cobyla/200evals/1000shots/final.npy',np.asarray(result))
max_len = max([len(i) for i in results_evals])
evals = np.zeros((len(results_evals),max_len))
for i in range(len(results_evals)):
    num = len(results_evals[i])
    evals[i,num:] = results_evals[i]
np.save('data/optimizers/cobyla/200evals/1000shots/evals.npy',np.asarray(evals))

