import numpy as np

from tqdm import tqdm

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


pairing =  PairingHamiltonian(n,l,h_pq,h_pqrs)

#options = {'count_states':False,'shots':1000,'print':False}
options = {'shots':1000,
           'optimization_level':1,
           'device':'ibmq_london',
           'layout':[1,0,2,3],
           'noise_model':True,
           'basis_gates':True,
           'coupling_map':True,
           'print':False}
#model = VQE(pairing,Minimizer('Powell'),'RYRZ',ansatz_depth=10,options=options)



# Get grid
grid_points = 50
x = np.linspace(0,2*np.pi,grid_points)
y = np.linspace(0,2*np.pi,grid_points)
params = []
for i,xi in enumerate(x):
    for j,yi in enumerate(y):
        params.append([xi,yi])
Es = np.zeros(grid_points*grid_points)

vqe = VQE(pairing,Minimizer('Powell'),'RYRZ',options=options)

og_params = vqe.theta
assert len(params) == grid_points*grid_points
i = 0
for theta in tqdm(params):
    #print('{} / {}'.format(i+1,grid_points*grid_points))
    Es[i] = vqe.expval(theta)
    i += 1

np.save('data/brute/ryrz/grid.npy',x)
np.save('data/brute/ryrz/parameters.npy',np.asarray(params))
np.save('data/brute/ryrz/values1000.npy',Es)
