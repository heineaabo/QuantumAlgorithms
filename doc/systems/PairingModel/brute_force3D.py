import numpy as np
import matplotlib.pyplot as plt 

import sys
sys.path.append('../..')
sys.path.append('../../../../QuantumCircuitOptimizer')
from vqe import *
from ucc import UnitaryCoupledCluster
#from spsa import SPSA

from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian

from tqdm import tqdm


l = 5
n = 2

delta = 1 # Level spacing
g = 1     # Interaction strength
#g /= 4

# Matrix elements
h_pq = np.identity(l)
for p in range(l):
	h_pq[p,p] *= delta*(p - (p%2))/2
	
h_pqrs = np.zeros((l,l,l,l))
for p in range(0,l-1,2):
	for r in range(0,l-1,2):
		h_pqrs[p,p+1,r,r+1] = -0.5*g

# Prepare circuit list
pairing = SecondQuantizedHamiltonian(n,l)
pairing.set_integrals(h_pq,h_pqrs)
pairing.get_circuit()
circuit_list = pairing.to_circuit_list(ptype='vqe')

ansatz = UnitaryCoupledCluster(n,l,'D',1)
og_params = ansatz.new_parameters(pairing.h,
                                  pairing.v)
print(len(og_params))

grid_points = 15
x = np.linspace(0,2*np.pi,grid_points)
y = np.linspace(0,2*np.pi,grid_points)
z = np.linspace(0,2*np.pi,grid_points)
params = []
for i,xi in enumerate(x):
    for j,yi in enumerate(y):
        for k,zi in enumerate(z):
            params.append([xi,yi,zi])
Es = np.zeros(grid_points*grid_points*grid_points)

theta = og_params
vqe = VQE(n_qubits = l,
    ansatz = ansatz,
    circuit_list = circuit_list,
    shots = 1000,
    ancilla=0,
    max_energy=False,
    prnt=False)

assert len(params) == grid_points*grid_points*grid_points
i = 0
for theta in tqdm(params):
    Es[i] = vqe.expval(theta)
    i += 1

np.save('data/brute/3D/grid.npy',x)
np.save('data/brute/3D/parameters.npy',np.asarray(params))
np.save('data/brute/3D/values1000.npy',Es)

