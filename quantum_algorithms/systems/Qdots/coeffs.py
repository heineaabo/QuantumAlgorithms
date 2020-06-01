import numpy as np
from matrix import get_qdot_matrix

import sys
sys.path.append('../..')
from fci import FCI
sys.path.append('../../optimizers')
sys.path.append('../../../../QuantumCircuitOptimizer')
from quantum_circuit import SecondQuantizedHamiltonian

from optimizer import Minimizer

from vqe import VQE

l = 8     # number of spin orbitals / number of qubits
n = 4     # Number of occupied spin orbitals

omega = 1

h,v,E = get_qdot_matrix(n,l,omega)

eig,sta = FCI(n,l,h,v,ret_all=True)
ground = sta[:,0]
ground[np.absolute(ground) < 1e-8 ] = 0
#print(ground)

from itertools import combinations
states = []
states_int = []
for state in combinations(range(l),n):
    states_int.append([orb for orb in state])
states_str = [['0' for j in range(l)] for i in states_int]
for state,inds in zip(states_str,states_int):
    for ind in inds:
        state[ind] = '1'
    states.append(''.join(state))
for i,j in zip(ground,states):
    if i != 0:
        print(j,i)
