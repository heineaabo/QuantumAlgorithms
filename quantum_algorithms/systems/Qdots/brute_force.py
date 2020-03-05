import numpy as np
import matplotlib.pyplot as plt 

import sys
sys.path.append('../..')
sys.path.append('../../../../QuantumCircuitOptimizer')
from vqe import *
from ucc import UnitaryCoupledCluster
#from spsa import SPSA

from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian

from quantum_systems import ODQD
from quantum_systems.quantum_dots.one_dim.one_dim_potentials import HOPotential
from coupled_cluster import CCSD


l = 4
n = 2


# Generate system
omega = 1

grid_length = 5
num_grid_points = 1001

odho = ODQD(n,l,grid_length,num_grid_points)
odho.setup_system(potential=HOPotential(omega),add_spin=True)

one_body = odho.h
one_body[np.absolute(one_body) < 1e-8 ] = 0
two_body = odho.u
two_body[np.absolute(two_body) < 1e-8 ] = 0

# Coupled Cluster
print('Reference energy:',odho.compute_reference_energy())
ccsd = CCSD(odho,verbose=False)
ccsd.compute_ground_state()
print('ECCSD =',ccsd.compute_energy())


# Prepare circuit list
H2 = SecondQuantizedHamiltonian(n,l)
H2.set_integrals(one_body,two_body,anti_symmetric=True)
H2.get_circuit()
circuit_list = H2.to_circuit_list(ptype='vqe')

ansatz = UnitaryCoupledCluster(n,l,'S',1)
og_params = ansatz.new_parameters(H2.h,
                                  H2.v)
print(ansatz.num_S)
print(ansatz.num_D)
print(len(og_params))

grid_points = 500
x = np.linspace(0,2*np.pi,grid_points)
params = [[i] for i in x]
Es = np.zeros(grid_points)

theta = og_params
vqe = VQE(n_qubits = l,
    ansatz = ansatz,
    circuit_list = circuit_list,
    shots = 1000,
    ancilla=0,
    max_energy=False,
    prnt=False)

for i,theta in enumerate(params):
    print('{} / {}'.format(i+1,grid_points))
    Es[i] = vqe.expval(theta)

np.save('data/brute/x.npy',x)
np.save('data/brute/values1000.npy',Es)

