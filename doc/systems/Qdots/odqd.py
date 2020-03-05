import numpy as np
import sys
sys.path.append('../..')
sys.path.append('../../../../QuantumCircuitOptimizer')
from vqe import *
import numpy as np
from ucc import UnitaryCoupledCluster

from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian
from quantum_circuit.utils import molecular2sec_quant
import matplotlib.pyplot as plt

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

#one_body,two_body = molecular2sec_quant(one_body,two_body)
#print(one_body)
#print(two_body)
#for i in range(two_body.shape[0]):
#    for j in range(two_body.shape[1]):
#        for k in range(two_body.shape[2]):
#            for l in range(two_body.shape[3]):
#                if not np.isclose(two_body[i,j,k,l],0):
#                    print(i,j,k,l,two_body[i,j,k,l])


# Variational Quantum Eigensolver
uccsd = UnitaryCoupledCluster(n,l,'SD',1)
theta = uccsd.new_parameters()

odho_vqe = SecondQuantizedHamiltonian(n,l)
odho_vqe.set_integrals(one_body,two_body,anti_symmetric=True)
odho_vqe.get_circuit()
circuit_list = odho_vqe.to_circuit_list(ptype='vqe')

vqe = VQE(n_qubits = l,
        ansatz = uccsd,
        circuit_list = circuit_list,
        shots = 500,
        ancilla=0,
        prnt=True,
        max_energy=False)
print('Starting minimization')
t = vqe.optimize_classical(theta,method='Powell')


