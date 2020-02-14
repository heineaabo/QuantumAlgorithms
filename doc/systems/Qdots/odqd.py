import numpy as np
import sys
sys.path.append('../..')
from vqe import *
import numpy as np
from ucc import UnitaryCoupledCluster

from quantum_circuit import QuantumCircuit,Hamiltonian
from quantum_circuit.utils import molecular2sec_quant

from quantum_systems import ODQD
from quantum_systems.quantum_dots.one_dim.one_dim_potentials import HOPotential
from coupled_cluster import CCSD

import matplotlib.pyplot as plt

l = 4
n = 2


# Generate system
omega = 0.25

grid_length = 5
num_grid_points = 1001

odho = ODQD(n,l,grid_length,num_grid_points)
odho.setup_system(potential=HOPotential(omega),add_spin=False)

one_body = odho.h
one_body[np.absolute(one_body) < 1e-8 ] = 0
two_body = odho.u
two_body[np.absolute(two_body) < 1e-8 ] = 0
#print(np.real(one_body))
#print(np.real(two_body))


# Coupled Cluster
#print('Reference energy:',odho.compute_reference_energy())
#ccsd = CCSD(odho,verbose=False)
#ccsd.compute_ground_state()
#print('ECCSD =',ccsd.compute_energy())

one_body,two_body = molecular2sec_quant(one_body,two_body)
#print(one_body)
#print(two_body)


# Variational Quantum Eigensolver
uccsd = UnitaryCoupledCluster(n,l,'SD',1)
theta = uccsd.new_parameters()

odho_vqe = Hamiltonian(n,l)
odho_vqe.set_integrals(one_body,two_body)
odho_vqe.get_circuit()
circuit_list = odho_vqe.circuit_list
print(odho_vqe)

vqe = VQE(n_qubits = l,
        ansatz = uccsd,
        circuit_list = circuit_list,
        shots = 500,
        ancilla=0,
        prnt=True,
        max_energy=False)

t = vqe.optimize_classical(theta,method='Powell')


