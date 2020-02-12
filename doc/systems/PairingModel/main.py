import numpy as np

from quantum_circuit import QuantumCircuit,Hamiltonian

import sys
sys.path.append('../..')

from vqe import VQE
#from Qoperator import *
from ucc import UnitaryCoupledCluster
from spsa import SPSA


l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals
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

UCCD = UnitaryCoupledCluster(n,l,'D',1)
theta = UCCD.new_parameters()

Pairing = Hamiltonian(n,l)
Pairing.set_integrals(h_pq,h_pqrs)
Pairing.get_circuit()
circuit_list = Pairing.to_circuit_list(ptype='vqe')

model = VQE(n_qubits = l,
    ansatz = UCCD,
    circuit_list = circuit_list,
    shots = 500,
    ancilla=0,
    max_energy=False)
#model.optimize_classical(theta,method='Powell')
optimization = SPSA(model.expval,
                    theta)
optimization.run()

