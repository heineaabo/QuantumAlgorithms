import numpy as np
from tqdm import tqdm

import sys
sys.path.append('../..')
sys.path.append('../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian

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
theta[0] = 5.829889373194686 # Hardcode good parameter

Pairing = SecondQuantizedHamiltonian(n,l)
Pairing.set_integrals(h_pq,h_pqrs)
Pairing.get_circuit()
circuit_list = Pairing.to_circuit_list(ptype='vqe')

all_shots = np.linspace(1000,100000,100)
Es = np.zeros(100)
legal = np.zeros(100)
illegal = np.zeros(100)

i = 0
for shot in tqdm(all_shots):
    model = VQE(n_qubits = l,
        ansatz = UCCD,
        circuit_list = circuit_list,
        shots = 500, #int(shot),
        ancilla=0,
        seed = None,
        prnt=False,
        max_energy=False,
        count_states=True)
    Es[i] = model.expval(theta)
    legal[i] = model.legal
    illegal[i] = model.illegal
    i += 1

np.save('data/shots/shots.npy',all_shots)
np.save('data/shots/values.npy',Es)
np.save('data/shots/legal.npy',legal)
np.save('data/shots/illegal.npy',illegal)
