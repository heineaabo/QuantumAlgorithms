import numpy as np


import sys
sys.path.append('../..')
sys.path.append('../../rest')
from Qoperator2 import Hamiltonian as hamil
sys.path.append('../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian

from vqe import VQE
from ucc import UnitaryCoupledCluster
from spsa import SPSA
from quantum_gradient_descent import QuantumGradientDescent

from qiskit.aqua.components.optimizers import SPSA as qkspsa


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
new_theta = UCCD.new_parameters(h_pq,h_pqrs)

Pairing = SecondQuantizedHamiltonian(n,l)
Pairing.set_integrals(h_pq,h_pqrs)
Pairing.get_circuit()
circuit_list = Pairing.to_circuit_list(ptype='vqe')
#print('new')
for i in circuit_list: print(i)

#p2 = hamil(l)
#cir = p2.get_circuits(h_pq,h_pqrs)
#print('old')
#for i in cir: print(i)


model = VQE(n_qubits = l,
    ansatz = UCCD,
    circuit_list = circuit_list,
    shots = 500,
    ancilla=0,
    prnt=True,
    max_energy=False)
#model.optimize_classical(new_theta,method='Cobyla')
#model.optimize_gradient(new_theta)
#options = {'feedback':1,'grad_avg':5}
#optimization = SPSA(model.expval,
#                    new_theta,
#                    min_change=0.1,
#                    noise_var = 0.01,
#                    options=options)
#optimization()
#method = 'SPSA'
#np.save('data/'+method+'.npy', model.energies)

#optimize = QuantumGradientDescent(model.expval,new_theta)
#optimize()

opt = qkspsa()
opt.optimize(len(new_theta),model.expval,initial_point=new_theta)

