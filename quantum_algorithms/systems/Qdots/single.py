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

l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals

omega = 1

h,v,E = get_qdot_matrix(n,l,omega)


print(E,FCI(n,l,h,v))

qdot = SecondQuantizedHamiltonian(n,l,h,v,anti_symmetric=True)
qdot.group_paulis()
options = {'shots':10000,
           'print':True,
           'device':'ibmq_london',
           'noise_model':True,
           'meas_fit':True,
           'coupling_map':True,
           'layout':[1,0,2,3],
           'basis_gates':True}

model = VQE(qdot,Minimizer('Cobyla'),'RYPAIRING',ansatz_depth=1,options=options)

theta = model.optimize()
print(model.get_state_coeff(theta))
print(model.get_mean(theta))

