import numpy as np
from tqdm import tqdm

import sys
sys.path.append('../..')
from matrix import get_pairing_matrix
sys.path.append('../../../..')
from vqe import VQE
from fci import FCI
sys.path.append('../../../../optimizers')
from optimizer import Minimizer
sys.path.append('../../../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian


l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals
delta = 1 # Level spacing
g = 1



device = ''
if len(sys.argv) > 1:
    device = sys.argv[1]
else:
    device = 'london'

ansatz = 'RYPAIRING'


shots = 1000
h_pq,h_pqrs = get_pairing_matrix(n,l,delta,g)
pairing =  PairingHamiltonian(n,l,h_pq,h_pqrs)
pairing.group_paulis()
option1 = {'shots':shots,
           'optimization_level':1,
           #'seed':1,
           'print':False}
option2 = {'shots':shots,
           'optimization_level':1,
           'device':'ibmq_{}'.format(device),
           'layout':[1,0,2,3],
           'noise_model':True,
           'basis_gates':True,
           'coupling_map':True,
           #'seed':1,
           'meas_fit':False,
           'print':False}
option3 = {'shots':shots,
           'optimization_level':1,
           'device':'ibmq_{}'.format(device),
           'layout':[1,0,2,3],
           'noise_model':True,
           'basis_gates':True,
           'coupling_map':True,
           #'seed':1,
           'print':False}
optionz = [option1,option2,option3]
model = VQE(pairing,Minimizer('Cobyla',tol=1e-05,disp=False),ansatz,options=optionz[0])
theta = model.optimize()
mean,var = model.get_mean(theta)
coeff = model.get_state_coeff(theta)
fciE,fciV = FCI(n,l,h_pq,h_pqrs,ret_all=True) 
print(fciE[0])
print(mean,var)

fciCoef = np.power(fciV[:,0],2)
print(fciCoef)
print(coeff)
print(np.linalg.norm(fciCoef-coeff))
