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




device = ''
if len(sys.argv) > 1:
    device = sys.argv[1]
else:
    device = 'london'

ansatz = 'RYPAIRING'


gs = np.arange(-2,2,0.05)
data = np.zeros((len(gs),3))
fci  = np.zeros_like(gs)
shots = 1000
i = 0
for g in tqdm(gs):
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
    for o,options in enumerate(optionz):
        model = VQE(pairing,Minimizer('Cobyla',tol=1e-04,disp=False),ansatz,options=options)
        theta = model.optimize()
        data[i,o],_ = model.get_mean(theta)
    fci[i] = FCI(n,l,h_pq,h_pqrs) 
    i += 1
np.save('{}_{}_{}shots.npy'.format(device,ansatz,shots),data)
np.save('fci.npy',fci)
