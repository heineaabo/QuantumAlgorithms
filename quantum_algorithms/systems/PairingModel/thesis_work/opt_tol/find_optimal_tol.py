import numpy as np

import sys
sys.path.append('../../../..')
from vqe import VQE
sys.path.append('../../../../optimizers')
from optimizer import Minimizer
sys.path.append('../../../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian


l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals
delta = 1 # Level spacing
g = 1     # Interaction strength


# Matrix elements
h_pq = np.identity(l)
for p in range(l):
	h_pq[p,p] *= delta*(p - (p%2))/2
	
h_pqrs = np.zeros((l,l,l,l))
for p in range(0,l-1,2):
	for r in range(0,l-1,2):
		h_pqrs[p,p+1,r,r+1] = -0.5*g



pairing =  PairingHamiltonian(n,l,h_pq,h_pqrs)


theta = [5.829889373194686] # Hardcode good parameter

device = ''
if len(sys.argv) > 1:
    device = sys.argv[1]
else:
    device = 'london'

data = np.zeros((2,3,25))

shotz = [1000,10000]
for s,shots in enumerate(shotz):
    print('For shots:',shots)
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
               'meas_fitter':False,
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
        print('For option',o+1)
        for i in range(25):
            model = VQE(pairing,Minimizer('Powell'),'RYPAIRING',options=options)
            data[s,o,i] = model.expval(theta)
np.save('variance_{}_{}t.npy'.format(device,25),data)
