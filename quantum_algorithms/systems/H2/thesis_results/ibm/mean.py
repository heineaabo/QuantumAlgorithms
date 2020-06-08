import numpy as np
import sys


import numpy as np

import sys
sys.path.append('../..')
from matrix import get_h2_matrix
sys.path.append('../../../..')
sys.path.append('../../../../optimizers')
sys.path.append('../../../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian,CircuitList,PauliString,X
from fci import FCI

from vqe import VQE
from optimizers import Minimizer

# IBM PREP
import qiskit as qk
qk.IBMQ.load_account()
provider = qk.IBMQ.get_provider('ibm-q')
#qcomp = provider.get_backend('ibmq_essex')
qcomp = provider.get_backend('ibmq_london')

n = 2
l = 4

R = 2.2
Rstr = '220'

h,v,Enuc,E = get_h2_matrix(R)
h2 = SecondQuantizedHamiltonian(n,l,h,v,nuclear_repulsion=Enuc,add_spin=True)
print('Real FCI:',E['fci'])
print('HF:',E['hf'])

h2.group_paulis(qwc=True,gc=True)


E = np.load('calculated/E_{}.npy'.format(Rstr))

t = np.load('calculated/t_{}.npy'.format(Rstr))


theta = t[-1]



options_q = {'shots':4000,
           #'seed':1,
           'optimization_level':1,
           'backend':qcomp,
           'ibmq':True,
           #'device':'ibmq_london',
           #'meas_fit':True,
           #'layout':[0,3,2,1],  
           'layout':[0,1,2,3],  
           #'layout':[1,0,2,3],  
           #'noise_model':True,
           'basis_gates':True,
           #'coupling_map':True,
           'print':True}
#model = VQE(h2,Minimizer('Cobyla',tol=1e-05),'RYPAIRING',options=options_i)
modelq = VQE(h2,Minimizer('Cobyla',tol=1e-06,max_iter=70),'RYPAIRING',options=options_q)
m,v = modelq.get_mean(theta)
np.save('calculated/mean_{}.npy'.format(Rstr),m)
np.save('calculated/var_{}.npy'.format(Rstr),v)
