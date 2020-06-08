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

np.save('ibmq_energies.npy',np.array([]))
np.save('ibmq_theta.npy',np.array([]))

l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals
delta = 1 # Level spacing
g = 1     # Interaction strength

h,v,Enuc,E = get_h2_matrix(R)
h2 = SecondQuantizedHamiltonian(n,l,h,v,nuclear_repulsion=Enuc,add_spin=True)
print('Real FCI:',E['fci']
print('HF:',E['hf']

Efci = FCI(n,l,h,v)
print('FCI energy :',Efci+Enuc,Efci-Enuc)

h2.group_paulis(qwc=True,gc=True)
#options = {'count_states':False,'shots':1000,'print':True}
options_i = {'shots':2000,'optimization_level':1,
           #'noise_model':True,
           #'meas_fit':True,
           # 'device':'ibmq_london',
           'print':True}
options_q = {'shots':2000,
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
modelq = VQE(h2,Minimizer('Cobyla',tol=1e-05,max_iter=70),'RYPAIRING',options=options_q)
theta = modelq.optimize()
print(theta)
#print(model.get_mean(theta))
#print('simulation')
#print(model.expval(theta))
#print('Real quantum')
#print(modelq.expval(theta))

