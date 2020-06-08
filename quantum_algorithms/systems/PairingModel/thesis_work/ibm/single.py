import numpy as np

import sys
sys.path.append('../..')
from matrix import get_pairing_matrix
sys.path.append('../../../..')
sys.path.append('../../../../optimizers')
sys.path.append('../../../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian
from fci import FCI

from vqe import VQE
from optimizers import Minimizer

# IBM PREP
import qiskit as qk
qk.IBMQ.load_account()
provider = qk.IBMQ.get_provider('ibm-q')
qcomp = provider.get_backend('ibmq_essex')
#qcomp = provider.get_backend('ibmq_london')



l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals
delta = 1 # Level spacing
g = 1     # Interaction strength

h,v = get_pairing_matrix(n,l,delta,g)

Efci = FCI(n,l,h,v)
print('FCI energy :',Efci)

import time
t1 = time.time()
pairing =  PairingHamiltonian(n,l,h,v)
t2 = time.time()
print('Time:',t2-t1)

pairing.group_paulis(qwc=True,gc=True)



#theta = [5.829889373194686] # Hardcode good parameter
#theta = [-0.4]


#options = {'count_states':False,'shots':1000,'print':True}
#options = {'shots':10000,'print':True}
options = {'shots':1000,
           #'seed':1,
           'optimization_level':3,
           'backend':qcomp,
           'ibmq':True,
           #'device':'ibmq_london',
           #'meas_fit':True,
           #'layout':[0,3,2,1],  
           #'layout':[1,0,2,3],  
           #'noise_model':True,
           'basis_gates':True,
           #'coupling_map':True,
           'print':True}
model = VQE(pairing,Minimizer('Cobyla',tol=1e-05),'UCCD',options=options)
#print(len(model.circuit_list))
#theta = model.optimize(theta)
#print(theta)
#print(model.get_mean(theta))
print(model.expval(theta))

