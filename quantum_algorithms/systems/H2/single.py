import numpy as np
from h2getter import get_H2

import sys
sys.path.append('../..')
from vqe import VQE
sys.path.append('../../optimizers')
from optimizer import Minimizer
sys.path.append('../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian


l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals

h,v,Enuc,fci = get_H2(1.0)

print('FCI energy :',fci)


import time
t1 = time.time()
h2 = SecondQuantizedHamiltonian(n,l,h,v,nuclear_repulsion=Enuc,add_spin=True)
t2 = time.time()
print('Time:',t2-t1)

h2.group_paulis(qwc=True,gc=True)


#options = {'count_states':False,'shots':1000,'print':True}
#options = {'shots':10000,'print':True}
options = {'shots':10000,
           #'seed':3,
           #'optimization_level':1,
           #'device':'ibmq_essex',
           #'meas_fit':True,
           #'layout':[0,1,2,3],  
           #'noise_model':True,
           #'basis_gates':True,
           #'coupling_map':True,
           'print':True}
model = VQE(h2,Minimizer('Cobyla'),'RYPAIRING',ansatz_depth=1,options=options)

theta = model.optimize()
print(model.get_state_coeff(theta))
print(model.get_mean(theta))

