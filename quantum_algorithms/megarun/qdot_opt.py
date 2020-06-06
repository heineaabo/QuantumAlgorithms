import numpy as np

import sys
sys.path.append('..')
from vqe import VQE
from fci import FCI
from systems import get_h2_matrix,get_qdot_matrix
sys.path.append('../optimizers')
from optimizer import Minimizer,QKSPSA
sys.path.append('../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian


l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals

shots = 8192
option_ideal = {'shots':shots,'print':False}
option_noisy_RY = {'shots':shots,'print':False,
            'device':'ibmq_london','noise_model':True,'basis_gates':True,
            'coupling_map':True,'layout':[1,0,2,3],'meas_fit':True}
option_noisy_UCC = {'shots':shots,'print':False,
            'device':'ibmq_london','noise_model':True,'basis_gates':True,
            'coupling_map':True,'layout':[0,3,2,1],'meas_fit':True}

omega = 1
h,v,E = get_qdot_matrix(n,l,omega)
print('Energies:',E)
h2 = SecondQuantizedHamiltonian(n,l,h,v,anti_symmetric=True)
h2.group_paulis(qwc=True,gc=True)

ansatze = ['UCCD','UCCDr','UCCSD','UCCSDr']
#ansatze = ['RYPAIRING','UCCD','UCCDr','UCCSD','UCCSDr']
methods = ['cobyla','powell','nelder-mead','SPSA']
for i,ansatz in enumerate(ansatze):
    for j,method in enumerate(methods):
        print('Now: {} with {}'.format(ansatz,method))
        if method == 'SPSA':
            optim = QKSPSA()
        else:
            #optim = Minimizer(method,disp=False,adapt=True)
            optim = Minimizer(method,tol=1/(100*shots),disp=False)
        model = VQE(h2,optim,ansatz,options=option_ideal)
        theta = model.optimize()
        np.save('qdot/{}/E_conv_w{}_{}.npy'.format(method,omega,ansatz),np.asarray(model.energies))
        if method == 'powell' and ansatz in ['RYPAIRING','UCCD','UCCDr']:
            theta = [theta]
        mean,var = model.get_mean(theta)
        print(mean)
        np.save('qdot/{}/E_w{}_final_{}.npy'.format(method,omega,ansatz),np.asarray(mean))


