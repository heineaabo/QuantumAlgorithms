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

R = 0.74
h,v,Enuc,E = get_h2_matrix(R)
h2 = SecondQuantizedHamiltonian(n,l,h,v,nuclear_repulsion=Enuc,add_spin=True)
ev,vv = FCI(n,l,h2.h,h2.v,ret_all=True)
fci_coeffs[i] = vv[:,0].real
FCIs[i] = E['fci']
FCIs[i] = ev[0].real
HFs[i] = E['hf']
CCs[i] = E['ccsd']


h2.group_paulis(qwc=True,gc=True)

ansatze = ['RYPAIRING','UCCD','UCCDr','UCCSD','UCCSDr']
methods = ['Cobyla','Powell','Nelder-Mead','SPSA']
for i,ansatz in enumerate(ansatze):
    for j,method in enumerate(methods):
        if method == 'SPSA':
            optim = QKSPSA()
        else:
            optim = Minimizer(method,tol=1/(100*shots),disp=False)
    model = VQE(h2,optim,ansatz,options=option_ideal)
    np.save('{}/E_conv_R{}_{}.npy'.format(method,074,ansatz),np.asarray(model.energies))


