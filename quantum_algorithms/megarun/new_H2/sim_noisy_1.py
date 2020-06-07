import numpy as np
from tqdm import tqdm

import sys
sys.path.append('../..')
from vqe import VQE
from fci import FCI
from systems import get_h2_matrix,get_qdot_matrix
sys.path.append('../../optimizers')
from optimizer import Minimizer,QKSPSA
sys.path.append('../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian


l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals

shots = 8192
option_noisy_RY = {'shots':shots,'print':False,
            'device':'ibmq_london','noise_model':True,'basis_gates':True,
            'coupling_map':True,'layout':[1,0,2,3],'meas_fit':True}
option_noisy_RYRZ = {'shots':shots,'print':False,
            'device':'ibmq_london','noise_model':True,'basis_gates':True,
            'meas_fit':True}
option_noisy_UCC = {'shots':shots,'print':False,
            'device':'ibmq_london','noise_model':True,'basis_gates':True,
            'coupling_map':True,'layout':[0,3,2,1],'meas_fit':True}


########### H2

#Rs = np.arange(0.3,3.0,0.02)
Rs = np.arange(0.3,3.0,0.02)
Rs = Rs[[i for i in range(0,len(Rs),5)]]

FCIs = np.zeros_like(Rs)
HFs = np.zeros_like(Rs)
CCs = np.zeros_like(Rs)
fci_coeffs = np.zeros((len(Rs),6)) 

# RY
E_RY = np.zeros_like(Rs)
var_RY = np.zeros_like(Rs)
coeff_RY = np.zeros((len(Rs),6)) 
# UCC
E_UCCD = np.zeros_like(Rs)
var_UCCD = np.zeros_like(Rs)
coeff_UCCD = np.zeros((len(Rs),6)) 
E_UCCDr = np.zeros_like(Rs)
var_UCCDr = np.zeros_like(Rs)
coeff_UCCDr = np.zeros((len(Rs),6)) 

i = 0
for R in tqdm(Rs):
    R = np.round(R,4)
    h,v,Enuc,E = get_h2_matrix(R)
    h2 = SecondQuantizedHamiltonian(n,l,h,v,nuclear_repulsion=Enuc,add_spin=True)
    ev,vv = FCI(n,l,h2.h,4*h2.v,ret_all=True)
    fci_coeffs[i] = vv[:,0].real
    FCIs[i] = E['fci']
    #FCIs[i] = ev[0].real
    HFs[i] = E['hf']
    CCs[i] = E['ccsd']


    h2.group_paulis(qwc=True,gc=True)

    #### RY
    # Noisy
    model = VQE(h2,Minimizer('Cobyla',tol=1/(1000*shots),disp=False),'RYPAIRING',options=option_noisy_RY)
    theta = model.optimize()
    E_RY[i],var_RY[i] = model.get_mean(theta,N=10000,M=10)
    coeff_RY[i] = model.get_state_coeff(theta)

    #### UCC
    ## Noisy
    model = VQE(h2,Minimizer('Cobyla',tol=1/(1000*shots),disp=False),'UCCD',options=option_noisy_UCC)
    theta = model.optimize()
    E_UCCD[i],var_UCCD[i] = model.get_mean(theta,N=10000,M=10)
    coeff_UCCD[i] = model.get_state_coeff(theta)
    ## Noisy
    model = VQE(h2,Minimizer('Cobyla',tol=1/(1000*shots),disp=False),'UCCDr',options=option_noisy_UCC)
    theta = model.optimize()
    E_UCCDr[i],var_UCCDr[i] = model.get_mean(theta,N=10000,M=10)
    coeff_UCCDr[i] = model.get_state_coeff(theta)

    i += 1

np.save('bonds.npy',Rs)
np.save('fci.npy',FCIs)
np.save('hf.npy',HFs)
np.save('cc.npy',CCs)
np.save('fci_coeff.npy',fci_coeffs)
## RY
np.save('noisy/single_param/E_RY.npy',E_RY)
np.save('noisy/single_param/var_RY.npy',var_RY)
np.save('noisy/single_param/coeff_RY.npy',coeff_RY)
## UCC
np.save('noisy/single_param/E_UCCD.npy',E_UCCD)
np.save('noisy/single_param/var_UCCD.npy',var_UCCD)
np.save('noisy/single_param/coeff_UCCD.npy',coeff_UCCD)
# r
np.save('noisy/single_param/E_UCCDr.npy',E_UCCDr)
np.save('noisy/single_param/var_UCCDr.npy',var_UCCDr)
np.save('noisy/single_param/coeff_UCCDr.npy',coeff_UCCDr)
