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

# FIVE PARAM
# RYRZ
E_RYRZ = np.zeros_like(Rs)
var_RYRZ = np.zeros_like(Rs)
coeff_RYRZ = np.zeros((len(Rs),6)) 

# UCCSD
E_UCCSD = np.zeros_like(Rs)
var_UCCSD = np.zeros_like(Rs)
coeff_UCCSD = np.zeros((len(Rs),6)) 
# r
E_UCCSDr = np.zeros_like(Rs)
var_UCCSDr = np.zeros_like(Rs)
coeff_UCCSDr = np.zeros((len(Rs),6)) 

i = 0
for R in tqdm(Rs):
    R = np.round(R,4)
    h,v,Enuc,E = get_h2_matrix(R)
    h2 = SecondQuantizedHamiltonian(n,l,h,v,nuclear_repulsion=Enuc,add_spin=True)

    h2.group_paulis(qwc=True,gc=True)

    #### RYRZ
    ## Noisy
    model = VQE(h2,Minimizer('Cobyla',tol=1/(10000*shots),disp=False),'RYRZ',ansatz_depth=2,options=option_noisy_RYRZ)
    theta = model.optimize()
    E_RYRZ[i],var_RYRZ[i] = model.get_mean(theta,N=10000,M=10)
    coeff_RYRZ[i] = model.get_state_coeff(theta)

    #### UCC
    ## Noisy
    #model = VQE(h2,Minimizer('Cobyla',tol=1/(1000*shots),disp=False),'UCCSD',options=option_noisy_UCC)
    #theta = model.optimize()
    #E_UCCSD[i],var_UCCSD[i] = model.get_mean(theta,N=10000,M=10)
    #coeff_UCCSD[i] = model.get_state_coeff(theta)
    ## Noisy
    model = VQE(h2,Minimizer('Cobyla',tol=1/(10000*shots),disp=False),'UCCSDr',options=option_noisy_UCC)
    theta = model.optimize()
    E_UCCSDr[i],var_UCCSDr[i] = model.get_mean(theta,N=10000,M=10)
    coeff_UCCSDr[i] = model.get_state_coeff(theta)

    i += 1

## RYRZ
np.save('noisy/five_param/E_RYRZ_s10k.npy',E_RYRZ)
np.save('noisy/five_param/var_RYRZ_s10k.npy',var_RYRZ)
np.save('noisy/five_param/coeff_RYRZ_s10k.npy',coeff_RYRZ)
## UCC
#np.save('noisy/five_param/E_UCCSD.npy',E_UCCSD)
#np.save('noisy/five_param/var_UCCSD.npy',var_UCCSD)
#np.save('noisy/five_param/coeff_UCCSD.npy',coeff_UCCSD)
np.save('noisy/five_param/E_UCCSDr_s10k.npy',E_UCCSDr)
np.save('noisy/five_param/var_UCCSDr_s10k.npy',var_UCCSDr)
np.save('noisy/five_param/coeff_UCCSDr_s10k.npy',coeff_UCCSDr)
