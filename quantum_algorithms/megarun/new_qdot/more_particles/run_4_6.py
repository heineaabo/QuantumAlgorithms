import numpy as np
from tqdm import tqdm

import sys
sys.path.append('../../..')
from vqe import VQE
from fci import FCI
from systems import get_qdot_matrix
sys.path.append('../../../optimizers')
from optimizer import Minimizer,QKSPSA
sys.path.append('../../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian


l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals

shots = 8192
option_ideal = {'shots':shots,'print':False}
option_noisy_RY = {'shots':shots,'print':False,
            'device':'ibmq_london','noise_model':True,'basis_gates':True,
            'coupling_map':True,'layout':[1,0,2,3],'meas_fit':True}
option_noisy_RYRZ = {'shots':shots,'print':False,
            'device':'ibmq_london','noise_model':True,'basis_gates':True,
            'meas_fit':True}
option_noisy_UCC = {'shots':shots,'print':False,
            'device':'ibmq_london','noise_model':True,'basis_gates':True,
            'coupling_map':True,'layout':[0,3,2,1],'meas_fit':True}


######## QDOT
omegas = np.arange(0.1,3.0,0.02)
omegas = [] #omegas[[i for i in range(0,len(omegas),5)]]

FCIs = np.zeros_like(omegas)
HFs = np.zeros_like(omegas)
fci_coeffs = np.zeros((len(omegas),6)) 
#
#Es_ideal_RY = np.zeros_like(omegas)
#VARs_ideal_RY = np.zeros_like(omegas)
#vqe_coeffs_ideal_RY = np.zeros((len(omegas),6)) 
#Es_noisy_RY = np.zeros_like(omegas)
#VARs_noisy_RY = np.zeros_like(omegas)
#vqe_coeffs_noisy_RY = np.zeros((len(omegas),6)) 

# RYRZ
Es_ideal_RYRZ = np.zeros_like(omegas)
VARs_ideal_RYRZ = np.zeros_like(omegas)
vqe_coeffs_ideal_RYRZ = np.zeros((len(omegas),6)) 
Es_noisy_RYRZ = np.zeros_like(omegas)
VARs_noisy_RYRZ = np.zeros_like(omegas)
vqe_coeffs_noisy_RYRZ = np.zeros((len(omegas),6)) 

# UCC
#Es_ideal_UCCr = np.zeros_like(omegas)
#VARs_ideal_UCCr = np.zeros_like(omegas)
#vqe_coeffs_ideal_UCCr = np.zeros((len(omegas),6)) 
#Es_noisy_UCCr = np.zeros_like(omegas)
#VARs_noisy_UCCr = np.zeros_like(omegas)
#vqe_coeffs_noisy_UCCr = np.zeros((len(omegas),6)) 
Es_ideal_UCC = np.zeros_like(omegas)
VARs_ideal_UCC = np.zeros_like(omegas)
vqe_coeffs_ideal_UCC = np.zeros((len(omegas),6)) 
Es_noisy_UCC = np.zeros_like(omegas)
VARs_noisy_UCC = np.zeros_like(omegas)
vqe_coeffs_noisy_UCC = np.zeros((len(omegas),6)) 

i = 0
for omega in tqdm(omegas):
    h,v,E = get_qdot_matrix(n,l,omega)
    ev,vv = FCI(n,l,h,v,ret_all=True)
    FCIs[i] = ev[0]
    HFs[i] = E[0].real
    fci_coeffs[i] = vv[:,0].real
    h2 = SecondQuantizedHamiltonian(n,l,h,v,anti_symmetric=True)

    h2.group_paulis(qwc=True,gc=True)

    ### UCC
    # Ideal
    model = VQE(h2,Minimizer('Cobyla',tol=1/(1000*shots),disp=False),'UCCSDr',options=option_ideal)
    theta = model.optimize()
    Es_ideal_UCC[i],VARs_ideal_UCC[i] = model.get_mean(theta,N=10000,M=10)
    vqe_coeffs_ideal_UCC[i] = model.get_state_coeff(theta)

    i += 1

np.save('omegas.npy',omegas)
np.save('fci.npy',FCIs)
np.save('hf.npy',HFs)
np.save('fci_coeff.npy',fci_coeffs)
## UCCSDr
np.save('UCCSDr_E_i.npy',Es_ideal_UCC)
np.save('UCCSDr_var_i.npy',VARs_ideal_UCC)
np.save('UCCSDr_coeff_i.npy',vqe_coeffs_ideal_UCC)

