import numpy as np
from tqdm import tqdm

import sys
sys.path.append('..')
from vqe import VQE
from fci import FCI
from systems import get_h2_matrix,get_qdot_matrix
sys.path.append('../optimizers')
from optimizer import Minimizer,QKSPSA
sys.path.append('../../../QuantumCircuitOptimizer')
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
omegas = [1] #np.arange(0.1,3.0,0.02)
#omegas = omegas[[i for i in range(0,len(omegas),5)]]

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
Es_ideal_UCCr = np.zeros_like(omegas)
VARs_ideal_UCCr = np.zeros_like(omegas)
vqe_coeffs_ideal_UCCr = np.zeros((len(omegas),6)) 
Es_noisy_UCCr = np.zeros_like(omegas)
VARs_noisy_UCCr = np.zeros_like(omegas)
vqe_coeffs_noisy_UCCr = np.zeros((len(omegas),6)) 
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
    ### RY
    # Ideal
    #model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'RYPAIRING',options=option_ideal)
    #theta = model.optimize()
    #Es_ideal_RY[i],VARs_ideal_RY[i] = model.get_mean(theta,N=10000,M=10)
    #vqe_coeffs_ideal_RY[i] = model.get_state_coeff(theta)
    ## Noisy
    #model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'RYPAIRING',options=option_noisy_RY)
    #theta = model.optimize()
    #Es_noisy_RY[i],VARs_noisy_RY[i] = model.get_mean(theta,N=10000,M=10)
    #vqe_coeffs_noisy_RY[i] = model.get_state_coeff(theta)
    
    #### RYRZ
    # Ideal
    model = VQE(h2,Minimizer('Cobyla',tol=1/(1000*shots),disp=False),'RYRZ',ansatz_depth=2,options=option_ideal)
    theta = model.optimize()
    Es_ideal_RYRZ[i],VARs_ideal_RYRZ[i] = model.get_mean(theta,N=10000,M=10)
    vqe_coeffs_ideal_RYRZ[i] = model.get_state_coeff(theta)
    ## Noisy
    model = VQE(h2,Minimizer('Cobyla',tol=1/(1000*shots),disp=False),'RYRZ',ansatz_depth=2,options=option_noisy_RYRZ)
    theta = model.optimize()
    Es_noisy_RYRZ[i],VARs_noisy_RYRZ[i] = model.get_mean(theta,N=10000,M=10)
    vqe_coeffs_noisy_RYRZ[i] = model.get_state_coeff(theta)

    ### UCC
    # Ideal
    model = VQE(h2,Minimizer('Cobyla',tol=1/(1000*shots),disp=False),'UCCSD',options=option_ideal)
    theta = model.optimize()
    Es_ideal_UCC[i],VARs_ideal_UCC[i] = model.get_mean(theta,N=10000,M=10)
    vqe_coeffs_ideal_UCC[i] = model.get_state_coeff(theta)
    # r
    model = VQE(h2,Minimizer('Cobyla',tol=1/(1000*shots),disp=False),'UCCDr',options=option_ideal)
    theta = model.optimize()
    Es_ideal_UCCr[i],VARs_ideal_UCCr[i] = model.get_mean(theta,N=10000,M=10)
    vqe_coeffs_ideal_UCCr[i] = model.get_state_coeff(theta)
    ## Noisy
    model = VQE(h2,Minimizer('Cobyla',tol=1/(1000*shots),disp=False),'UCCSD',options=option_noisy_UCC)
    theta = model.optimize()
    Es_noisy_UCC[i],VARs_noisy_UCC[i] = model.get_mean(theta,N=10000,M=10)
    vqe_coeffs_noisy_UCC[i] = model.get_state_coeff(theta)
    #r
    model = VQE(h2,Minimizer('Cobyla',tol=1/(1000*shots),disp=False),'UCCDr',options=option_noisy_UCC)
    theta = model.optimize()
    Es_noisy_UCCr[i],VARs_noisy_UCCr[i] = model.get_mean(theta,N=10000,M=10)
    vqe_coeffs_noisy_UCCr[i] = model.get_state_coeff(theta)

    i += 1

np.save('new_qdot/omegas.npy',omegas)
np.save('new_qdot/fci.npy',FCIs)
np.save('new_qdot/hf.npy',HFs)
np.save('new_qdot/fci_coeff.npy',fci_coeffs)
## RYRZ
np.save('new_qdot/RYRZ_E_i.npy',Es_ideal_RYRZ)
np.save('new_qdot/RYRZ_var_i.npy',VARs_ideal_RYRZ)
np.save('new_qdot/RYRZ_coeff_i.npy',vqe_coeffs_ideal_RYRZ)
np.save('new_qdot/RYRZ_E_n.npy',Es_noisy_RYRZ)
np.save('new_qdot/RYRZ_var_n.npy',VARs_noisy_RYRZ)
np.save('new_qdot/RYRZ_coeff_n.npy',vqe_coeffs_noisy_RYRZ)
## UCCSD
np.save('new_qdot/UCCSD_E_i.npy',Es_ideal_UCC)
np.save('new_qdot/UCCSD_var_i.npy',VARs_ideal_UCC)
np.save('new_qdot/UCCSD_coeff_i.npy',vqe_coeffs_ideal_UCC)
np.save('new_qdot/UCCSD_E_n.npy',Es_noisy_UCC)
np.save('new_qdot/UCCSD_var_n.npy',VARs_noisy_UCC)
np.save('new_qdot/UCCSD_coeff_n.npy',vqe_coeffs_noisy_UCC)
## UCCDr
np.save('new_qdot/UCCDr_E_i.npy',Es_ideal_UCCr)
np.save('new_qdot/UCCDr_var_i.npy',VARs_ideal_UCCr)
np.save('new_qdot/UCCDr_coeff_i.npy',vqe_coeffs_ideal_UCCr)
np.save('new_qdot/UCCDr_E_n.npy',Es_noisy_UCCr)
np.save('new_qdot/UCCDr_var_n.npy',VARs_noisy_UCCr)
np.save('new_qdot/UCCDr_coeff_n.npy',vqe_coeffs_noisy_UCCr)

