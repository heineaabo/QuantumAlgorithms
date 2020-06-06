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


########### H2

#Rs = np.arange(0.3,3.0,0.02)
Rs = np.arange(0.3,3.0,0.02)
Rs = Rs[[i for i in range(0,len(Rs),5)]]

FCIs = np.zeros_like(Rs)
HFs = np.zeros_like(Rs)
CCs = np.zeros_like(Rs)
fci_coeffs = np.zeros((len(Rs),6)) 

# RY
#Es_ideal_RY = np.zeros_like(Rs)
#VARs_ideal_RY = np.zeros_like(Rs)
#vqe_coeffs_ideal_RY = np.zeros((len(Rs),6)) 
#Es_noisy_RY = np.zeros_like(Rs)
#VARs_noisy_RY = np.zeros_like(Rs)
#vqe_coeffs_noisy_RY = np.zeros((len(Rs),6)) 

# RYRZ
Es_ideal_RYRZ = np.zeros_like(Rs)
VARs_ideal_RYRZ = np.zeros_like(Rs)
vqe_coeffs_ideal_RYRZ = np.zeros((len(Rs),6)) 
#Es_noisy_RYRZ = np.zeros_like(Rs)
#VARs_noisy_RYRZ = np.zeros_like(Rs)
#vqe_coeffs_noisy_RYRZ = np.zeros((len(Rs),6)) 

# UCC
Es_ideal_UCCr = np.zeros_like(Rs)
VARs_ideal_UCCr = np.zeros_like(Rs)
vqe_coeffs_ideal_UCCr = np.zeros((len(Rs),6)) 
Es_noisy_UCCr = np.zeros_like(Rs)
VARs_noisy_UCCr = np.zeros_like(Rs)
vqe_coeffs_noisy_UCCr = np.zeros((len(Rs),6)) 
#Es_ideal_UCC = np.zeros_like(Rs)
#VARs_ideal_UCC = np.zeros_like(Rs)
#vqe_coeffs_ideal_UCC = np.zeros((len(Rs),6)) 
#Es_noisy_UCC = np.zeros_like(Rs)
#VARs_noisy_UCC = np.zeros_like(Rs)
#vqe_coeffs_noisy_UCC = np.zeros((len(Rs),6)) 

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
    #model = VQE(h2,Minimizer('Cobyla',tol=1/(1000*shots),disp=False),'RYRZ',ansatz_depth=2,options=option_ideal)
    ###model = VQE(h2,QKSPSA(max_iter=100),'RYRZ',ansatz_depth=2,options=option_ideal)
    #theta = model.optimize()
    #Es_ideal_RYRZ[i],VARs_ideal_RYRZ[i] = model.get_mean(theta,N=10000,M=10)
    #vqe_coeffs_ideal_RYRZ[i] = model.get_state_coeff(theta)
    #print('RYRZ number of parameters:',len(theta))
    #print('RYRZ number of evauliations:',model.evals)
    #print('Energy:',Es_ideal_RYRZ[i])
    ## Noisy
    #model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'RYRZ',ansatz_depth=2,options=option_noisy_RYRZ)
    #theta = model.optimize()
    #Es_noisy_RYRZ[i],VARs_noisy_RYRZ[i] = model.get_mean(theta,N=10000,M=10)
    #vqe_coeffs_noisy_RYRZ[i] = model.get_state_coeff(theta)

    #### UCC
    # Ideal
    #model = VQE(h2,Minimizer('Cobyla',tol=1/(1000*shots),disp=False),'UCCSD',options=option_ideal)
    ###model = VQE(h2,QKSPSA(max_iter=100),'UCCSD',options=option_ideal)
    #theta = model.optimize()
    #Es_ideal_UCC[i],VARs_ideal_UCC[i] = model.get_mean(theta,N=10000,M=10)
    #vqe_coeffs_ideal_UCC[i] = model.get_state_coeff(theta)
    # r
    model = VQE(h2,Minimizer('Cobyla',tol=1/(1000*shots),disp=False),'UCCSDr',options=option_ideal)
    ###model = VQE(h2,QKSPSA(max_iter=100),'UCCSDr',options=option_ideal)
    theta = model.optimize()
    Es_ideal_UCCr[i],VARs_ideal_UCCr[i] = model.get_mean(theta,N=10000,M=10)
    vqe_coeffs_ideal_UCCr[i] = model.get_state_coeff(theta)
    ## Noisy
    model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'UCCSDr',options=option_noisy_UCC)
    theta = model.optimize()
    Es_noisy_UCCr[i],VARs_noisy_UCCr[i] = model.get_mean(theta,N=10000,M=10)
    vqe_coeffs_noisy_UCCr[i] = model.get_state_coeff(theta)

    i += 1

##np.save('H2/bonds.npy',Rs)
##np.save('H2/fci.npy',FCIs)
##np.save('H2/hf.npy',HFs)
##np.save('H2/cc.npy',CCs)
##np.save('H2/fci_coeff.npy',fci_coeffs)
## RY
##np.save('H2/cobyla/E_ideal_RY.npy',Es_ideal_RY)
##np.save('H2/cobyla/var_ideal_RY.npy',VARs_ideal_RY)
##np.save('H2/cobyla/coeff_ideal_RY.npy',vqe_coeffs_ideal_RY)
##np.save('H2/cobyla/E_noisy_RY.npy',Es_noisy_RY)
##np.save('H2/cobyla/var_noisy_RY.npy',VARs_noisy_RY)
##np.save('H2/cobyla/coeff_noisy_RY.npy',vqe_coeffs_noisy_RY)
## RYRZ
#np.save('H2/five_param/E_ideal_RYRZ.npy',Es_ideal_RYRZ)
#np.save('H2/five_param/var_ideal_RYRZ.npy',VARs_ideal_RYRZ)
#np.save('H2/five_param/coeff_ideal_RYRZ.npy',vqe_coeffs_ideal_RYRZ)
##np.save('qdot/cobyla/E_noisy_RYRZ.npy',Es_noisy_RYRZ)
##np.save('qdot/cobyla/var_noisy_RYRZ.npy',VARs_noisy_RYRZ)
##np.save('qdot/cobyla/coeff_noisy_RYRZ.npy',vqe_coeffs_noisy_RYRZ)
## UCC
##np.save('H2/five_param/E_ideal_UCCSD.npy',Es_ideal_UCC)
##np.save('H2/five_param/var_ideal_UCCSD.npy',VARs_ideal_UCC)
##np.save('H2/five_param/coeff_ideal_UCCSD.npy',vqe_coeffs_ideal_UCC)
np.save('H2/five_param/E_ideal_UCCSDr.npy',Es_ideal_UCCr)
np.save('H2/five_param/var_ideal_UCCSDr.npy',VARs_ideal_UCCr)
np.save('H2/five_param/coeff_ideal_UCCSDr.npy',vqe_coeffs_ideal_UCCr)
np.save('H2/five_param/E_noisy_UCCSDr.npy',Es_noisy_UCCr)
np.save('H2/five_param/var_noisy_UCCSDr.npy',VARs_noisy_UCCr)
np.save('H2/five_param/coeff_noisy_UCCSDr.npy',vqe_coeffs_noisy_UCCr)
##np.save('H2/cobyla/E_ideal_UCCSD_spsa.npy',Es_ideal_UCC)
##np.save('H2/cobyla/var_ideal_UCCSD_spsa.npy',VARs_ideal_UCC)
##np.save('H2/cobyla/coeff_ideal_UCCSD_spsa.npy',vqe_coeffs_ideal_UCC)
##np.save('H2/cobyla/E_noisy_UCCSD_full.npy',Es_noisy_UCC)
##np.save('H2/cobyla/var_noisy_UCCSD_full.npy',VARs_noisy_UCC)
##np.save('H2/cobyla/coeff_noisy_UCCSD_full.npy',vqe_coeffs_noisy_UCC)
    

######## QDOT
#omegas = np.arange(0.1,3.0,0.02)
#
#FCIs = np.zeros_like(omegas)
#HFs = np.zeros_like(omegas)
#fci_coeffs = np.zeros((len(omegas),6)) 
##
##Es_ideal_RY = np.zeros_like(omegas)
##VARs_ideal_RY = np.zeros_like(omegas)
##vqe_coeffs_ideal_RY = np.zeros((len(omegas),6)) 
##Es_noisy_RY = np.zeros_like(omegas)
##VARs_noisy_RY = np.zeros_like(omegas)
##vqe_coeffs_noisy_RY = np.zeros((len(omegas),6)) 
#
## RYRZ
#Es_ideal_RYRZ = np.zeros_like(omegas)
#VARs_ideal_RYRZ = np.zeros_like(omegas)
#vqe_coeffs_ideal_RYRZ = np.zeros((len(omegas),6)) 
##Es_noisy_RYRZ = np.zeros_like(omegas)
##VARs_noisy_RYRZ = np.zeros_like(omegas)
##vqe_coeffs_noisy_RYRZ = np.zeros((len(omegas),6)) 
#
## UCC
#Es_ideal_UCCr = np.zeros_like(omegas)
#VARs_ideal_UCCr = np.zeros_like(omegas)
#vqe_coeffs_ideal_UCCr = np.zeros((len(omegas),6)) 
#Es_ideal_UCC = np.zeros_like(omegas)
#VARs_ideal_UCC = np.zeros_like(omegas)
#vqe_coeffs_ideal_UCC = np.zeros((len(omegas),6)) 
##Es_noisy_UCC = np.zeros_like(omegas)
##VARs_noisy_UCC = np.zeros_like(omegas)
##vqe_coeffs_noisy_UCC = np.zeros((len(omegas),6)) 
#
#i = 0
#for omega in tqdm(omegas):
#    h,v,E = get_qdot_matrix(n,l,omega)
#    ev,vv = FCI(n,l,h,v,ret_all=True)
#    FCIs[i] = ev[0]
#    HFs[i] = E[0].real
#    fci_coeffs[i] = vv[:,0].real
#    h2 = SecondQuantizedHamiltonian(n,l,h,v,anti_symmetric=True)
#
#    h2.group_paulis(qwc=True,gc=True)
#    ### RY
#    # Ideal
#    #model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'RYPAIRING',options=option_ideal)
#    #theta = model.optimize()
#    #Es_ideal_RY[i],VARs_ideal_RY[i] = model.get_mean(theta,N=10000,M=10)
#    #vqe_coeffs_ideal_RY[i] = model.get_state_coeff(theta)
#    ## Noisy
#    #model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'RYPAIRING',options=option_noisy_RY)
#    #theta = model.optimize()
#    #Es_noisy_RY[i],VARs_noisy_RY[i] = model.get_mean(theta,N=10000,M=10)
#    #vqe_coeffs_noisy_RY[i] = model.get_state_coeff(theta)
#    
#    #### RYRZ
#    # Ideal
#    #model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'RYRZ',ansatz_depth=2,options=option_ideal)
#    model = VQE(h2,QKSPSA(max_iter=100),'RYRZ',ansatz_depth=2,options=option_ideal)
#    theta = model.optimize()
#    Es_ideal_RYRZ[i],VARs_ideal_RYRZ[i] = model.get_mean(theta,N=10000,M=10)
#    vqe_coeffs_ideal_RYRZ[i] = model.get_state_coeff(theta)
#    ## Noisy
#    #model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'RYRZ',ansatz_depth=2,options=option_noisy_RYRZ)
#    #theta = model.optimize()
#    #Es_noisy_RYRZ[i],VARs_noisy_RYRZ[i] = model.get_mean(theta,N=10000,M=10)
#    #vqe_coeffs_noisy_RYRZ[i] = model.get_state_coeff(theta)
#
#    ### UCC
#    # Ideal
#    #model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'UCCSD',options=option_ideal)
#    model = VQE(h2,QKSPSA(max_iter=100),'UCCSD',options=option_ideal)
#    theta = model.optimize()
#    Es_ideal_UCC[i],VARs_ideal_UCC[i] = model.get_mean(theta,N=10000,M=10)
#    vqe_coeffs_ideal_UCC[i] = model.get_state_coeff(theta)
#    # r
#    model = VQE(h2,QKSPSA(max_iter=100),'UCCSDr',options=option_ideal)
#    theta = model.optimize()
#    Es_ideal_UCCr[i],VARs_ideal_UCCr[i] = model.get_mean(theta,N=10000,M=10)
#    vqe_coeffs_ideal_UCCr[i] = model.get_state_coeff(theta)
#    ## Noisy
#    #model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'UCCSD',options=option_noisy_UCC)
#    #theta = model.optimize()
#    #Es_noisy_UCC[i],VARs_noisy_UCC[i] = model.get_mean(theta,N=10000,M=10)
#    #vqe_coeffs_noisy_UCC[i] = model.get_state_coeff(theta)
#
#    i += 1
#
###np.save('qdot/omegas.npy',omegas)
###np.save('qdot/fci.npy',FCIs)
###np.save('qdot/hf.npy',HFs)
###np.save('qdot/fci_coeff.npy',fci_coeffs)
### RY
###np.save('qdot/cobyla/E_ideal_RY.npy',Es_ideal_RY)
###np.save('qdot/cobyla/var_ideal_RY.npy',VARs_ideal_RY)
###np.save('qdot/cobyla/coeff_ideal_RY.npy',vqe_coeffs_ideal_RY)
###np.save('qdot/cobyla/E_noisy_RY.npy',Es_noisy_RY)
###np.save('qdot/cobyla/var_noisy_RY.npy',VARs_noisy_RY)
###np.save('qdot/cobyla/coeff_noisy_RY.npy',vqe_coeffs_noisy_RY)
### RYRZ
#np.save('qdot/SPSA/E_ideal_RYRZ.npy',Es_ideal_RYRZ)
#np.save('qdot/SPSA/var_ideal_RYRZ.npy',VARs_ideal_RYRZ)
#np.save('qdot/SPSA/coeff_ideal_RYRZ.npy',vqe_coeffs_ideal_RYRZ)
###np.save('qdot/cobyla/E_noisy_RYRZ.npy',Es_noisy_RYRZ)
###np.save('qdot/cobyla/var_noisy_RYRZ.npy',VARs_noisy_RYRZ)
###np.save('qdot/cobyla/coeff_noisy_RYRZ.npy',vqe_coeffs_noisy_RYRZ)
### UCC
#np.save('qdot/SPSA/E_ideal_UCCSD.npy',Es_ideal_UCC)
#np.save('qdot/SPSA/var_ideal_UCCSD.npy',VARs_ideal_UCC)
#np.save('qdot/SPSA/coeff_ideal_UCCSD.npy',vqe_coeffs_ideal_UCC)
#np.save('qdot/SPSA/E_ideal_UCCSDr.npy',Es_ideal_UCCr)
#np.save('qdot/SPSA/var_ideal_UCCSDr.npy',VARs_ideal_UCCr)
#np.save('qdot/SPSA/coeff_ideal_UCCSDr.npy',vqe_coeffs_ideal_UCCr)
##np.save('qdot/cobyla/E_noisy_UCCSD_spsa.npy',Es_noisy_UCC)
##np.save('qdot/cobyla/var_noisy_UCCSD_spsa.npy',VARs_noisy_UCC)
##np.save('qdot/cobyla/coeff_noisy_UCCSD_spsa.npy',vqe_coeffs_noisy_UCC)
##
###np.save('qdot/cobyla/E_ideal_UCC.npy',Es_ideal_UCC)
###np.save('qdot/cobyla/var_ideal_UCC.npy',VARs_ideal_UCC)
###np.save('qdot/cobyla/coeff_ideal_UCC.npy',vqe_coeffs_ideal_UCC)
###np.save('qdot/cobyla/E_noisy_UCC.npy',Es_noisy_UCC)
###np.save('qdot/cobyla/var_noisy_UCC.npy',VARs_noisy_UCC)
###np.save('qdot/cobyla/coeff_noisy_UCC.npy',vqe_coeffs_noisy_UCC)
    
