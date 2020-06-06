import numpy as np
from tqdm import tqdm

import sys
sys.path.append('..')
from vqe import VQE
from fci import FCI
from ansatz import RYRZ,UnitaryCoupledCluster
from systems import get_h2_matrix,get_qdot_matrix
sys.path.append('../optimizers')
from optimizer import Minimizer,QKSPSA
sys.path.append('../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian


l = 6     # number of spin orbitals / number of qubits
n = 4     # Number of occupied spin orbitals

shots = 5000
option_ideal = {'shots':shots,'print':True}
option_noisy_RYRZ = {'shots':shots,'print':False,
            'device':'ibmq_london','noise_model':True,'basis_gates':True,
            'meas_fit':True}
option_noisy_UCC = {'shots':shots,'print':False,
            'device':'ibmq_london','noise_model':True,'basis_gates':True,
            'meas_fit':True}

######## QDOT
omegas = np.arange(0.1,2.5,0.1)

FCIs = np.zeros_like(omegas)
HFs = np.zeros_like(omegas)
CCs = np.zeros_like(omegas)
fci_coeffs = np.zeros((len(omegas),15)) 

# RYRZ
Es_ideal_RYRZ = np.zeros_like(omegas)
VARs_ideal_RYRZ = np.zeros_like(omegas)
vqe_coeffs_ideal_RYRZ = np.zeros((len(omegas),15)) 

# UCC
Es_ideal_UCC = np.zeros_like(omegas)
VARs_ideal_UCC = np.zeros_like(omegas)
vqe_coeffs_ideal_UCC = np.zeros((len(omegas),15)) 

i = 0
for omega in tqdm(omegas):
    h,v,E = get_qdot_matrix(n,l,omega)
    ev,vv = FCI(n,l,h,v,ret_all=True)
    FCIs[i] = ev[0]
    HFs[i] = E[0].real
    CCs[i] = E[1].real
    fci_coeffs[i] = vv[:,0].real
    h2 = SecondQuantizedHamiltonian(n,l,h,v,anti_symmetric=True)

    h2.group_paulis()
    
    #### RYRZ
    # Ideal
    model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'RYRZ',ansatz_depth=2,options=option_noisy_RYRZ)
    theta = model.optimize()
    Es_ideal_RYRZ[i],VARs_ideal_RYRZ[i] = model.get_mean(theta,N=10000,M=10)
    vqe_coeffs_ideal_RYRZ[i] = model.get_state_coeff(theta)

    ### UCC
    # Ideal
    model = VQE(h2,Minimizer('Cobyla',tol=1e-05,disp=False),'UCCSDr',options=option_noisy_UCC)
    theta = model.optimize()
    Es_ideal_UCC[i],VARs_ideal_UCC[i] = model.get_mean(theta,N=10000,M=10)
    vqe_coeffs_ideal_UCC[i] = model.get_state_coeff(theta)

    i += 1

np.save('qdot/n4l6/omegas.npy',omegas)
np.save('qdot/n4l6/fci.npy',FCIs)
np.save('qdot/n4l6/hf.npy',HFs)
np.save('qdot/n4l6/cc.npy',CCs)
np.save('qdot/n4l6/fci_coeff.npy',fci_coeffs)
## RYRZ
#np.save('qdot/n4l6/E_ideal_RYRZ.npy',Es_ideal_RYRZ)
#np.save('qdot/n4l6/var_ideal_RYRZ.npy',VARs_ideal_RYRZ)
#np.save('qdot/n4l6/coeff_ideal_RYRZ.npy',vqe_coeffs_ideal_RYRZ)
## UCC
np.save('qdot/n4l6/E_ideal_UCCSD.npy',Es_ideal_UCC)
np.save('qdot/n4l6/var_ideal_UCCSD.npy',VARs_ideal_UCC)
np.save('qdot/n4l6/coeff_ideal_UCCSD.npy',vqe_coeffs_ideal_UCC)
    
