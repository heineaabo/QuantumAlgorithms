import numpy as np
from tqdm import tqdm

import sys
sys.path.append('../..')
from matrix import get_qdot_matrix
sys.path.append('../../../..')
from vqe import VQE
from fci import FCI
sys.path.append('../../../../optimizers')
from optimizer import Minimizer
sys.path.append('../../../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian


l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals

shots = 1000
option_ideal = {'shots':shots,'print':False}
option_noisy = {'shots':shots,'print':False,
            'device':'ibmq_london','noise_model':True,'basis_gates':True,
            'coupling_map':True,'layout':[1,0,2,3],'meas_fit':True}

omegas = np.arange(0.1,3.0,0.1)

FCIs = np.zeros_like(omegas)
HFs = np.zeros_like(omegas)
fci_coeffs = np.zeros((len(omegas),6)) 
Es_ideal = np.zeros_like(omegas)
VARs_ideal = np.zeros_like(omegas)
vqe_coeffs_ideal = np.zeros((len(omegas),6)) 
Es_noisy = np.zeros_like(omegas)
VARs_noisy = np.zeros_like(omegas)
vqe_coeffs_noisy = np.zeros((len(omegas),6)) 

i = 0
for omega in tqdm(omegas):
    h,v,E = get_qdot_matrix(n,l,omega)
    ev,vv = FCI(n,l,h,v,ret_all=True)
    FCIs[i] = ev[0]
    HFs[i] = E[0].real
    fci_coeffs[i] = vv[:,0].real
    h2 = SecondQuantizedHamiltonian(n,l,h,v,anti_symmetric=True)

    h2.group_paulis(qwc=True,gc=True)
    # Ideal
    model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'UCCSD',options=option_ideal)
    theta = model.optimize()
    Es_ideal[i],VARs_ideal[i] = model.get_mean(theta,N=10000,M=10)
    vqe_coeffs_ideal[i] = model.get_state_coeff(theta)
    # Noisy
    model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'UCCSD',options=option_noisy)
    theta = model.optimize()
    Es_noisy[i],VARs_noisy[i] = model.get_mean(theta,N=10000,M=10)
    vqe_coeffs_noisy[i] = model.get_state_coeff(theta)
    i += 1

np.save('UCC/bonds.npy',omegas)
np.save('UCC/fci.npy',FCIs)
np.save('UCC/hf.npy',HFs)
np.save('UCC/cobyla/E_ideal.npy',Es_ideal)
np.save('UCC/cobyla/var_ideal.npy',VARs_ideal)
np.save('UCC/cobyla/coeff_ideal.npy',vqe_coeffs_ideal)
np.save('UCC/cobyla/E_noisy.npy',Es_noisy)
np.save('UCC/cobyla/var_noisy.npy',VARs_noisy)
np.save('UCC/cobyla/coeff_noisy.npy',vqe_coeffs_noisy)
    
    
