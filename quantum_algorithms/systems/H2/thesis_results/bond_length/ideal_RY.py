import numpy as np
from tqdm import tqdm

import sys
sys.path.append('../..')
from matrix import get_h2_matrix
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

Rs = np.arange(0.3,2.3,0.1)

FCIs = np.zeros_like(Rs)
HFs = np.zeros_like(Rs)
CCs = np.zeros_like(Rs)
fci_coeffs = np.zeros((len(Rs),6)) 
Es_ideal = np.zeros_like(Rs)
VARs_ideal = np.zeros_like(Rs)
vqe_coeffs_ideal = np.zeros((len(Rs),6)) 
Es_noisy = np.zeros_like(Rs)
VARs_noisy = np.zeros_like(Rs)
vqe_coeffs_noisy = np.zeros((len(Rs),6)) 

i = 0
for R in tqdm(Rs):
    h,v,Enuc,E = get_h2_matrix(R)
    h2 = SecondQuantizedHamiltonian(n,l,h,v,nuclear_repulsion=Enuc,add_spin=True)
    ev,vv = FCI(n,l,h2.h,h2.v,ret_all=True)
    fci_coeffs[i] = vv[:,0].real
    FCIs[i] = E['fci']
    HFs[i] = E['hf']
    CCs[i] = E['ccsd']


    h2.group_paulis(qwc=True,gc=True)
    # Ideal
    model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'RYPAIRING',options=option_ideal)
    theta = model.optimize()
    Es_ideal[i],VARs_ideal[i] = model.get_mean(theta,N=10000,M=10)
    vqe_coeffs_ideal[i] = model.get_state_coeff(theta)
    # Noisy
    model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'RYPAIRING',options=option_noisy)
    theta = model.optimize()
    Es_noisy[i],VARs_noisy[i] = model.get_mean(theta,N=10000,M=10)
    vqe_coeffs_noisy[i] = model.get_state_coeff(theta)
    i += 1

np.save('RY/bonds.npy',Rs)
np.save('RY/fci.npy',FCIs)
np.save('RY/hf.npy',HFs)
np.save('RY/cc.npy',CCs)
np.save('RY/cobyla/E_ideal.npy',Es_ideal)
np.save('RY/cobyla/var_ideal.npy',VARs_ideal)
np.save('RY/cobyla/coeff_ideal.npy',vqe_coeffs_ideal)
np.save('RY/cobyla/E_noisy.npy',Es_noisy)
np.save('RY/cobyla/var_noisy.npy',VARs_noisy)
np.save('RY/cobyla/coeff_noisy.npy',vqe_coeffs_noisy)
    
    
