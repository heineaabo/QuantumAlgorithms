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


l = 6     # number of spin orbitals / number of qubits
n = 4     # Number of occupied spin orbitals

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
omegas = omegas[[i for i in range(0,len(omegas),5)]]

FCIs = np.zeros_like(omegas)
HFs = np.zeros_like(omegas)
CCs = np.zeros_like(omegas)
fci_coeffs = np.zeros((len(omegas),15)) 

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

    h2.group_paulis(qwc=True,gc=True)
    if i == 0: print(h2._circuit_list)

    ### UCC
    # Ideal
    #model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'UCCSD',options=option_ideal)
    #theta = model.optimize()
    #Es_ideal_UCC[i],VARs_ideal_UCC[i] = model.get_mean(theta,N=10000,M=10)
    #vqe_coeffs_ideal_UCC[i] = model.get_state_coeff(theta)

    i += 1

np.save('omegas.npy',omegas)
np.save('fci.npy',FCIs)
np.save('hf.npy',HFs)
np.save('cc.npy',CCs)
np.save('fci_coeff.npy',fci_coeffs)
## UCCSDr
#np.save('UCCSD_E_1.npy',Es_ideal_UCC)
#np.save('UCCSD_var_1.npy',VARs_ideal_UCC)
#np.save('UCCSD_coeff_1.npy',vqe_coeffs_ideal_UCC)

