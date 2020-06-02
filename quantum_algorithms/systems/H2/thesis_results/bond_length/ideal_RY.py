import numpy as np
from tqdm import tqdm

import sys
sys.path.append('../..')
from h2getter import get_H2
sys.path.append('../../../..')
from vqe import VQE
sys.path.append('../../../../optimizers')
from optimizer import Minimizer
sys.path.append('../../../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian


l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals

shots = 1000
options = {'shots':shots,'print':False}

Rs = np.arange(0.3,2.3,0.1)

FCIs = np.zeros_like(Rs)
HFs = np.zeros_like(Rs)
CCs = np.zeros_like(Rs)
Es = np.zeros_like(Rs)
VARs = np.zeros_like(Rs)
vqe_coeffs = np.zeros((len(Rs),6)) 
fci_coeffs = np.zeros((len(Rs),6)) 

i = 0
for R in tqdm(Rs):
    h,v,Enuc,E = get_H2(R)
    ev,vv = FCI(h,v)
    fci_coeffs[i] = vv[:,0]
    FCIs[i] = E['fci']
    HFs[i] = E['hf']
    CCs[i] = E['ccsd']

    h2 = SecondQuantizedHamiltonian(n,l,h,v,nuclear_repulsion=Enuc,add_spin=True)

    h2.group_paulis(qwc=True,gc=True)

    model = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'RYPAIRING',options=options)
    theta = model.optimize()
    Es[i],VARs[i] = model.get_mean(theta,N=10000,M=10)
    vqe_coeffs[i] = model.get_state_coeffs(theta)
    i += 1

np.save('RY/bonds.npy',Rs)
np.save('RY/fci.npy',FCIs)
np.save('RY/hf.npy',HFs)
np.save('RY/cc.npy',CCs)
np.save('RY/cobyla/E.npy',Es)
np.save('RY/cobyla/var.npy',VARs)
    
    
