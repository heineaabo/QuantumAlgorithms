import numpy as np

import sys
sys.path.append('../../..')
from matrix import get_h2_matrix
sys.path.append('../../../../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian,CircuitList,PauliString,X

l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals

Rs = [0.4,0.74,1.0,1.6,2.2,2.8]

fci = np.zeros_like(Rs)
hf  = np.zeros_like(Rs)
vqe_factors  = np.zeros_like(Rs)

for i,R in enumerate(Rs):
    h,v,Enuc,E = get_h2_matrix(R)
    fci[i] = E['fci']
    hf[i] = E['hf']
    
    h2 = SecondQuantizedHamiltonian(n,l,h,v,nuclear_repulsion=Enuc,add_spin=True)
    for p in h2.circuit_list('vqe'):
        if len(p) == 0:
            print('factor:',p.factor)
            vqe_factors[i] = p.factor

np.save('fci.npy',fci)
np.save('hf.npy',hf)
np.save('vqe_factors.npy',vqe_factors)

