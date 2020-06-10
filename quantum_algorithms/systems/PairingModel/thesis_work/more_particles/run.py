import numpy as np
from tqdm import tqdm

import sys
sys.path.append('../..')
from matrix import get_pairing_matrix
sys.path.append('../../../..')
from vqe import VQE
sys.path.append('../../../../optimizers')
from optimizer import Minimizer
sys.path.append('../../../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian


l = 8     # number of spin orbitals / number of qubits
n = 4     # Number of occupied spin orbitals
delta = 1 # Level spacing


gs = np.arange(-2,2,0.4)
shots = 8000
#ansatze = ['RYRZ','UCCD']
ansatze = ['RYRZ']
depth = 1
#for g in gs:
for ansatz in ansatze:
    if ansatz == 'RYRZ':
        depth = 2
    else:
        depth = 1
    print('For ',ansatz,'with depth',depth)
    i = 0
    fci  = np.zeros_like(gs)
    hf  = np.zeros_like(gs)
    data = np.zeros_like(gs)
    coeff = np.zeros((len(gs),70))
    for g in tqdm(gs):
        h_pq,h_pqrs,Efci,Ehf = get_pairing_matrix(n,l,delta,g,w_E=True)
        fci[i] = Efci
        hf[i] = Ehf
        #print('FCI:',fci[i])
        #print('HF:',Ehf)
        pairing = PairingHamiltonian(n,l,h_pq,h_pqrs)
        #pairing = SecondQuantizedHamiltonian(n,l,h_pq,h_pqrs)
        pairing.group_paulis()
        #print('Num measures:',len(pairing.circuit_list('vqe')))
        options = {'shots':shots,
                   #'optimization_level':1,
                   #'seed':1,
                   'print':False}
        model = VQE(pairing,Minimizer('Cobyla',tol=1e-05,disp=False),ansatz,ansatz_depth=depth,options=options)
        theta = model.optimize()
        data[i],_ = model.get_mean(theta)
        coeff[i] = model.get_state_coeff(theta)
        i += 1
    np.save('E_{}_{}shots.npy'.format(ansatz,shots),data)
    np.save('coeff_{}_{}shots.npy'.format(ansatz,shots),coeff)
    np.save('fci.npy',fci)

