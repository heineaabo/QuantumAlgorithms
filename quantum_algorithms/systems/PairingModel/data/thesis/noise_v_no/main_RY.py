import numpy as np
from tqdm import tqdm

import sys
sys.path.append('../../../')
from matrix import pairing_mat_elems
sys.path.append('../../../../..')
from vqe import VQE
sys.path.append('../../../../../optimizers')
from optimizers import SPSA
sys.path.append('../../../../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian

n = 2
l = 4
delta = 1
g = 1

h,v = pairing_mat_elems(n,l,delta,g)


pairing =  PairingHamiltonian(n,l,h,v)


points = 500
thetas = [np.array([i]) for i in np.arange(1,2*np.pi+1,2*np.pi/points)]


option1 = {'shots':1000,
           'optimization_level':1,
           #'seed':1,
           'print':False}
option2 = {'shots':1000,
           'optimization_level':1,
           'device':'ibmq_essex',
           'layout':[1,0,2,3],
           'noise_model':True,
           'basis_gates':True,
           'coupling_map':True,
           #'seed':1,
           'meas_fitter':False,
           'print':False}
option3 = {'shots':1000,
           'optimization_level':1,
           'device':'ibmq_essex',
           'layout':[1,0,2,3],
           'noise_model':True,
           'basis_gates':True,
           'coupling_map':True,
           #'seed':1,
           'print':False}
options = [option1,option2,option3]
data = np.zeros((3,points))
for j in range(1,3):
    i = 0
    for theta in tqdm(thetas):
        vqe = VQE(pairing,
                  SPSA(),
                  'RY',
                  options=options[j])
        data[j,i] = vqe.expval(theta)
        i += 1
#np.save('normal.npy',data[0])
np.save('ESSEX_noisy_no_meas_fit.npy',data[1])
np.save('ESSEX_noisy_with_meas_fit.npy',data[2])
