import numpy as np
from tqdm import tqdm

import sys
sys.path.append('../..')
from vqe import VQE
from fci import FCI
from systems import get_h2_matrix,get_qdot_matrix
sys.path.append('../../optimizers')
from optimizer import Minimizer,QKSPSA
sys.path.append('../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian


l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals

shots = 1
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


R = 1.0
h,v,Enuc,E = get_h2_matrix(R)
h2 = SecondQuantizedHamiltonian(n,l,h,v,nuclear_repulsion=Enuc,add_spin=True)
h2.group_paulis(qwc=True,gc=True)

# Ideal
ry = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'RYPAIRING',options=option_ideal)
ryrz = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'RYRZ',ansatz_depth=2,options=option_ideal)
uccd = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'UCCD',options=option_ideal)
uccdr = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'UCCDr',options=option_ideal)
uccsd = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'UCCSD',options=option_ideal)
uccsdr = VQE(h2,Minimizer('Cobyla',tol=1/(100*shots),disp=False),'UCCSDr',options=option_ideal)

# RY
ry.print_circ_info()
# RY
ryrz.print_circ_info()
# UCCD
uccd.print_circ_info()
# UCCDr
uccdr.print_circ_info()
# UCCSD
uccsd.print_circ_info()
# UCCSDr
uccsdr.print_circ_info()

