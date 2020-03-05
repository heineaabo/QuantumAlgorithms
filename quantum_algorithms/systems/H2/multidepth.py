import numpy as np
import matplotlib.pyplot as plt 

import sys
sys.path.append('../..')
sys.path.append('../../../../QuantumCircuitOptimizer')
from vqe import *
from ucc_nonopt import UnitaryCoupledCluster
#from spsa import SPSA

from openfermion.hamiltonians import MolecularData
from openfermionpsi4 import run_psi4

from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian

def get_H2_info(R,
                       basis='sto-3g',
                       multiplicity=1,
                       charge=0):
    geometry = [['H',[0,0,0]],
                ['H',[0,0,R]]]
    h2_molecule = MolecularData(geometry,basis,multiplicity,charge)

    h2_molecule = run_psi4(h2_molecule,
                            run_mp2=True,
                            run_cisd=True,
                            run_ccsd=True,
                            run_fci=True)
    one_body = h2_molecule.one_body_integrals
    two_body = h2_molecule.two_body_integrals
    Enuc = h2_molecule.nuclear_repulsion
    
    energies = {}
    energies['fci'] = h2_molecule.fci_energy
    energies['hf'] = h2_molecule.hf_energy
    energies['cisd'] = h2_molecule.cisd_energy
    energies['ccsd'] = h2_molecule.ccsd_energy
    #print('electron:',h2_molecule.n_electrons)

    return one_body, two_body, Enuc, energies

l = 4
n = 2

max_iter = 1000
max_evals = 1000

R = 0.7

print('For bond length {}'.format(R))
# Get psi4 
h_pq,h_pqrs,Enuc,energies = get_H2_info(R)
# Save energies
cisd = energies['cisd']
ccsd = energies['ccsd']
fci = energies['fci']
hf = energies['hf']

# Prepare circuit list
H2 = SecondQuantizedHamiltonian(n,l)
H2.set_integrals(h_pq,h_pqrs,Enuc,add_spin=True)
H2.get_circuit()
circuit_list = H2.to_circuit_list(ptype='vqe')

max_depth = 10
depths = np.arange(max_depth)
Es = np.zeros(max_depth)
state_count = np.zeros((max_depth,2)) # [i,0] -> Legal states, [i,1] -> illegal states

for i,depth in enumerate(depths):
    print('{} / {}'.format(i+1,max_depth))
    ansatz = UnitaryCoupledCluster(n,l,'D',depth=depth)
    og_params = ansatz.new_parameters(H2.h,
                                      H2.v)

    theta = og_params
    vqe = VQE(n_qubits = l,
        ansatz = ansatz,
        circuit_list = circuit_list,
        shots = 1000,
        ancilla=0,
        max_energy=False,
        prnt=False,
        count_states=True)
    Es[i] = vqe.expval(og_params)
    state_count[i,0] = vqe.legal
    state_count[i,1] = vqe.illegal

np.save('data/UCCdepth/depths.npy',depths)
np.save('data/UCCdepth/legal_unopt.npy',state_count[:,0])
np.save('data/UCCdepth/illegal_unopt.npy',state_count[:,1])

