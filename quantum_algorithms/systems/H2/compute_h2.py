import numpy as np
import matplotlib.pyplot as plt 

import sys
from openfermion.hamiltonians import MolecularData
from openfermionpsi4 import run_psi4

def get_H2_info(R,
                       basis='sto-3g',
                       multiplicity=1,
                       charge=0):
    geometry = [['H',[0,0,0]],
                ['H',[0,0,R]]]
    h2_molecule = MolecularData(geometry,basis,multiplicity,charge)

    h2_molecule = run_psi4(h2_molecule,
                            run_mp2=True,
                            run_ccsd=True,
                            run_fci=True)
    one_body = h2_molecule.one_body_integrals
    two_body = h2_molecule.two_body_integrals
    Enuc = h2_molecule.nuclear_repulsion
    
    energies = {}
    energies['fci'] = h2_molecule.fci_energy
    energies['hf'] = h2_molecule.hf_energy
    energies['ccsd'] = h2_molecule.ccsd_energy
    #print('electron:',h2_molecule.n_electrons)

    return one_body, two_body, Enuc, energies

l = 4
n = 2

bonds = np.linspace(0.1,3.0,59)

ccsd = np.zeros(len(bonds))
hf = np.zeros(len(bonds))
fci = np.zeros(len(bonds))

for i,R in enumerate(bonds):
    R = np.round(R,3)
    print('For bond length {}'.format(R))
    # Get psi4 
    h_pq,h_pqrs,Enuc,energies = get_H2_info(R)
    np.save('matrix_elements/h_R{}'.format(R),h_pq)
    np.save('matrix_elements/v_R{}'.format(R),h_pqrs)
    # Save energies
    ccsd[i] = energies['ccsd']
    fci[i] = energies['fci']
    hf[i] = energies['hf']

np.save('data/pre_calc/bonds.npy',bonds)
np.save('data/pre_calc/ccsd.npy',ccsd)
np.save('data/pre_calc/hf.npy',hf)
np.save('data/pre_calc/fci.npy',fci)


