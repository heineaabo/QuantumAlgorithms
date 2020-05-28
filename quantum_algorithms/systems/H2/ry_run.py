import numpy as np
from h2getter import get_H2

import sys
sys.path.append('../..')
from vqe import VQE
sys.path.append('../../optimizers')
from optimizers import QuantumGradientDescent,Minimizer,QKSPSA
sys.path.append('../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian,PairingHamiltonian

#h,v,Enuc,fci = get_H2(0.7)
def get_H2_info(R,
                       basis='sto-3g',
                       multiplicity=1,
                       charge=0):
    from openfermion.hamiltonians import MolecularData
    from openfermionpsi4 import run_psi4
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

h,v,Enuc,E = get_H2_info(0.7)
fci = E['fci']
hf = E['hf']


print('HF  energy:',hf)
print('FCI energy:',fci)

mol = SecondQuantizedHamiltonian(2,4,h,v,nuclear_repulsion=Enuc,add_spin=True)

options = {'shots':1000,
           'optimization_level':1,
           #'device':'ibmq_london',
           #'layout':[1,0,2,3],
           #'noise_model':True,
           #'basis_gates':True,
           #'coupling_map':True,
           #'seed':1,
           'print':True}
#options = {'shots':10000,'print':True}
alg = VQE(mol,Minimizer('Cobyla'),'RY',options=options)
alg.optimize()
#alg.expval([0,0,0,0,0])


