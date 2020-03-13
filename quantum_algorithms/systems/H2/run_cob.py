import numpy as np
import matplotlib.pyplot as plt 

import sys
sys.path.append('../..')
sys.path.append('../../../../QuantumCircuitOptimizer')
from vqe import *
#from spsa import SPSA
from qiskit.aqua.components.optimizers import COBYLA,SPSA,POWELL,NELDER_MEAD

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

def get_avg(L, num):
    L = np.asarray(L[len(L)-num:])
    return np.mean(L),np.std(L)


l = 4
n = 2

max_iter = 1 #200
max_evals = 1 # 200

bonds = np.linspace(0.5,1.5,11)

result = []
results_evals = []

for i,R in enumerate(bonds):
    print('For bond length {}'.format(R))
    # Get psi4 
    h_pq,h_pqrs,Enuc,energies = get_H2_info(R)
    # Save energies
    ccsd[i] = energies['ccsd']
    fci[i] = energies['fci']
    hf[i] = energies['hf']

    # Prepare circuit list
    H2 = SecondQuantizedHamiltonian(n,l,h_pq,h_pqrs,
                                    nuclear_repulsion=Enuc,add_spin=True)

    options = {'shots':1000}

    for method in methods:
        print('    - now running {}'.format(method),end='')
        vqe = VQE(H2,options=options)
        theta = vqe.theta
        optimizer = COBYLA(maxiter=max_evals)
        theta,E,_ = optimizer.optimize(len(theta),vqe.expval,initial_point=theta) 
        result.append(E)
        results_evals.append(vqe.energies)
        print(', {} function evaluations'.format(vqe.evals))
        
np.save('data/optimizers/cobyla/200evals/1000shots/final.npy',np.asarray(result))
max_len = max([len(i) for i in results_evals])
evals = np.zeros((len(results_evals),max_len))
for i in range(len(results_evals)):
    num = len(results_evals[i])
    evals[i,num:] = results_evals[i]
np.save('data/optimizers/cobyla/200evals/1000shots/evals.npy',np.asarray(evals))

