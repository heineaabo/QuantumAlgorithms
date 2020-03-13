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

def get_avg(L, num):
    L = np.asarray(L[len(L)-num:])
    return np.mean(L),np.std(L)


l = 4
n = 2

max_iter = 100
max_evals = 100

methods = ['Powell','Nelder-Mead','Cobyla','SPSA']
bonds = np.linspace(0.5,1.5,11)

cisd = np.zeros(11)
ccsd = np.zeros(11)
fci = np.zeros(11)
hf = np.zeros(11)

result = {'Powell':[],
          'Nelder-Mead':[],
          'Cobyla':[],
          'SPSA':[]}
result_avg = {'Powell':[],
              'Nelder-Mead':[],
              'Cobyla':[],
              'SPSA':[]}
result_std = {'Powell':[],
              'Nelder-Mead':[],
              'Cobyla':[],
              'SPSA':[]}

for i,R in enumerate(bonds):
    print('For bond length {}'.format(R))
    # Get psi4 
    h_pq,h_pqrs,Enuc,energies = get_H2_info(R)
    # Save energies
    cisd[i] = energies['cisd']
    ccsd[i] = energies['ccsd']
    fci[i] = energies['fci']
    hf[i] = energies['hf']

    # Prepare circuit list
    H2 = SecondQuantizedHamiltonian(n,l,h_pq,h_pqrs,add_spin=True)

    options = {'shots':500}

    for method in methods:
        print('    - now running {}'.format(method),end='')
        vqe = VQE(H2,options=options)
        theta = vqe.theta
        optimizer=None
        if method == 'Powell':
            optimizer = POWELL(maxfev=max_evals)
        if method == 'Nelder-Mead':
            optimizer = NELDER_MEAD(maxiter=max_iter,maxfev=max_evals)
        if method == 'Cobyla':
            optimizer = COBYLA(maxiter=max_evals)
        if method == 'SPSA':
            optimizer = SPSA(max_trials=int(max_evals/2))
        theta,E,_ = optimizer.optimize(len(theta),vqe.expval,initial_point=theta) 
        result[method].append(E)
        avg,std = get_avg(vqe.energies,10)
        result_avg[method].append(avg)
        result_std[method].append(std)
        print(', {} function evaluations'.format(vqe.evals))
        

[np.save('data/opt100evals/{}500.npy'.format(key),np.asarray(item)) for key,item in result.items()]
[np.save('data/opt100evals/{}500_avg10.npy'.format(key),np.asarray(item)) for key,item in result_avg.items()]
[np.save('data/opt100evals/{}500_std10.npy'.format(key),np.asarray(item)) for key,item in result_std.items()]
np.save('data/opt100evals/cisd.npy',cisd)
np.save('data/opt100evals/ccsd.npy',ccsd)
np.save('data/opt100evals/fci.npy',fci)
np.save('data/opt100evals/hf.npy',hf)

