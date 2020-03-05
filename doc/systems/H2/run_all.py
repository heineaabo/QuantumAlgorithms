import numpy as np
import matplotlib.pyplot as plt 

import sys
sys.path.append('../..')
sys.path.append('../../../../QuantumCircuitOptimizer')
from vqe import *
from ucc import UnitaryCoupledCluster
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

l = 4
n = 2

max_iter = 1000
max_evals = 1000

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
    H2 = SecondQuantizedHamiltonian(n,l)
    H2.set_integrals(h_pq,h_pqrs,Enuc,add_spin=True)
    H2.get_circuit()
    circuit_list = H2.to_circuit_list(ptype='vqe')

    ansatz = UnitaryCoupledCluster(n,l,'D',1)
    og_params = ansatz.new_parameters(H2.h,
                                      H2.v)
    for method in methods:
        print('    - now running {}'.format(method),end='')
        theta = og_params
        vqe = VQE(n_qubits = l,
            ansatz = ansatz,
            circuit_list = circuit_list,
            shots = 1000,
            ancilla=0,
            max_energy=False,
            prnt=False)
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
        print(', {} function evaluations'.format(vqe.evals))
        

[np.save('data/{}.npy'.format(key),np.asarray(item)) for key,item in result.items()]
np.save('data/cisd.npy',cisd)
np.save('data/ccsd.npy',ccsd)
np.save('data/fci.npy',fci)
np.save('data/hf.npy',hf)

