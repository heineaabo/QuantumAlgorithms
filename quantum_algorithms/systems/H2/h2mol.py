import numpy as np
import matplotlib.pyplot as plt 

import sys
sys.path.append('../..')
sys.path.append('../../../../QuantumCircuitOptimizer')
from vqe import *
from ucc import UnitaryCoupledCluster
from spsa import SPSA

from openfermion.hamiltonians import MolecularData
from openfermionpsi4 import run_psi4

from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian
from quantum_circuit.utils import molecular2sec_quant

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
    print('electron:',h2_molecule.n_electrons)

    return one_body, two_body, Enuc, energies

def print_mat(mat):
    if mat.shape == (2,2,2,2):
        p,q,r,s = mat.shape
        for i in range(p):
            for j in range(q):
                for k in range(r):
                    #for l in range(s):
                    print('|{:0.4f}, {:0.4f}|'.format(mat[i,j,k,0],
                                                      mat[i,j,k,1]),end='')
                print('')
            print('\n')
    if mat.shape == (4,4,4,4):
        p,q,r,s = mat.shape
        for i in range(p):
            for j in range(q):
                for k in range(r):
                    #for l in range(s):
                    print('|{:0.4f}, {:0.4f}, {:0.4f}, {:0.4f}|'.format(mat[i,j,k,0],
                                                                    mat[i,j,k,1],
                                                                    mat[i,j,k,2],
                                                                    mat[i,j,k,3]),end='')
                print('')
            print('\n')

def print_non_zero(mat):
    p,q,r,s = mat.shape
    for i in range(p):
        for j in range(q):
            for k in range(r):
                for l in range(s):
                    if not np.isclose(mat[i,j,k,l],0):
                        print('{:0.4f} -> {}, {}, {}, {}'.format(mat[i,j,k,l],i,j,k,l))

l = 4
n = 2


h_pq,h_pqrs,Enuc,energies = get_H2_info(0.7)
print('HF energy:',energies['hf'])
print('CISD energy:',energies['cisd'])
print('CCSD energy:',energies['ccsd'])
print('FCI energy:',energies['fci'])
print('Nuclear repulsion:',Enuc)
print('diff',energies['fci']+Enuc)
print_mat(h_pqrs)
print_non_zero(h_pqrs)
h_pq,h_pqrs = molecular2sec_quant(h_pq,h_pqrs)
#print(h_pq.shape)
print_mat(h_pqrs)
#print(h_pqrs.shape)
print_non_zero(h_pqrs)

#H2 = Hamiltonian(n,l)
#H2.set_integrals(h_pq,h_pqrs,Enuc)
#H2.get_circuit()
#circuit_list = H2.circuit_list
#for i in circuit_list: print(i)
#
#ansatz = UnitaryCoupledCluster(n,l,'SD',1)
#theta = ansatz.new_parameters()
#
#
#VQE_H2 = VQE(n_qubits = l,
#    ansatz = ansatz,
#    circuit_list = circuit_list,
#    shots = 500,
#    ancilla=0,
#    max_energy=False,
#    prnt=False)
#
##t = VQE_H2.optimize_classical(theta,method='Powell')
#
#options = {'feedback':1,'grad_avg':5}
#optimization = SPSA(VQE_H2.expval,
#                    theta,
#                    min_change=0.1,
#                    noise_var = 0.01,
#                    options=options)
#optimization.run()
#


