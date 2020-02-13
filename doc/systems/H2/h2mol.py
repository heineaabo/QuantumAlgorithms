import numpy as np
import matplotlib.pyplot as plt 

import sys
sya.path.append('..')
from vqe import *
from Qoperator import *
from ucc import UnitaryCoupledCluster

from openfermion.hamiltonians import MolecularData
from openfermionpsi4 import run_psi4

from quantum_circuit import QuantumCircuit,Hamiltonian
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
    
    fci_energy = h2_molecule.fci_energy

    return one_body, two_body, Enuc, fci_energy,


l = 4
n = 2


h_pq,h_pqrs,Enuc,fci = get_H2_info(0.6)
print('FCI energy:',fci)

h_pq,h_pqrs = molecular2sec_quant(h_pq,h_pqrs)

H2 = Hamiltonian(n,l)
H2.set_integrals(h_pq,h_pqrs,Enuc)
H2.get_circuit()
circuit_list = H2.to_circuit_list(ptype='vqe')

ansatz = UnitaryCoupledCluster(n,l,'SD',1)
theta = ansatz.new_parameters()


VQE_H2 = VQE(n_qubits = l,
    ansatz = ansatz,
    circuit_list = circuit_list,
    shots = 500,
    ancilla=0,
    max_energy=False)

t = VQE_H2.optimize_classical(theta,method='Powell')















#Rs = np.load('systems/H2/data/Rs.npy')
#energies = np.zeros_like(Rs)
#FCIs = np.load('systems/H2/data/fci_energies.npy')
#for i,R in enumerate(Rs):
#    fci = FCIs[i]
#
#    h_pq = np.load('systems/H2/data/R{}_one.npy'.format(i))
#    h_pqrs = np.load('systems/H2/data/R{}_two.npy'.format(i))
#    #print(h_pq)
#    #print(h_pqrs)
#    #print(h_pqrs[0,1,1,0])
#
#    h_pq, h_pqrs = molecular2sec_quant(h_pq,h_pqrs)
#
#    Model = Hamiltonian(l)
#    circuit_list = Model.get_circuits(h_pq,h_pqrs)
#    for i in circuit_list: print(i)
#
#    UCCS = UnitaryCoupledCluster(n,l,'SD',1)
#    theta = UCCS.new_parameters()
#
#    VQE_H2 = VQE(n_qubits = l,
#        ansatz = UCCS,#ansatz,
#        circuit_list = circuit_list,
#        shots = 500,
#        ancilla=0,
#        max_energy=False)
#
#    t = VQE_H2.optimize_classical(theta,method='Cobyla')
#    energies[i] = VQE_H2.expval(t)
#
#plt.plot(Rs,FCIs,label='FCI')
#plt.plot(Rs,energies,label='VQE')
#plt.legend()
#plt.xlabel('Bond length [Å]')
#plt.xlabel('Energy [Hartrees]')
#plt.title('Ground state energy of H2 molecule')
#plt.show()
