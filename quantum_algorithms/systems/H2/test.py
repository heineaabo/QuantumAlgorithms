from openfermion.hamiltonians import MolecularData
from openfermionpsi4 import run_psi4
import numpy as np

import sys
import sys
sys.path.append('../..')
from vqe import *
from ucc import UnitaryCoupledCluster
from spsa import SPSA

sys.path.append('../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian
from quantum_circuit.utils import molecular2sec_quant

from openfermion.hamiltonians import MolecularData
from openfermion.transforms import get_fermion_operator, get_sparse_operator, jordan_wigner
from openfermion.utils import get_ground_state
from openfermionpsi4 import run_psi4

R = 1.5

geometry = [['H',[0,0,0]],
            ['H',[0,0,R]]]
molecule = MolecularData(geometry,'sto-3g',1,0)
molecule = run_psi4(molecule,run_fci=True)
h_pq = molecule.one_body_integrals
h_pqrs = molecule.two_body_integrals
repulsion = molecule.nuclear_repulsion
fci_energy = molecule.fci_energy
print(fci_energy)

# OpenFermion Hamiltonian
molecular_hamiltonian = molecule.get_molecular_hamiltonian()
fermion_hamiltonian = get_fermion_operator(molecular_hamiltonian)
qubit_hamiltonian = jordan_wigner(fermion_hamiltonian)
qubit_hamiltonian.compress()
#print('The Jordan-Wigner Hamiltonian in canonical basis follows:\n{}'.format(qubit_hamiltonian))
circuit_list = qubit_hamiltonian


#print(overlap)
#print(h_pq)
#
#for p in range(4):
#    for q in range(4):
#        for i in range(2):
#            h_pq[p,q] -= 0.5*h_pqrs[p,i,q,i]

#h_pq, h_pqrs = molecular2sec_quant(h_pq,h_pqrs)

Model = SecondQuantizedHamiltonian(2,4)
Model.set_integrals(h_pq,h_pqrs,repulsion,add_spin=True)
Model.get_circuit()
circuit_list2 = Model.to_circuit_list(ptype='openfermion')
print('OepnFerm')
for i in circuit_list: print(i)
print('New')
for i in circuit_list2: print(i)

#print('Open Fermion:')
#for i in circuit_list:
#    print(i)
#print('\n')
#print('My:')
#for i in circuit_list2:
#    print(i)
#print(h_pq)
#print(h_pqrs)
#for i in circuit_list:
#    for j in circuit_list2:
#        if i[1:] == j[1:]:
#            print('{:5f} {:5f} {:5f} {:5f} {}'.format(i[0],j[0].real,j[0].real*2,j[0].real*4,i[1:]))

#one = np.load('QS_one0714.npy')
#two = np.load('QS_two0714.npy')
#print(one,two)
#
#h2_qs = Hamiltonian(4)
##one_two = molecular2sec_quant(one,two)
#circuit_list3 = h2_qs.get_circuits(one,two)


