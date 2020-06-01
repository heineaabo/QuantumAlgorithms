from openfermion.hamiltonians import MolecularData
from openfermionpsi4 import run_psi4
import numpy as np

import sys
import sys
sys.path.append('../..')
from vqe import *
sys.path.append('../../ansatz')
from ucc import UnitaryCoupledCluster

sys.path.append('../../../../QuantumCircuitOptimizer')
from quantum_circuit import QuantumCircuit,SecondQuantizedHamiltonian
from quantum_circuit.utils import molecular2sec_quant

from openfermion.hamiltonians import MolecularData
from openfermion.transforms import get_fermion_operator, get_sparse_operator, jordan_wigner
from openfermion.utils import get_ground_state
from openfermionpsi4 import run_psi4

R = 0.7

geometry = [['H',[0,0,0]],
            ['H',[0,0,R]]]
molecule = MolecularData(geometry,'sto-3g',1,0)
molecule = run_psi4(molecule,run_fci=True)
h_pq = molecule.one_body_integrals
h_pqrs = molecule.two_body_integrals
Enuc = molecule.nuclear_repulsion
fci_energy = molecule.fci_energy
print(fci_energy)

# OpenFermion Hamiltonian
molecular_hamiltonian = molecule.get_molecular_hamiltonian()
fermion_hamiltonian = get_fermion_operator(molecular_hamiltonian)
qubit_hamiltonian = jordan_wigner(fermion_hamiltonian)
qubit_hamiltonian.compress()
#print('The Jordan-Wigner Hamiltonian in canonical basis follows:\n{}'.format(qubit_hamiltonian))
circuit_list = qubit_hamiltonian

#h_pq, h_pqrs = molecular2sec_quant(h_pq,h_pqrs)

Model = SecondQuantizedHamiltonian(2,4,h_pq,h_pqrs,nuclear_repulsion=Enuc,add_spin=True,anti_symmetric=False)
#Model.get_circuit()
circuit_list2 = Model.circuit_list('vqe')
print('OepnFerm')
for i in circuit_list: print(i)
print('New')
for i in circuit_list2: print(i)
