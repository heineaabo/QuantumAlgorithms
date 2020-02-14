import numpy as np
import qiskit as qk

# PATH
import sys
sys.path.append('../..')

from ucc import *
from vqe import *

from openfermion.hamiltonians import MolecularData
from openfermion.transforms import get_fermion_operator, get_sparse_operator, jordan_wigner
from openfermion.utils import get_ground_state
from openfermionpsi4 import run_psi4

from quantum_circuit import QuantumCircuit,Hamiltonian
from quantum_circuit.utils import molecular2sec_quant

# LiH atom
diatomic_bond_length = 1.45
geometry = [('Li', (0., 0., 0.)), ('H', (0., 0., diatomic_bond_length))]
basis = 'sto-3g'
multiplicity = 1
molecule = MolecularData(geometry, basis, multiplicity, description=str(diatomic_bond_length))
#molecule.load()
molecule = run_psi4(molecule,
                    run_scf=True,
                    run_fci=True)
h_pq = molecule.one_body_integrals
h_pqrs = molecule.two_body_integrals
nuc_rep = molecule.nuclear_repulsion
fci = molecule.fci_energy
hf = molecule.hf_energy

# OpenFermion Hamiltonian
active_space_start = 1
active_space_stop = 3
molecular_hamiltonian = molecule.get_molecular_hamiltonian(
    occupied_indices=range(active_space_start),
    active_indices=range(active_space_start, active_space_stop))
fermion_hamiltonian = get_fermion_operator(molecular_hamiltonian)
qubit_hamiltonian = jordan_wigner(fermion_hamiltonian)
qubit_hamiltonian.compress()
#print('The Jordan-Wigner Hamiltonian in canonical basis follows:\n{}'.format(qubit_hamiltonian))

#sparse_hamiltonian = get_sparse_operator(qubit_hamiltonian)
#energy, state = get_ground_state(sparse_hamiltonian)
#print('Ground state energy is {} Hartree.\n'.format(energy))
#print('Ground state energy with FCI is {} Hartree.\n'.format(molecule.fci_energy))

print(nuc_rep)
print('FCI:',fci)
l = molecule.n_electrons*2
n = molecule.n_electrons
print('system size: {}, {}'.format(n,l))
h_pq[np.absolute(h_pq) < 1e-8] = 0.
print(h_pq)
print(h_pqrs.shape)
sys.exit()

LiH = Hamiltonian(n,l)
one,two = molecular2sec_quant(h_pq,h_pqrs)
LiH.set_integrals(h_pq,h_pqrs)#,nuclear_repulsion=nuc_rep)
LiH.get_circuit()
circuit_list = LiH.circuit_list
print(LiH)
print(qubit_hamiltonian)

ansatz = UnitaryCoupledCluster(n,l,'SD',1)
theta = ansatz.new_parameters()

vqe_LiH = VQE(n_qubits = l,
    ansatz = ansatz,
    circuit_list = circuit_list,
    shots = 500,
    ancilla=0,
    max_energy=False,
    prnt=True)

t = vqe_LiH.optimize_classical(theta,method='Powell')
