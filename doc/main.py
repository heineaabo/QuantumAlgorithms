from vqe import *
from Qoperator import *
import numpy as np


def ansatz(theta,qb,qc,cb):
	t = theta.reshape(2,2,2,2)
	n_qubits = 4
	n_fermi = 2
			
	for i in range(n_fermi,n_qubits):
			qc.x(qb[i])

	for i in range(n_fermi):
		for j in range(i+1,n_fermi):
			for a in range(n_fermi,n_qubits):
				for b in range(a+1,n_qubits):
					for k in range(i+1,j):
						qc.z(qb[k])
					for l in range(a+1,b):
						qc.z(qb[l])
					#TERM1
					qc.h(qb[i])
					qc.h(qb[j])
					qc.h(qb[a])
					qc.rz(np.pi/2,qb[b])
					qc.h(qb[b])

					qc.cx(qb[i],qb[n_qubits])
					qc.cx(qb[j],qb[n_qubits])
					qc.cx(qb[a],qb[n_qubits])
					qc.cx(qb[b],qb[n_qubits])

					qc.rz(-t[i,j,a-n_fermi,b-n_fermi]/8,qb[n_qubits])

					qc.cx(qb[i],qb[n_qubits])
					qc.cx(qb[j],qb[n_qubits])
					qc.cx(qb[a],qb[n_qubits])
					qc.cx(qb[b],qb[n_qubits])


					qc.h(qb[i])
					qc.h(qb[j])
					qc.h(qb[a])
					qc.h(qb[b])
					qc.rz(-np.pi/2,qb[b])

					#TERM2

					qc.h(qb[i])
					qc.h(qb[j])
					qc.rz(np.pi/2,qb[a])
					qc.h(qb[a])
					qc.h(qb[b])

					qc.cx(qb[i],qb[n_qubits])
					qc.cx(qb[j],qb[n_qubits])
					qc.cx(qb[a],qb[n_qubits])
					qc.cx(qb[b],qb[n_qubits])

					qc.rz(-t[i,j,a-n_fermi,b-n_fermi]/8,qb[n_qubits])

					qc.cx(qb[i],qb[n_qubits])
					qc.cx(qb[j],qb[n_qubits])
					qc.cx(qb[a],qb[n_qubits])
					qc.cx(qb[b],qb[n_qubits])


					qc.h(qb[i])
					qc.h(qb[j])
					qc.h(qb[a])
					qc.rz(-np.pi/2,qb[a])
					qc.h(qb[b])


					#TERM3

					qc.h(qb[i])
					qc.rz(np.pi/2,qb[j])
					qc.h(qb[j])
					qc.h(qb[a])
					qc.h(qb[b])

					qc.cx(qb[i],qb[n_qubits])
					qc.cx(qb[j],qb[n_qubits])
					qc.cx(qb[a],qb[n_qubits])
					qc.cx(qb[b],qb[n_qubits])

					qc.rz(t[i,j,a-n_fermi,b-n_fermi]/8,qb[n_qubits])

					qc.cx(qb[i],qb[n_qubits])
					qc.cx(qb[j],qb[n_qubits])
					qc.cx(qb[a],qb[n_qubits])
					qc.cx(qb[b],qb[n_qubits])


					qc.h(qb[i])
					qc.h(qb[j])
					qc.rz(-np.pi/2,qb[j])
					qc.h(qb[a])
					qc.h(qb[b])

					#TERM4

					qc.h(qb[i])
					qc.rz(np.pi/2,qb[j])
					qc.h(qb[j])
					qc.rz(np.pi/2,qb[a])
					qc.h(qb[a])
					qc.rz(np.pi/2,qb[b])
					qc.h(qb[b])

					qc.cx(qb[i],qb[n_qubits])
					qc.cx(qb[j],qb[n_qubits])
					qc.cx(qb[a],qb[n_qubits])
					qc.cx(qb[b],qb[n_qubits])

					qc.rz(-t[i,j,a-n_fermi,b-n_fermi]/8,qb[n_qubits])

					qc.cx(qb[i],qb[n_qubits])
					qc.cx(qb[j],qb[n_qubits])
					qc.cx(qb[a],qb[n_qubits])
					qc.cx(qb[b],qb[n_qubits])


					qc.h(qb[i])
					qc.h(qb[j])
					qc.rz(-np.pi/2,qb[j])
					qc.h(qb[a])
					qc.rz(-np.pi/2,qb[a])
					qc.h(qb[b])
					qc.rz(-np.pi/2,qb[b])

					#TERM5

					qc.rz(np.pi/2,qb[i])
					qc.h(qb[i])
					qc.h(qb[j])
					qc.h(qb[a])
					qc.h(qb[b])

					qc.cx(qb[i],qb[n_qubits])
					qc.cx(qb[j],qb[n_qubits])
					qc.cx(qb[a],qb[n_qubits])
					qc.cx(qb[b],qb[n_qubits])

					qc.rz(t[i,j,a-n_fermi,b-n_fermi]/8,qb[n_qubits])

					qc.cx(qb[i],qb[n_qubits])
					qc.cx(qb[j],qb[n_qubits])
					qc.cx(qb[a],qb[n_qubits])
					qc.cx(qb[b],qb[n_qubits])


					qc.h(qb[i])
					qc.rz(-np.pi/2,qb[i])
					qc.h(qb[j])
					qc.h(qb[a])
					qc.h(qb[b])

					#TERM6

					qc.rz(np.pi/2,qb[i])
					qc.h(qb[i])
					qc.h(qb[j])
					qc.rz(np.pi/2,qb[a])
					qc.h(qb[a])
					qc.rz(np.pi/2,qb[b])
					qc.h(qb[b])

					qc.cx(qb[i],qb[n_qubits])
					qc.cx(qb[j],qb[n_qubits])
					qc.cx(qb[a],qb[n_qubits])
					qc.cx(qb[b],qb[n_qubits])

					qc.rz(-t[i,j,a-n_fermi,b-n_fermi]/8,qb[n_qubits])

					qc.cx(qb[i],qb[n_qubits])
					qc.cx(qb[j],qb[n_qubits])
					qc.cx(qb[a],qb[n_qubits])
					qc.cx(qb[b],qb[n_qubits])


					qc.h(qb[i])
					qc.rz(-np.pi/2,qb[i])
					qc.h(qb[j])
					qc.h(qb[a])
					qc.rz(-np.pi/2,qb[a])
					qc.h(qb[b])
					qc.rz(-np.pi/2,qb[b])

					#TERM7

					qc.rz(np.pi/2,qb[i])
					qc.h(qb[i])
					qc.rz(np.pi/2,qb[j])
					qc.h(qb[j])
					qc.h(qb[a])
					qc.rz(np.pi/2,qb[b])
					qc.h(qb[b])

					qc.cx(qb[i],qb[n_qubits])
					qc.cx(qb[j],qb[n_qubits])
					qc.cx(qb[a],qb[n_qubits])
					qc.cx(qb[b],qb[n_qubits])

					qc.rz(t[i,j,a-n_fermi,b-n_fermi]/8,qb[n_qubits])

					qc.cx(qb[i],qb[n_qubits])
					qc.cx(qb[j],qb[n_qubits])
					qc.cx(qb[a],qb[n_qubits])
					qc.cx(qb[b],qb[n_qubits])


					qc.h(qb[i])
					qc.rz(-np.pi/2,qb[i])
					qc.h(qb[j])
					qc.rz(-np.pi/2,qb[j])
					qc.h(qb[a])
					qc.h(qb[b])
					qc.rz(-np.pi/2,qb[b])

					#TERM8

					qc.rz(np.pi/2,qb[i])
					qc.h(qb[i])
					qc.rz(np.pi/2,qb[j])
					qc.h(qb[j])
					qc.rz(np.pi/2,qb[a])
					qc.h(qb[a])
					qc.h(qb[b])

					qc.cx(qb[i],qb[n_qubits])
					qc.cx(qb[j],qb[n_qubits])
					qc.cx(qb[a],qb[n_qubits])
					qc.cx(qb[b],qb[n_qubits])

					qc.rz(t[i,j,a-n_fermi,b-n_fermi]/8,qb[n_qubits])

					qc.cx(qb[i],qb[n_qubits])
					qc.cx(qb[j],qb[n_qubits])
					qc.cx(qb[a],qb[n_qubits])
					qc.cx(qb[b],qb[n_qubits])


					qc.h(qb[i])
					qc.rz(-np.pi/2,qb[i])
					qc.h(qb[j])
					qc.rz(-np.pi/2,qb[j])
					qc.h(qb[a])
					qc.rz(-np.pi/2,qb[a])
					qc.h(qb[b])
	return(qb,qc,cb)


n_qubits = 4
n_fermi = 2
delta = 1
g = 1
g /= 4
theta = 2*np.pi*np.random.randn(16)
h_pq = np.identity(n_qubits)
for p in range(n_qubits):
	h_pq[p,p] *= delta*(p - (p%2))/2


	
h_pqrs = np.zeros((n_qubits,n_qubits,n_qubits,n_qubits))
for p in range(0,n_qubits-1,2):
	for r in range(0,n_qubits-1,2):
		h_pqrs[p,p+1,r,r+1] = -0.5*g

Pairing = Hamiltonian(n_qubits)
circuit_list = Pairing.get_circuits(h_pq,h_pqrs)


pairing_vqe = VQE(n_qubits = n_qubits,
				ansatz = ansatz,
				circuit_list = circuit_list,
				shots = 500,
				ancilla=1,
				max_energy=False)

pairing_vqe.optimize_classical(theta,method='Powell')

