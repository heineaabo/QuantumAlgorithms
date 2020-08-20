import numpy as np
import matplotlib.pyplot as plt
import qiskit as qk

from quantum_algorithms import QPE
from quantum_circuit import SecondQuantizedHamiltonian


l = 4     # number of spin orbitals / number of qubits
n = 2     # Number of occupied spin orbitals
delta = 1 # Level spacing
g = 1     # Interaction strength

# Matrix elements
h_pq = np.identity(l)
for p in range(l):
	h_pq[p,p] *= delta*(p - (p%2))/2
	
h_pqrs = np.zeros((l,l,l,l))
for p in range(0,l-1,2):
	for r in range(0,l-1,2):
		h_pqrs[p,p+1,r,r+1] = -0.5*g


def pairing_ansatz(qc,qb,n_simulation):
    for qbit in range(0,n_simulation,2):
        qc.h(qb[qbit])
        qc.cx(qb[qbit],qb[qbit+1])
    return qc

n_work = 6
Emax = 2

t = 0.49*np.pi
rho = 10

print('Setting up Hamiltonian')
pairing_model = SecondQuantizedHamiltonian(n,l,h_pq,h_pqrs,exp=True)
print('Setting up QPE class')
qpe = QPE(pairing_model,pairing_ansatz,n_work,Emax,t=t,rho=rho)
print('Running...')
qpe.estimate()

plt.plot(qpe.x,qpe.y)
plt.show()
