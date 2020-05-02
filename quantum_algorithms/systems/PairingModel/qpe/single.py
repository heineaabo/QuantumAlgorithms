import numpy as np
import sys
sys.path.append('../../..')
sys.path.append('../../../../../QuantumCircuitOptimizer')
from qpe import QPE
from quantum_circuit import SecondQuantizedHamiltonian

def pairing_ansatz(qc,qb,n_simulation):
    for qbit in range(0,n_simulation,2):
        qc.h(qb[qbit])
        qc.cx(qb[qbit],qb[qbit+1])
    return qc


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

pairing = SecondQuantizedHamiltonian(n,l,h_pq,h_pqrs,exp=True)

Emax = 2

n_work = 6

t = 0.49*np.pi
n = 100

estimator = QPE(pairing,pairing_ansatz,n_work,Emax,t=t,n=n)
estimator.estimate()

data = np.stack((estimator.x,estimator.y))
np.save('data/qpe_n2_l4_rho100.npy',data)


