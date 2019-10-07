import qiskit as qk
import numpy as np
from scipy.optimize import minimize
from Qoperator import *
qk.IBMQ.load_accounts()




class VQE_UCC:
	def __init__(self,n_qubits,n_fermi,circuit_list,shots=1000):
		self.n_qubits = n_qubits
		self.n_fermi = n_fermi
		self.shots=shots
		self.delta=delta
		self.circuit_list = circuit_list

	def wavefunction_ansatz(self,qc,qb,cb,t,measure=False):

		n_qubits = self.n_qubits
		n_fermi = self.n_fermi
		t = t.reshape(self.n_fermi,self.n_fermi,self.n_qubits - self.n_fermi,self.n_qubits - self.n_fermi)
		
		
		for i in range(n_fermi,n_qubits):
			qc.x(qb[i])

		for i in range(n_fermi):
			for j in range(i+1,n_fermi):
				for a in range(n_fermi,n_qubits):
					for b in range(a+1,n_qubits):
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
		if measure:
			qc.measure(qb,cb)
			job = qk.execute(qc, backend = qk.Aer.get_backend('qasm_simulator'), shots=self.shots)
			result = job.result()
			self.result = result.get_counts(qc)

	def calculate_energy(self,t,circuit_list):
		n_nodes = self.n_fermi
		n_qubits = self.n_qubits
		
		E = 0
		for circuit in circuit_list:
			qb = qk.QuantumRegister(n_qubits+1)
			cb = qk.ClassicalRegister(n_qubits+1)
			qc = qk.QuantumCircuit(qb,cb)
			self.wavefunction_ansatz(qc,qb,cb,t)
			indices = []
			count = 0
			for i in range(n_qubits):
				operation = circuit.get(i).op
				if operation == '':
					count += 1
					continue
				exec('qc.{}(qb[{}])'.format(operation,i))
				if operation == 'x':
					qc.h(qb[i])
				if operation == 'y':
					qc.u1(-np.pi/2,qb[i])
					qc.h(qb[i])
				indices.append(i)
			

			qc.measure(qb,cb)
			job = qk.execute(qc, backend = qk.Aer.get_backend('qasm_simulator'), shots=self.shots)
			result = job.result()
			self.result = result.get_counts(qc)
			H0 = 0
			for key, value in self.result.items():
				key1 = key[::-1]
				eigenval = 1
				if count == n_qubits:
					eigenval = 1
				else:
					for ind in indices:
						eigenval *= 1 if key1[ind] == '0' else -1
				H0 += eigenval*value
			H0 /= self.shots
			E += H0*circuit.factor.real

		qb = qk.QuantumRegister(n_qubits+1)
		cb = qk.ClassicalRegister(n_qubits+1)
		qc = qk.QuantumCircuit(qb,cb)
		self.wavefunction_ansatz(qc,qb,cb,t,measure=True)
		print(self.result)
		print(E)
		return(E)

	def energy(self,t):
		n_nodes = self.n_fermi
		n_qubits = self.n_qubits
		circuit_list = self.circuit_list
		E = 0
		for circuit in circuit_list:
			qb = qk.QuantumRegister(n_qubits+1)
			cb = qk.ClassicalRegister(n_qubits+1)
			qc = qk.QuantumCircuit(qb,cb)
			self.wavefunction_ansatz(qc,qb,cb,t)
			indices = []
			for i in range(n_qubits):
				operation = circuit.get(i).op
				if operation == '':
					continue
				exec('qc.{}(qb[{}])'.format(operation,i))
				if operation == 'x':
					qc.h(qb[i])
				if operation == 'y':
					qc.u1(-np.pi/2,qb[i])
					qc.h(qb[i])
				indices.append(i)

			qc.measure(qb,cb)
			job = qk.execute(qc, backend = qk.Aer.get_backend('qasm_simulator'), shots=self.shots)
			result = job.result()
			self.result = result.get_counts(qc)
			H0 = 0
			for key, value in self.result.items():
				key1 = key[::-1]
				eigenval = 1
				for ind in indices:
					eigenval *= 1 if key1[ind] == '0' else -1
				H0 += eigenval*value
			H0 /= self.shots
			E += H0*circuit.factor.real

		qb = qk.QuantumRegister(n_qubits+1)
		cb = qk.ClassicalRegister(n_qubits+1)
		qc = qk.QuantumCircuit(qb,cb)
		self.wavefunction_ansatz(qc,qb,cb,t,measure=True)
		print('----------')
		print(self.result)
		print(E)
		print('Complete energy:')
		self.calculate_energy(t,circuit_list2)
		return(E)


	def non_gradient_optimization(self):
		t = np.random.rand(self.n_fermi*self.n_fermi*(self.n_qubits - self.n_fermi)*(self.n_qubits - self.n_fermi))

		res = minimize(self.energy, t, method='COBYLA', options={'disp': True},tol=1e-12)
		print(res.x)
		
		print(self.result)


n_qubits = 4
n_fermi = 2
delta = 1
g = -1
h_pq = np.identity(n_qubits)
for p in range(n_qubits):
    h_pq[p,p] *= delta*(p - (p%2))/2
    
h_pqrs = np.zeros((n_qubits,n_qubits,n_qubits,n_qubits))
for p in range(0,n_qubits-1,2):
    for r in range(0,n_qubits-1,2):
        h_pqrs[p,p+1,r,r+1] = -0.5*g

Pairing = Hamiltonian(n_qubits)

circuit_list = Pairing.get_circuits(h_pq,h_pqrs)
for oplist in circuit_list:
	print(oplist.factor)
	for i in range(n_qubits):
		print('qb[{}]'.format(i),oplist.get(i).op)


print('---------')
Pairing = Hamiltonian(n_qubits)
circuit_list2 = Pairing.get_circuits(h_pq,h_pqrs,remove_identity=False)



test = VQE_UCC(n_qubits=n_qubits,n_fermi=n_fermi,circuit_list=circuit_list)
test.non_gradient_optimization()


						
