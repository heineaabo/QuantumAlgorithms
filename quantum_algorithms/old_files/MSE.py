import qiskit as qk
import numpy as np
from scipy.optimize import minimize
from Qoperator import *
qk.IBMQ.load_accounts()
						


class VQE_MaxCut:
	def __init__(self,w=None,X=None,y=None,shots=100,n_nodes = None):
		self.shots = shots
		if not w is None:
			self.n_nodes = w.shape[0]
		self.w = w
		self.X = X
		self.y = y
		if not X is None:
			self.n_nodes = X.shape[1]
		self.preds = []
		self.loss = []

	def wavefunction_ansatz(self,qc,qb,theta,X = None,ansatz = None):
		if ansatz is None:
			for qubit in range(self.n_nodes):
				qc.ry(theta[qubit],qb[qubit])
					
			for qubit in range(self.n_nodes-1):
				qc.cx(qb[qubit],qb[qubit+1])
		if ansatz == 'mse':
			theta  = theta.reshape(X.shape[0],X.shape[0])


			for qubit in range(self.n_nodes):

				qc.ry(np.dot(X,theta[qubit,:]),qb[qubit])


					
			for qubit in range(self.n_nodes-1):
				qc.cx(qb[qubit],qb[qubit+1])

			

	def show_ansatz(self,theta):
		n_nodes = self.n_nodes
		qb = qk.QuantumRegister(n_nodes)
		cb = qk.ClassicalRegister(n_nodes)
		qc = qk.QuantumCircuit(qb,cb)
		for qubit in range(self.n_nodes):
			qc.ry(theta[qubit],qb[qubit])
					
		for qubit in range(self.n_nodes-1):
			qc.cx(qb[qubit],qb[qubit+1])
		qc.measure(qb,cb)
		job = qk.execute(qc, backend = qk.Aer.get_backend('qasm_simulator'), shots=self.shots)
		result = job.result()
		result = result.get_counts(qc)
		return(result)

	def energy(self,theta):
		n_nodes = self.n_nodes
		E = 0
		for i in range(n_nodes):
			for j in range(i+1,n_nodes):
				qb = qk.QuantumRegister(n_nodes)
				cb = qk.ClassicalRegister(n_nodes)
				qc = qk.QuantumCircuit(qb,cb)
				

				self.wavefunction_ansatz(qc,qb,theta)
				
				qc.z(qb[i])
				qc.z(qb[j])

				
				qc.measure(qb,cb)
				job = qk.execute(qc, backend = qk.Aer.get_backend('qasm_simulator'), shots=self.shots)
				result = job.result()
				self.result = result.get_counts(qc)
				E_terms = 0
				for key, value in self.result.items():
					key1 = key[::-1]
					e1 = 1 if key1[i] == '0' else -1
					e2 = 1 if key1[j] == '0' else -1
					eigenval = e1*e2
					E_terms += eigenval*value
				E_terms /= self.shots
				E += E_terms*self.w[i,j]
		print(E)
		print(self.show_ansatz(theta))

		return(E)

	def MSE(self,theta):
		n_nodes = self.n_nodes
		E = 0
		y = self.y
		X = self.X
		for i in range(y.shape[0]):
			term1 = 0
			
			term2 = 0
			for j in range(n_nodes):
				qb = qk.QuantumRegister(n_nodes)
				cb = qk.ClassicalRegister(n_nodes)
				qc = qk.QuantumCircuit(qb,cb)


				self.wavefunction_ansatz(qc,qb,theta,X[i,:],ansatz='mse')

				qc.z(qb[j])

				qc.measure(qb,cb)
				job = qk.execute(qc, backend = qk.Aer.get_backend('qasm_simulator'), shots=self.shots)
				result = job.result()
				self.result = result.get_counts(qc)
				E_terms = 0
				for key, value in self.result.items():
					key1 = key[::-1]
					eigenval = 1 if key1[j] == '0' else -1
					E_terms += eigenval*value
				E_terms /= self.shots
				term1 += E_terms*(1/(2**(j+1)))
				


				for k in range(n_nodes):
					qb = qk.QuantumRegister(n_nodes)
					cb = qk.ClassicalRegister(n_nodes)
					qc = qk.QuantumCircuit(qb,cb)
					

					self.wavefunction_ansatz(qc,qb,theta,X[i,:],ansatz='mse')
					
					qc.z(qb[j])
					qc.z(qb[k])

					
					qc.measure(qb,cb)
					job = qk.execute(qc, backend = qk.Aer.get_backend('qasm_simulator'), shots=self.shots)
					result = job.result()
					self.result = result.get_counts(qc)
					E_term = 0
					for key, value in self.result.items():
						
						key1 = key[::-1]
						
						e1 = 1 if key1[j] == '0' else -1
						e2 = 1 if key1[k] == '0' else -1
						eigenval = e1*e2
						E_terms += eigenval*value
					E_terms /= self.shots
					term2 += E_terms*(0.25/(2**(j+k+2)))
			E += (y[i] - 0.5*(1 - 2**(-self.n_nodes)))*term1  + term2

		print('Energy: ',E)
		loss = 0
		y_pred = np.zeros(y.shape[0])
		for i in range(y.shape[0]):
			y_pred[i] = self.predict(X[i,:],theta)
			print('Pred:',y_pred[i],'Actual:',y[i],'rel_err:',np.abs(y_pred[i] - y[i])/y[i])

		
		print('Loss:',np.mean((y - y_pred)**2))
		self.loss.append(np.mean((y - y_pred)**2))
		self.preds.append(y_pred)
		print('Theta:',theta)
		print('------------')
		return(E)


	def predict(self,X,theta):
		n_nodes = self.n_nodes
		qb = qk.QuantumRegister(n_nodes)
		cb = qk.ClassicalRegister(n_nodes)
		qc = qk.QuantumCircuit(qb,cb)
		
		
		theta  = theta.reshape(X.shape[0],X.shape[0])
		for qubit in range(self.n_nodes):

			qc.ry(np.dot(X,theta[qubit,:]),qb[qubit])


				
		for qubit in range(self.n_nodes-1):
			qc.cx(qb[qubit],qb[qubit+1])

		


		qc.measure(qb,cb)
		job = qk.execute(qc, backend = qk.Aer.get_backend('qasm_simulator'), shots=self.shots)
		result = job.result()
		self.result = result.get_counts(qc)
		"""
		key = max(self.result, key=self.result.get)
		key = key[::-1]
		y_pred = 0
		for i,bit in enumerate(key):
			y_pred += int(bit)/2**(i+1)

		"""
		y_pred = 0
		for key,value in self.result.items():
			key1 = key[::-1]
			pred = 0
			for i,bit in enumerate(key1):
				pred += int(bit)/(2**(i+1))
			y_pred += pred*value
		y_pred /= self.shots
		return(y_pred)


	def non_gradient_optimization(self):
		thetas = []
		E = []
		for n_thetas in range(10):
			theta = 2*np.pi*np.random.randn(self.n_nodes)
			thetas.append(theta)
			E.append(self.energy(theta))
		E = np.array(E)
		thetas = np.array(thetas)
		theta_arg = np.argmin(E)
		theta = thetas[theta_arg]
		print(E)

		res = minimize(self.energy, theta, method='powell', options={'disp': True},tol=1e-12)
		print(res.x)
		print(self.show_ansatz(res.x))
		print(self.result)

	def non_gradient_optimization_mse(self):
		theta = 2*np.pi*np.random.randn(self.n_nodes*self.n_nodes)
		print('theta shape', theta.shape)
		print(theta)
		res = minimize(self.MSE, theta, method='L-BFGS-B', options={'disp': True},tol=1e-12)
		print(res.x)
		print(self.energy(res.x))
		print(self.result)

	def gradient_optimization(self,theta):
		n_nodes = self.n_nodes + 1 #Add ancilla qubit
		theta_grad = np.zeros(theta.shape[0])
		for idx,theta_i in enumerate(theta):
			E = 0
			for i in range(n_nodes-1):
				for j in range(i+1,n_nodes-1):
					qb = qk.QuantumRegister(n_nodes)
					cb = qk.ClassicalRegister(n_nodes)
					qc = qk.QuantumCircuit(qb,cb)
					qc.h(qb[n_nodes-1])

					for qubit in range(n_nodes-1):
						
						if idx == qubit:
							qc.cy(qb[n_nodes-1],qb[qubit])
						qc.ry(theta_i,qb[qubit])

					
					for qubit in range(n_nodes-2):
						qc.cx(qb[qubit],qb[qubit+1])
						
					qc.cz(qb[n_nodes-1],qb[i])
					qc.cz(qb[n_nodes-1],qb[j])

					qc.h(qb[n_nodes-1])
					qc.rx(np.pi/2,qb[n_nodes-1])

					qc.measure(qb,cb)
					job = qk.execute(qc, backend = qk.Aer.get_backend('qasm_simulator'), shots=self.shots)
					result = job.result()
					self.result = result.get_counts(qc)
					E_terms = 0
					shots = 0
					for key, value in self.result.items():
						key1 = key[::-1]
						if key1[-1] == '0':
							shots += value
					
					prob = shots/self.shots
					
					prob = -2*prob + 1

					E += prob*self.w[i,j]
			theta_grad[idx] = E
		return(theta_grad)

	def gradient_optimization_mse(self,theta):
		X = self.X
		y = self.y
		print(X.shape)
		n_nodes = self.n_nodes + 1 #Add ancilla qubit
		theta_grad = np.zeros(theta.shape[0])
		for idx,theta_i in enumerate(theta):
			for y_idx in range(y.shape[0]):
				E = 0
				term1 = 0
				for i in range(n_nodes-1):
					qb = qk.QuantumRegister(n_nodes)
					cb = qk.ClassicalRegister(n_nodes)
					qc = qk.QuantumCircuit(qb,cb)

					qc.h(qb[n_nodes-1])

					for qubit in range(n_nodes-1):
						qc.ry(X[y_idx,qubit],qb[qubit])
						
					for qubit in range(n_nodes-2):
						qc.cx(qb[qubit],qb[qubit+1])

					for qubit in range(n_nodes-1):
						
						if idx == qubit:
							qc.cy(qb[n_nodes-1],qb[qubit])
						qc.ry(theta_i,qb[qubit])
					
					for qubit in range(n_nodes-2):
						qc.cx(qb[qubit],qb[qubit+1])

					qc.cz(qb[n_nodes-1],qb[i])

					qc.h(qb[n_nodes-1])
					qc.rx(np.pi/2,qb[n_nodes-1])


					qc.measure(qb,cb)
					job = qk.execute(qc, backend = qk.Aer.get_backend('qasm_simulator'), shots=self.shots)
					result = job.result()
					self.result = result.get_counts(qc)
					E_terms = 0
					shots = 0
					for key, value in self.result.items():
						key1 = key[::-1]
						if key1[-1] == '0':
							shots += value
					
					prob = shots/self.shots
					
					prob = -2*prob + 1

					term1 += prob*(1/(2**(i+1)))

					for j in range(n_nodes-1):
						qb = qk.QuantumRegister(n_nodes)
						cb = qk.ClassicalRegister(n_nodes)
						qc = qk.QuantumCircuit(qb,cb)
						qc.h(qb[n_nodes-1])

						for qubit in range(n_nodes-1):
							qc.ry(X[y_idx,qubit],qb[qubit])
							
							if idx == qubit:
								qc.cy(qb[n_nodes-1],qb[qubit])
							qc.ry(theta_i,qb[qubit])

						
						for qubit in range(n_nodes-2):
							qc.cx(qb[qubit],qb[qubit+1])
							
						qc.cz(qb[n_nodes-1],qb[i])
						qc.cz(qb[n_nodes-1],qb[j])

						qc.h(qb[n_nodes-1])
						qc.rx(np.pi/2,qb[n_nodes-1])

						qc.measure(qb,cb)
						job = qk.execute(qc, backend = qk.Aer.get_backend('qasm_simulator'), shots=self.shots)
						result = job.result()
						self.result = result.get_counts(qc)
						E_terms = 0
						shots = 0
						for key, value in self.result.items():
							key1 = key[::-1]
							if key1[-1] == '0':
								shots += value
						
						prob = shots/self.shots
						
						prob = -2*prob + 1

						E += prob*(0.25/(2**(i+j+2)))
				theta_grad[idx] += (y[y_idx] - 0.5*(1 - 2**(-(n_nodes-1))))*term1 + E
			
		return(theta_grad)


	def gradient_descent(self,iters=1000,learning_rate=1e-5):
		theta =  2*np.pi*np.random.randn(self.n_nodes)

		grad2 = np.zeros_like(theta)
		for i in range(iters):
			grads = self.gradient_optimization(theta)
			
			theta = theta - learning_rate*grads
			print(self.energy(theta))
			print(self.result)
			print('      ')
			print('--------')
			grad2 = grads

		return(theta)

	def gradient_descent_mse(self,iters=1000,learning_rate=1e-1):
		theta =  2*np.pi*np.random.randn(self.n_nodes)
		for i in range(iters):
			grads = self.gradient_optimization_mse(theta)
			print(grads)
			y_pred = np.zeros(self.y.shape[0])
			theta = theta - learning_rate*grads
			for j in range(self.y.shape[0]):
				y_pred[j] = self.predict(self.X[j,:],theta)
				print('Pred:',y_pred[j],'Actual:',self.y[j],'rel_err:',np.abs(y_pred[j] - self.y[j])/self.y[j])

		return(theta)

	def prediction_print(self,theta,X,y):
	    y_pred = np.zeros(y.shape[0])
	    print('{:<4} {:<12} {:<12} {:<12}'.format('i','Prediction','Actual','Rel_error'))
	    print('{:-^42s}'.format('-'))
	    for i in range(y.shape[0]):
	        y_pred[i] = self.predict(X[i,:],theta)
	        print('{:<4} {:<12f} {:<12f} {:<12f}'.format(i,y_pred[i],y[i],np.abs(y_pred[i] - y[i])/y[i]))
	        #print('Pred:',y_pred[i],'Actual:',y[i],'rel_err:',np.abs(y_pred[i] - y[i])/y[i])
	    print('Loss:',np.mean((y - y_pred)**2))
	    #self.loss.append(np.mean((y - y_pred)**2))
	    #self.preds.append(y_pred)
	    print('Theta:',theta)
	    print('{:-^42s}'.format('-'))

#w = np.array([[0,3,0,3,1],[3,0,2,0,0],[0,2,0,2,0],[3,0,2,0,3],[1,0,0,3,0]])

#w = np.array([[0,2,1,2],[2,0,3,2],[1,3,0,1],[2,2,1,0]])
w = np.array([[0,1,0,2],[1,0,4,3],[0,4,0,3],[2,3,3,0]])

test = VQE_MaxCut(w,shots=1000)
test.non_gradient_optimization()
print(test.gradient_descent(1000,1e-2))

"""
np.random.seed(42)

x = np.pi*np.random.uniform(size=(2,4))


y = 10*x[:,0] + 2*x[:,1] + 3*x[:,2] + 0.5*x[:,3]

print(y.shape)

ymin = np.min(y)
ymax = np.max(y)

y = y/(ymax + 10)

print(np.max(y))



test = VQE_MaxCut(X=x,y=y,shots=1000)
test.non_gradient_optimization_mse()
#test.gradient_descent_mse(iters=1000,learning_rate=1e-1)
"""






