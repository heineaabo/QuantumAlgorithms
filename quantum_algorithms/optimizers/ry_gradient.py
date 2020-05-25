import qiskit as qk
import numpy as np

class RyGradient:
    def __init__(self,
                 max_iters=200,
                 max_evals=200,
                 step=1, # Step length
                 tol=1e-08):
        self.max_iters = max_iters
        self.max_evals = max_evals
        self.step = step
        self.tol = tol
        self.vqe = None

    def __call__(self,theta):
        #self.energies = []
        qc,qb,qa,cb,ca = None,None,None,None,None
        new_theta = theta
        E_last = 1000
        for i in range(self.max_iters):
            for d_j in range(len(theta)):
                im = 0
                qb = qk.QuantumRegister(self.vqe.n_qubits)
                cb = qk.ClassicalRegister(self.vqe.n_qubits)
                qa = qk.QuantumRegister(1)
                ca = qk.ClassicalRegister(1)
                qc_ansatz = qk.QuantumCircuit(qb,cb)
                qc_ansatz.add_register(qa,ca)
                qc_ansatz = self.vqe.ansatz(theta,qc_ansatz,qb,qa)
                for pauli_string in self.vqe.circuit_list:
                    qc = qk.QuantumCircuit(qb,cb)
                    # Prepare ancilla qubit with Hadamard
                    qc.add_register(qa,ca)
                    qc.h(qa[0])
                    # Regular
                    qc = pauli_string.prepare(qc,qb,qa)
                    qc.h(qa[0])
                    qc.rx(-np.pi/2,qa[0])
                    qc = qc_ansatz + qc
                    im += self.ancilla_measure(pauli_string.factor,qc,qa,ca)
                new_theta[d_j] -= im*self.step
                if len(theta) > 1:
                    print('Updating theta {} -> {} - {} = {}'.format(d_j,theta[d_j],im,new_theta[d_j]))
            E_new = self.vqe.expval(new_theta)
            if len(theta) > 1:
                print('Iteration {} finished!'.format(i),new_theta)
                print('<E> =',E_new)
            else:
                print('<E> = {}, Grad = {}'.format(E_new,im*self.step))
            theta = new_theta
            if np.abs(E_new-E_last) < self.tol:
                print('Converged at iteration {} with <E> ={}'.format(i,E_new))
                break
            E_last = E_new
        return theta

    def ancilla_measure(self,factor,qc,qa,ca):
        qc.measure(qa,ca)
        job = qk.execute(qc, 
                        backend = self.vqe.backend, 
                        shots=self.vqe.shots,
                        seed_transpiler=self.vqe.seed,
                        seed_simulator=self.vqe.seed)
        result = job.result().get_counts(qc)
        im = 0
        for key, value in result.items():
            if key[0] == '0':
                im += 1 - 2*(value/self.vqe.shots)
            elif key[0] == '1':
                im += 2*(value/self.vqe.shots) - 1
        return factor*im

    def set_vqe(self,vqe):
        self.vqe = vqe
