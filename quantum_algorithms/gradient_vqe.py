from vqe import VQE

def optimize_gradient(self,
                      theta,
                      max_iters=200,
                      max_evals=200,
                      tol=1e-08):
    self.energies = []
    qc,qb,,qa,cb,ca = None,None,None,None,None
    new_theta = theta
    for i in range(max_iters):
        qa = QuantumRegister(1) # Ancilla qubit
        ca = ClassicalRegister(1)
        qc.add_register(qa)
        E_inter = self.expval(new_theta)
        print('<E> =',E_inter)
        for d_j in range(len(theta)):
            im = 0
            for pauli_string in self.circuit_list:
                factor = pauli_string[0].real
                qubit_list = []
                qb,qc,cb =self.initialize_circuit()
                qc = self.ansatz(theta,qc,qb,qa,d_j)
                for qubit,gate in pauli_string[1:]:
                    qc = self.add_gate(qubit,gate,qc,qb,qa)
                    qubit_list.append(qubit)
                im += self.ancilla_measure(factor,qc,qa,ca)
            print('gradient theta{}: {}'.format(d_j,im))
            new_theta[d_j] -= im 
        print('Iteration {} finished!'.format(i))
        print(len(theta))
    return theta
VQE.optimize_gradient = optimize_gradient

def ancilla_measure(self,factor,qc,qa,ca):
    qc.measure(qa,ca)
    job = qk.execute(qc, backend = qk.Aer.get_backend('qasm_simulator',
                     seed=self.seed), 
                     shots=self.shots)
    result = job.result().get_counts(qc)
    im = 0
    for key, value in result.items():
        if key == '0':
            im += 1 - 2*(value/self.shots)
        elif key == '1':
            im += 2*(value/self.shots) - 1
    return factor*im
VQE.ancilla_measure = ancilla_measure
