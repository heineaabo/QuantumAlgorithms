import numpy as np
import qiskit as qk
import sys



from ucc import *
from tools import print_state, get_num_parameters


n = int(sys.argv[1])
l = int(sys.argv[2])

qb = qk.QuantumRegister(l)
cb = qk.ClassicalRegister(l)
qc = qk.QuantumCircuit(qb,cb)

depth = 1

ansatz = UnitaryCoupledCluster(n,l,'SD',depth)
theta = ansatz.new_parameters()

qc = ansatz(theta,qc,qb)
print(qc.draw())#,interactive=True))

qc.measure(qb,cb)
job = qk.execute(qc, 
                backend = qk.Aer.get_backend('qasm_simulator'), 
                shots=500)
result = job.result().get_counts(qc)
E = 0
print_state(result,n,l)
