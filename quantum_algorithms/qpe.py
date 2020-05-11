import numpy as np
import qiskit as qk

from attributes import QuantumComputer
from algoritm import QuantumAlgorithm


class QPE(QuantumAlgorithm):
    """
    Quantum Phase Estimation algorithm. 
    TODO: Compatible with the QuantumCircuit class from QuantumCircuitOptimizer.
    """
    def __init__(self,
                 hamiltonian,
                 ansatz,
                 n_work,
                 Emax,
                 t = 0.5,
                 n = 100,
                 options={}):
        """
        Input:
            hamiltonian  (class) - Hamiltonian with system information
                                   and circuit_list.
            options      (dict)  - 
        """
        super().__init__(options)
        self.hamiltonian = hamiltonian
        self.n_simulation = hamiltonian.l
        self.n_work = n_work
        self.n_qubits = self.n_work+self.n_simulation
        self.circuit_list = hamiltonian.to_circuit_list(ptype='qpe')
        self.ansatz = ansatz

        self.qb_work = qk.QuantumRegister(self.n_work,'qW')
        self.cb_work = qk.ClassicalRegister(self.n_work,'cW')
        self.qb_simulation = qk.QuantumRegister(self.n_simulation,'qS')
        self.cb_simulation = qk.ClassicalRegister(self.n_simulation,'cS')

        self.qc = qk.QuantumCircuit(self.qb_work,self.cb_work,
                                    self.qb_simulation,self.cb_simulation)

        self.Emax = Emax
        self.t = t
        self.n = n
        self.dt = t/n

        #### Setup options
        self.options = options
        # For execution
        self.shots = 1000 if options.get('shots') == None\
                          else options.get('shots')
        self.seed = options.get('seed')
        if self.seed != None:
            from qiskit.aqua import aqua_globals
            aqua_globals.random_seed = self.seed
        # For Backend
        if options.get('backend') == None:
            self.options['backend'] = 'qasm_simulator' 
        self.backend = qk.Aer.get_backend(options['backend'])
        # For noise model and coupling map
        self.noise_model, self.coupling_map,self.basis_gates  = QuantumComputer(options.get('device'),
                                                               options.get('noise_model'),
                                                               options.get('coupling_map'),
                                                               options.get('basis_gates'))
        # GPU accelerated
        if options.get('gpu'):
            from qiskit_qcgpu_provider import QCGPUProvider
            Provider = QCGPUProvider()
            self.backend = Provider.get_backend(options['backend'])
        
    def estimate(self):
        self.prepare_work_register()
        self.qc = self.ansatz(self.qc,self.qb_simulation,self.n_simulation)
        self.evolve_state()
        #print('Circuit depth:',self.qc.depth())
        #print('Ops:',self.qc.count_ops())
        self.inverse_fourier_transform()
        self.measure()

    def prepare_work_register(self):
        for qbit in range(self.n_work):
            self.qc.h(self.qb_work[qbit])

    def evolve_state(self):
        for w in range(self.n_work):
            theta = (2**w)*self.dt
            for n in range(self.n):
                for gate_info in self.circuit_list:
                    gate = gate_info[0]
                    assert not np.iscomplex(gate_info[1])
                    factor = gate_info[1].real
                    qbit = []
                    targ = gate_info[-1]
                    qbit.append(targ+self.n_work) # add n_work to get correct qubit
                    if len(gate_info) == 4:
                        ctrl = gate_info[-2]
                        qbit.insert(0,ctrl+self.n_work) # add n_work to get correct qubit
                    if gate.char.lower() == 'rz':
                        self.qc.crz(theta*factor,
                                     self.qb_work[w],
                                     self.qb_simulation[targ])
                    elif gate.char.lower() == 'ph':
                        # Phase gate as u1*x*u1*x
                        self.qc.cu1(-factor*theta,
                                    self.qb_work[w],
                                    self.qb_simulation[targ])
                        self.qc.x(self.qb_simulation[targ])
                        self.qc.cu1(-factor*theta,
                                    self.qb_work[w],
                                    self.qb_simulation[targ])
                        self.qc.x(self.qb_simulation[targ])
                    else:
                        self.qc.append(gate.get_qiskit(),qbit,[])
                # Insert Emax term
                self.qc.cu1(theta*self.Emax,
                            self.qb_work[w],
                            self.qb_simulation[0])
                self.qc.x(self.qb_simulation[0])
                self.qc.cu1(theta*self.Emax,
                            self.qb_work[w],
                            self.qb_simulation[0])
                self.qc.x(self.qb_simulation[0])

    def inverse_fourier_transform(self,swap=True):
        if swap:
            for i in range(int(self.n_work/2)):
                self.qc.swap(self.qb_work[i],self.qb_work[self.n_work-i-1])
        for i in range(self.n_work):
            for j in range(i):
                self.qc.cu1(-(2*np.pi)/(2**(i+1-j)),self.qb_work[j],self.qb_work[i])
            self.qc.h(self.qb_work[i])

    def bits2phase(self,bitstr):
        bitstr = [int(i) for i in bitstr]
        phase = 0 
        for i,bit in enumerate(bitstr):
            phase += (2**(-i-1))*bit
        return phase

    def sort_results(self):
        x = [] # phase
        y = [] # hits
        psi = []
        for key,val in self.result.items():
            key = key[::-1]
            phi,eigenstate = key.split(' ')
            phi = phi[::-1]
            assert len(phi) == self.n_work
            assert len(eigenstate) == self.n_simulation
            psi.append(eigenstate)
            x.append(self.Emax - 2*np.pi*self.bits2phase(phi)/self.t)
            y.append(val)
        # Sort arrays
        x = np.array(x)
        y = np.array(y)
        idx = np.argsort(x)
        x = x[idx]
        y = y[idx]
        #psi = psi[idx]
        ## Check if same phase measured for different eigenstates
        x_ = []
        y_ = []
        for i,xi in enumerate(x):
            if i > 0:
                if xi == x_[-1]:
                    y_[-1] += y[i]
                else:
                    x_.append(xi)
                    y_.append(y[i])
            else:
                x_.append(xi)
                y_.append(y[i])
        self.x = x_
        self.y = y_
        return self

    def measure(self):
        self.qc.measure(self.qb_work,self.cb_work)
        self.qc.measure(self.qb_simulation,self.cb_simulation)
        job = qk.execute(self.qc, backend = self.backend, shots=self.shots)
        result = job.result().get_counts(self.qc)
        self.result = result
        self.sort_results()

