from attributes import QuantumComputer
import qiskit as qk


class QuantumAlgorithm:
    def __init__(self,options):
        """
        Quantum Algorithm superclass. Handles execution specifications.

            options:
                * seed         (int)  - Seed for simulator and transpiler.
                * shots        (int)  - Number of execution shots.
                * backend      (str)  - Name of backend.
                * device       (str)  - Name of IBMQ quantum computer.
                * noise_model  (bool) - Noise model of device.
                * coupling_map (bool) - Coupling map of device.
                * print        (bool) - Print for every function evaluation.
                * count_states (bool) - Count legal and illegal states for
                                        all function evaluations.
        """
        #### Setup options
        self.options = options
        # For execution
        self.shots = 1000 if options.get('shots') == None\
                          else options.get('shots')
        self.seed = options.get('seed')
        if self.seed != None:
            from qiskit.aqua import aqua_globals
            aqua_globals.random_seed = self.seed
        self.prnt = options.get('print')
        self.ancilla_measure = options.get('ancilla') if options.get('ancilla') != None else False

        # For Backend
        if options.get('backend') == None:
            self.options['backend'] = 'qasm_simulator' 
        self.backend = qk.Aer.get_backend(options['backend'])
        # For noise model, coupling map and basis gates
        if options.get('device') == None:
            self.noise_model, self.coupling_map, self.basis_gates = None,None,None
        else:
            device = QuantumComputer(options.get('device'))
            if options.get('noise_model') != None:
                self.noise_model = device.noise_model
                self.coupling_map = device.coupling_map
                self.basis_gates = device.basis_gates
        #self.noise_model, self.coupling_map, self.basis_gates  = QuantumComputer(options.get('device'),options.get('noise_model'),options.get('coupling_map'),options.get('basis_gates'))

        # GPU accelerated
        if options.get('gpu'):
            from qiskit_qcgpu_provider import QCGPUProvider
            Provider = QCGPUProvider()
            self.backend = Provider.get_backend(options['backend'])

    def measure(self,qc,qb,cb):
        """
        Get measurement results of circuit.
        NOTE: qb and cb can be lists, if multiple registers are to be measured.
        """
        if isinstance(qb,list) and isinstance(cb,list):
            assert len(qb) == len(cb)
            for q,c in zip(qb,cb):
                qc.measure(q,c)
        else:
            qc.measure(qb,cb)
        job = qk.execute(qc, 
                        backend = self.backend, 
                        shots=self.shots,
                        noise_model=self.noise_model,
                        coupling_map=self.coupling_map,
                        basis_gates=self.basis_gates,
                        seed_transpiler=self.seed,
                        seed_simulator=self.seed)
        return job.result().get_counts(qc)
