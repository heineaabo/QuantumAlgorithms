from .attributes import QuantumComputer
import qiskit as qk


class QuantumAlgorithm:
    def __init__(self,l,options):
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

        self.ibmq = False
        if options.get('ibmq') == True:
            print('Running on real quantum computer')
            self.ibmq = True
            self.backend = options['backend']
            from qiskit.tools.monitor import job_monitor
            self.monitor = job_monitor
            from attributes import get_measurement_fitter
            self.meas_fitter = get_measurement_fitter(l,
                                                     self.backend,
                                                     None,
                                                     self.shots)
            
        else:
            # For Backend
            if options.get('backend') == None:
                self.options['backend'] = 'qasm_simulator' 
            self.backend = qk.Aer.get_backend(options['backend'])
            # For noise model, coupling map and basis gates
            self.noise_model, self.coupling_map, self.basis_gates = None,None,None
            self.meas_fitter = None
            if options.get('device') != None:
                device = QuantumComputer(options.get('device'))
                if options.get('noise_model') != None:
                    self.noise_model = device.noise_model
                    # Create error mitigation fitter
                    if options.get('meas_fit') in [None,True]:
                        from attributes import get_measurement_fitter
                        self.meas_fitter = get_measurement_fitter(l,
                                                                 self.backend,
                                                                 device,
                                                                 self.shots)
                if options.get('coupling_map') != None:
                    self.coupling_map = device.coupling_map
                if options.get('basis_gates') != None:
                    self.basis_gates = device.basis_gates
            # Qubit layout, virtual to physical
        self.layout = options.get('layout')
        # Optimization level
        self.optimization_level= 1 if options.get('optimization_level')==None else options['optimization_level']

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
        if self.ibmq:
            job = qk.execute(qc, 
                            backend = self.backend, 
                            shots=self.shots,
                            optimization_level=self.optimization_level,
                            initial_layout=self.layout)
            self.monitor(job)
        else:
            job = qk.execute(qc, 
                            backend = self.backend, 
                            shots=self.shots,
                            optimization_level=self.optimization_level,
                            noise_model=self.noise_model,
                            coupling_map=self.coupling_map,
                            basis_gates=self.basis_gates,
                            initial_layout=self.layout,
                            seed_transpiler=self.seed,
                            seed_simulator=self.seed)
        return job.result().get_counts(qc)
