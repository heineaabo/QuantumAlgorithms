import pickle

#direct = '/Users/heine2307/Documents/Universitet/UiO/Master/GitHub/VQE/quantum_algorithms/attributes/'
direct = '/home/heineaabo/Documents/UiO/Master/VQE/quantum_algorithms/attributes/'

#def QuantumComputer(device,noise_model,coupling_map,basis_gates=False):
#        name = device
#        if device == None:
#            return None,None,None
#        else:
#            if noise_model:
#                noise_model = pickle.load(open(direct+'noise_models/'+device+'.pkl','rb'))
#            if coupling_map:
#                coupling_map = pickle.load(open(direct+'coupling_maps/'+device+'.pkl','rb'))
#            if basis_gates:
#                basis_gates  = pickle.load(open(direct+'basis_gates/'+device+'.pkl','rb'))
#            return noise_model,coupling_map,basis_gates


class QuantumComputer:
    def __init__(self,name):
        self.name = name
        self.noise_model = pickle.load(open(direct+'noise_models/'+name+'.pkl','rb'))
        self.coupling_map = pickle.load(open(direct+'coupling_maps/'+name+'.pkl','rb'))
        self.basis_gates  = pickle.load(open(direct+'basis_gates/'+name+'.pkl','rb'))
