import pickle

direct = '/Users/heine2307/Documents/Universitet/UiO/Master/GitHub/VQE/quantum_algorithms/attributes/'

def QuantumComputer(device,noise_model,coupling_map,basis_gates):
        name = device
        if device == None:
            return None,None,None
        else:
            if noise_model:
                noise_model = pickle.load(open(direct+'noise_models/'+device+'.pkl','rb'))
            if coupling_map:
                coupling_map = pickle.load(open(direct+'coupling_maps/'+device+'.pkl','rb'))
            if basis_gates:
                basis_gates  = pickle.load(open(direct+'basis_gates/'+device+'.pkl','rb'))
            return noise_model,coupling_map,basis_gates
