import pickle

def QuantumComputer(device,noise_model,coupling_map):
        name = device
        if device == None:
            return None,None
        else:
            if noise_model:
                noise_model = pickle.load(open('noise_models/'+device+'.pkl','r'))
            if coupling_map:
                coupling_map = pickle.load(open('coupling_maps/'+device+'.pkl','r'))
            return noise_model,coupling_map
