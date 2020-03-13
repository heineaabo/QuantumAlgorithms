
def get_H2(R):
    from numpy import load,linspace
    if not round(R,3) in linspace(0.1,3.0,59):
        raise ValueError('Choose allowed bond length!')
    h_pq = np.load('molecule/matrix_elements/h_{}.npy'.format(R))
    h_pqrs = np.load('/moleculematrix_elements/v_{}.npy'.format(R))
    Enuc = np.load('molecule/nuclear_repulsion/{}.npy'.format(R))
    return h_pq,h_pqrs,Enuc
