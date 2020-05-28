
def get_H2(R):
    from numpy import load,linspace,where
    if not round(R,3) in linspace(0.1,3.0,59):
        raise ValueError('Choose allowed bond length!')
    h_pq = load('molecule/matrix_elements/h_R{}.npy'.format(R))
    h_pqrs = load('molecule/matrix_elements/v_R{}.npy'.format(R))
    Enuc = load('molecule/nuclear_repulsion/R{}.npy'.format(R))

    #FCI 
    bonds = load('molecule/energies/bonds.npy')
    ind = where(bonds==R)[0][0]
    fci = load('molecule/energies/fci.npy')[ind]
    return h_pq,h_pqrs,Enuc,fci
