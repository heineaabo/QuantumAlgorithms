import numpy as np
def get_pairing_matrix(n,l,delta,g):
    """
    Pairing model matrix elements
    l     : Number of spin orbitals / number of qubits
    n     : Number of occupied spin orbitals
    delta : Level spacing
    g     : Interaction strength
    """
    h = np.identity(l)
    for p in range(l):
            h[p,p] *= delta*(p - (p%2))/2
            
    v = np.zeros((l,l,l,l))
    for p in range(0,l-1,2):
            for r in range(0,l-1,2):
                    v[p,p+1,r,r+1] = -0.5*g
    return h,v
