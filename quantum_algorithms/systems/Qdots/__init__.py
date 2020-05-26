from quantum_systems import ODQD
from quantum_systems.quantum_dots.one_dim.one_dim_potentials import HOPotential
from coupled_cluster import CCSD

def qdots_mat_elems(n,l,omega,grid_length=5,num_grid_points=1001):
    """
    One dimensional quantum dots matrix elements
    
    l     : Number of spin orbitals / number of qubits
    n     : Number of occupied spin orbitals
    omega : Angular frequency 
    """
    # Generate system
    odho = ODQD(n,l,grid_length,num_grid_points)
    odho.setup_system(potential=HOPotential(omega),add_spin=True)

    one_body = odho.h
    one_body[np.absolute(one_body) < 1e-8 ] = 0
    two_body = odho.u
    two_body[np.absolute(two_body) < 1e-8 ] = 0
    return one_body,two_body

