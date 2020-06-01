import numpy as np
from quantum_systems import ODQD
from quantum_systems.quantum_dots.one_dim.one_dim_potentials import HOPotential
from coupled_cluster import CCSD

def get_qdot_matrix(n,l,omega,grid_length=5,n_grid=1001):
    odho = ODQD(n,l,grid_length,n_grid)
    odho.setup_system(potential=HOPotential(omega),add_spin=True)

    one_body = odho.h
    one_body[np.absolute(one_body) < 1e-8 ] = 0
    two_body = odho.u
    two_body[np.absolute(two_body) < 1e-8 ] = 0


    # Coupled Cluster
    Eref = odho.compute_reference_energy()
    ccsd = CCSD(odho,verbose=False)
    ccsd.compute_ground_state()
    Ecc = ccsd.compute_energy()
    return one_body,two_body,[Eref,Ecc]

