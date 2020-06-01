
#def get_H2(R):
#    from numpy import load,linspace,where
#    if not round(R,3) in linspace(0.1,3.0,59):
#        raise ValueError('Choose allowed bond length!')
#    h_pq = load('molecule/matrix_elements/h_R{}.npy'.format(R))
#    h_pqrs = load('molecule/matrix_elements/v_R{}.npy'.format(R))
#    Enuc = load('molecule/nuclear_repulsion/R{}.npy'.format(R))
#
#    #FCI 
#    bonds = load('molecule/energies/bonds.npy')
#    ind = where(bonds==R)[0][0]
#    fci = load('molecule/energies/fci.npy')[ind]
#    return h_pq,h_pqrs,Enuc,fci


from openfermion.hamiltonians import MolecularData
from openfermionpsi4 import run_psi4
def get_h2_matrix(R,
        basis='sto-3g',
        multiplicity=1,
        charge=0):
    geometry = [['H',[0,0,0]],
                ['H',[0,0,R]]]
    h2_molecule = MolecularData(geometry,basis,multiplicity,charge)
    h2_molecule = run_psi4(h2_molecule,
                            run_mp2=True,
                            #run_hf=True,
                            run_ccsd=True,
                            run_fci=True)
    one_body = h2_molecule.one_body_integrals
    two_body = h2_molecule.two_body_integrals
    Enuc = h2_molecule.nuclear_repulsion
    
    energies = {}
    energies['fci'] = h2_molecule.fci_energy
    energies['hf'] = h2_molecule.hf_energy
    energies['ccsd'] = h2_molecule.ccsd_energy

    return one_body, two_body, Enuc, energies
