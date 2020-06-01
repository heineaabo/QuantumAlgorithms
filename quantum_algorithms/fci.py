import scipy.special as special
import numpy as np
import itertools as it

def FCI(n,l,h,v,ret_all=False):
    """
    n - Number of particles
    l - Number of spin orbitals
    h - One-body matrix elements
    v - Two-body matrix elements
    """
    n_SD = int(special.binom(l,n))
    H = np.zeros((n_SD,n_SD),dtype=complex)
    # Get occupied indices for all states
    S = configurations(n,l) 
    for row,bra in enumerate(S):
        for col,ket in enumerate(S):
            if np.sum(np.equal(bra,ket)) == bra.shape:
                # One-body contributions (Assumed diagonal)
                for p in bra:
                    H[row,col] += h[p,p]
                # Two-body contributions
                for p,q in it.combinations(bra,2):
                    H[row,col] += v[p,q,p,q]
            else:
                # Two-body contributions
                to_create = list(set(bra)-set(ket))
                to_annihilate = list(set(ket)-set(bra))
                eq = [i for i in bra if i not in to_create]
                # One orbital different
                if len(to_create) == 1 and len(to_annihilate) == 1:
                    for p in eq:
                        for c,a in zip(to_create,to_annihilate):
                            H[row,col] += v[p,c,p,a]
                # Two orbitals different
                elif len(to_create) == 2 and len(to_annihilate) == 2:
                    p,q = to_create
                    r,s = to_annihilate
                    H[row,col] += v[p,q,r,s]
    Es,Vs = np.linalg.eigh(H)
    if ret_all:
        return Es,Vs
    else:
        return Es[0]

def configurations(n,l):
    states = []
    for state in it.combinations(range(l),n):
        states.append([orb for orb in state])
    return np.asarray(states)

