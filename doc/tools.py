def get_num_parameters(n,l):
    S = n*(l-n)
    D = 0
    for i in range(n):
        for j in range(i+1,n):
            for a in range(n,l):
                for b in range(a+1,l):
                    D += 1
    return S, D

# PRINT STATE

from sty import fg, bg, ef, rs
from sympy.utilities.iterables import multiset_permutations

def make_legal(string):
    L = list(multiset_permutations(string))
    L_new = []
    for elem in L:
        L_new.append(''.join(elem))
    return L_new

def get_state_count(R,n,l):
    Hartree = '0'*n + '1'*(l-n)
    legal_states = make_legal(Hartree)#['0011','1100','0110','1001','0101','1010']
    legal = 0
    illegal = 0

    for i in R:
        # Add legal count
        if i in legal_states:
            legal += R[i]
        # Add illegal count and string
        else:
            illegal += R[i]
    return legal,illegal

def print_state(R,n,l,print_all_legal=False):
    # Legal states (6) 0011,1100,0110,1001,0101,1010
    # Illegal states (10)
    Hartree = '0'*n + '1'*(l-n)
    legal_states = make_legal(Hartree)#['0011','1100','0110','1001','0101','1010']
    legal_dict = {}
    for i in legal_states:
        legal_dict[i] = 0
    legal = 0
    illegal = 0
    string_legal = ''
    string_illegal = ''

    for i in R:
        # Add legal count
        if i in legal_states:
            legal_dict[i] = R[i]
            legal += R[i]
        # Add illegal count and string
        else:
            illegal += R[i]
            string_illegal += '{}\u2758{}\u27E9 '.format(R[i],i)
    # Add legal string
    for i in legal_dict:
        if print_all_legal:
            string_legal += '{}\u2758{}\u27E9 '.format(legal_dict[i],i)
        else:
            if legal_dict[i] != 0:
                string_legal += '{}\u2758{}\u27E9 '.format(legal_dict[i],i)
    # PRINT
    legal_count = fg.green + str(legal) + fg.rs
    illegal_count = fg.red + str(illegal) + fg.rs
    print('States after {} measurements!'.format(legal+illegal))
    print('')
    print('- Legal state count:',legal_count)
    print('')
    print('   \u21B3 '+string_legal)
    print('')
    if illegal > 0:
        print('- Illegal state count:',illegal_count)
        print('')
        if string_illegal != '':
            print('   \u21B3 '+string_illegal)
            print('')


# OpenFermion Hamiltonian to circuit_list
def from_OpenFermion(hamiltonian):
    circ = []
    for term in hamiltonian.terms:
        factor = hamiltonian.terms[term]
        action = [factor]
        if len(term) > 0:
            for i in range(len(term)):
                qbit = term[i][0]
                pauli = term[i][1].lower()
                action.append([qbit,pauli])
        circ.append(action)
    return circ


