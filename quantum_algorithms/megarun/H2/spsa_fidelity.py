import numpy as np
import matplotlib.pyplot as plt
from tikzplotlib import save

def fidelity(v,w):
    # v -> FCI
    # w -> VQE
    c = 0
    for i,j in zip(v,w):
        c += np.abs(i)*np.sqrt(j)
    return c**2

uccsd_i = np.load('SPSA/coeff_ideal_UCCSD.npy') 
uccsdr_i = np.load('SPSA/coeff_ideal_UCCSDr.npy') 
ryrz_i = np.load('SPSA/coeff_ideal_RYRZ.npy') 

fci = np.load('fci_coeff.npy')
x = np.load('bonds.npy')

#plt.plot(x,[fidelity(i,j) for i,j in zip(fci,uccsd_i)],label='UCCSD')
#plt.plot(x,[fidelity(i,j) for i,j in zip(fci,uccsdr_i)],label='UCCSDr')
plt.plot(x,[fidelity(i,j) for i,j in zip(fci,ryrz_i)],label='RYRZ')
plt.legend()

save('tex/spsa_fidelity.tex')

plt.show()

#k = 0
#for i,j in zip(ry_i,fci):
#    print('For R =',x[k])
#    print(i)
#    print(j)
#    print('')
#    k += 1
