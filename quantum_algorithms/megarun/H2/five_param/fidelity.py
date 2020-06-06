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

ryrz_i = np.load('RYRZ_good_coeff.npy') 
uccsd_i = np.load('UCCSDr_good_coeff.npy') 

fci = np.load('../fci_coeff.npy')
fci = fci[[i for i in range(0,len(fci),5)]]
x = np.load('../bonds.npy')
x = x[[i for i in range(0,len(x),5)]]

for a,b in zip(ryrz_i,uccsd_i):
    print(a)
    print(b)
    print('')


plt.plot(x,[fidelity(i,j) for i,j in zip(fci,uccsd_i)],label='UCCSD')
plt.plot(x,[fidelity(i,j) for i,j in zip(fci,ryrz_i)],label='RYRZ')
plt.legend()


#save('tex/fidelity.tex')

plt.show()
