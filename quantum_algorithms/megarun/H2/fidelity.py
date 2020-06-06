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

ry_i  = np.load('cobyla/coeff_ideal_RY.npy') 
ry_i = ry_i[[i for i in range(0,len(ry_i),5)]]
#ry_n  = np.load('cobyla/coeff_noisy_RY.npy') 
#ryrz_i = np.load('cobyla/coeff_ideal_RYRZ.npy') 
ucc_i = np.load('cobyla/coeff_ideal_UCC.npy') 
ucc_i = ucc_i[[i for i in range(0,len(ucc_i),5)]]
#ucc_n = np.load('cobyla/coeff_noisy_UCC.npy') 
#uccd_i = np.load('cobyla/coeff_ideal_UCC_full.npy') 
#uccd_n = np.load('cobyla/coeff_noisy_UCC_full.npy') 
#uccsd_i = np.load('cobyla/coeff_ideal_UCCSD.npy') 
#uccsd_n = np.load('cobyla/coeff_noisy_UCCSD.npy') 

fci = np.load('fci_coeff.npy')
fci = fci[[i for i in range(0,len(fci),5)]]
x = np.load('bonds.npy')
x = x[[i for i in range(0,len(x),5)]]

plt.plot(x,[fidelity(i,j) for i,j in zip(fci,ry_i)],'.',label='RY')
plt.plot(x,[fidelity(i,j) for i,j in zip(fci,ucc_i)],'+',label='UCCD')
#plt.plot(x,[fidelity(i,j) for i,j in zip(fci,uccd_i)],label='UCCD full')
#plt.plot(x,[fidelity(i,j) for i,j in zip(fci,uccsd_i)],label='UCCSD')
#plt.plot(x,[fidelity(i,j) for i,j in zip(fci,ryrz_i)],label='RYRZ')
plt.legend()
plt.ylim([-0.1,1.1])

save('tex/fidelity.tex')

plt.show()

#k = 0
#for i,j in zip(ry_i,fci):
#    print('For R =',x[k])
#    print(i)
#    print(j)
#    print('')
#    k += 1
