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

ry_iC  = np.load('cobyla/coeff_ideal_RY.npy') 
ry_iC = ry_iC[[i for i in range(0,len(ry_iC),5)]]
ucc_iC = np.load('cobyla/coeff_ideal_UCC.npy') 
ucc_iC = ucc_iC[[i for i in range(0,len(ucc_iC),5)]]

fciC = np.load('fci_coeff.npy')
fciC = fci[[i for i in range(0,len(fci),5)]]
x = np.load('bonds.npy')
x2 = x[[i for i in range(0,len(x),5)]]

plt.plot(x2,[fidelity(i,j) for i,j in zip(fciC,ry_iC)],'.',label='RY')
plt.plot(x2,[fidelity(i,j) for i,j in zip(fciC,ucc_iC)],'+',label='UCCD')
plt.legend()
plt.ylim([-0.1,1.1])


plt.show()


ry_i  = np.load('cobyla/E_ideal_RY.npy') 
ry_i = ry_i[[i for i in range(0,len(ry_i),5)]]
ucc_i = np.load('cobyla/E_ideal_UCC.npy') 
ucc_i = ucc_i[[i for i in range(0,len(ucc_i),5)]]

fci = np.load('cc.npy')
fci2 = fci[[i for i in range(0,len(fci),5)]]
hf = np.load('hf.npy')

plt.figure()
plt.fill_between(x,0,0.0016, alpha=0.2,label='Chemical accuracy')
plt.plot(x2,np.abs(ry_i-fci2),'.',label='RY')
plt.plot(x2,np.abs(ucc_i-fci2),'+',label='UCCD')
plt.legend()
plt.title('Absolute error')
plt.xlabel('Bond length [Å]')
plt.ylabel('Energy [Hartree]')

save('tex/h2_error_ideal_1.tex')

plt.show()
