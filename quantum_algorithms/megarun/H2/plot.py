import numpy as np
import matplotlib.pyplot as plt
from tikzplotlib import save

ry_i  = np.load('cobyla/E_ideal_RY.npy') 
ry_i = ry_i[[i for i in range(0,len(ry_i),5)]]
ry_n  = np.load('cobyla/E_noisy_RY.npy') 
ucc_i = np.load('cobyla/E_ideal_UCC.npy') 
ucc_i = ucc_i[[i for i in range(0,len(ucc_i),5)]]
ucc_n = np.load('cobyla/E_noisy_UCC.npy') 
uccd_i = np.load('cobyla/E_ideal_UCC_full.npy') 
uccd_i = uccd_i[[i for i in range(0,len(uccd_i),5)]]
uccd_n = np.load('cobyla/E_noisy_UCC_full.npy') 

x = np.load('bonds.npy')
x2 = x[[i for i in range(0,len(x),5)]]
fci = np.load('cc.npy')
fci2 = fci[[i for i in range(0,len(fci),5)]]
hf = np.load('hf.npy')


plt.figure()
plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.plot(x2,ry_i,'.',label='RY')
plt.plot(x2,ucc_i,'+',label='UCCD')
#plt.plot(x,ry_i,label='RY')
#plt.plot(x,ucc_i,label='UCCD')
#plt.plot(x,uccd_i,label='UCCD full')
plt.legend()
plt.title('Ideal simulation')
plt.xlabel('Bond length [Å]')
plt.ylabel('Energy [Hartree]')

save('tex/h2_bond_ideal_1.tex')

plt.figure()
plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.plot(x,ry_n,label='RY')
plt.plot(x,ucc_n,label='UCCD')
#plt.plot(x,uccd_n,label='UCCD full')
plt.legend()
plt.title('Noisy simulation')
plt.xlabel('Bond length [Å]')
plt.ylabel('Energy [Hartree]')

save('tex/h2_bond_noisy_1.tex')

plt.figure()
plt.fill_between(x,0,0.0016, alpha=0.2,label='Chemical accuracy')
plt.plot(x2,np.abs(ry_i-fci2),'.',label='RY')
plt.plot(x2,np.abs(ucc_i-fci2),'+',label='UCCD')
#plt.plot(x,np.abs(uccd_i-fci),label='UCCD full')
plt.legend()
plt.title('Absolute error')
plt.xlabel('Bond length [Å]')
plt.ylabel('Energy [Hartree]')

save('tex/h2_error_ideal_1.tex')

plt.show()
