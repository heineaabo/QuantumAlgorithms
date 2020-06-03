import numpy as np
import matplotlib.pyplot as plt
from tikzplotlib import save

ry_i  = np.load('cobyla/E_ideal_RY.npy') 
ry_n  = np.load('cobyla/E_noisy_RY.npy') 
ucc_i = np.load('cobyla/E_ideal_UCC.npy') 
ucc_n = np.load('cobyla/E_noisy_UCC.npy') 
uccd_i = np.load('cobyla/E_ideal_UCC_full.npy') 
uccd_n = np.load('cobyla/E_noisy_UCC_full.npy') 
uccsd_i = np.load('cobyla/E_ideal_UCCSD.npy') 
uccsd_n = np.load('cobyla/E_noisy_UCCSD.npy') 

x = np.load('bonds.npy')
fci = np.load('cc.npy')
hf = np.load('hf.npy')


plt.figure()
plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.plot(x,ry_i,label='RY')
plt.plot(x,ucc_i,label='UCCD')
#plt.plot(x,uccd_i,label='UCCD full')
plt.plot(x,uccsd_i,label='UCCSD')
plt.legend()
plt.title('Ideal simulation')
plt.xlabel('Bond length [Å]')
plt.ylabel('Energy [Hartree]')

save('tex/h2_bond.tex')

plt.figure()
plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.plot(x,ry_n,label='RY')
plt.plot(x,ucc_n,label='UCCD')
#plt.plot(x,uccd_n,label='UCCD full')
plt.plot(x,uccsd_n,label='UCCSD')
plt.legend()
plt.title('Noisy simulation')
plt.xlabel('Bond length [Å]')
plt.ylabel('Energy [Hartree]')

save('tex/h2_bond_noisy.tex')

plt.figure()
plt.fill_between(x,0,0.0016, alpha=0.2,label='Chemical accuracy')
plt.plot(x,np.abs(ry_i-fci),label='RY')
plt.plot(x,np.abs(ucc_i-fci),label='UCCD')
#plt.plot(x,np.abs(uccd_i-fci),label='UCCD full')
plt.plot(x,np.abs(uccsd_i-fci),label='UCCSD')
plt.legend()
plt.title('Absolute error')
plt.xlabel('Bond length [Å]')
plt.ylabel('Energy [Hartree]')

save('tex/h2_error.tex')

plt.show()
