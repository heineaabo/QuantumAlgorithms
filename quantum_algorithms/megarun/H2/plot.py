import numpy as np
import matplotlib.pyplot as plt

ry_i  = np.load('cobyla/E_ideal_RY.npy') 
ry_n  = np.load('cobyla/E_noisy_RY.npy') 
ucc_i = np.load('cobyla/E_ideal_UCC.npy') 
ucc_n = np.load('cobyla/E_noisy_UCC.npy') 

x = np.load('bonds.npy')
fci = np.load('fci.npy')
hf = np.load('hf.npy')


plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.plot(x,ry_i,label='RY')
plt.plot(x,ucc_i,label='UCCD')
plt.legend()
plt.title('Ideal simulation')
plt.xlabel('Bond length [Å]')
plt.ylabel('Energy [Hartree]')

plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.plot(x,ry_n,label='RY')
plt.plot(x,ucc_n,label='UCCD')
plt.legend()
plt.title('Noisy simulation')
plt.xlabel('Bond length [Å]')
plt.ylabel('Energy [Hartree]')

plt.show()
