import numpy as np
import matplotlib.pyplot as plt

ry_i  = np.load('cobyla/E_ideal_RY.npy') 
ry_n  = np.load('cobyla/E_noisy_RY.npy') 
ucc_i = np.load('cobyla/E_ideal_UCC.npy') 
ucc_n = np.load('cobyla/E_noisy_UCC.npy') 

x = np.load('omegas.npy')
#fci = np.load('cc.npy')
fci = np.load('fci.npy')
hf = np.load('hf.npy')


plt.figure()
plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.plot(x,ry_i,label='RY')
plt.plot(x,ucc_i,label='UCCD')
plt.legend()
plt.title('Ideal simulation')
plt.xlabel('Omega')
plt.ylabel('Energy [Hartree]')

plt.figure()
plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.plot(x,ry_n,label='RY')
plt.plot(x,ucc_n,label='UCCD')
plt.legend()
plt.title('Noisy simulation')
plt.xlabel('Omega')
plt.ylabel('Energy [Hartree]')

plt.figure()
plt.fill_between(x,0,0.0016, alpha=0.2)
plt.plot(x,np.abs(fci-ry_i),label='RY')
plt.plot(x,np.abs(fci-ucc_i),label='UCCD')
plt.legend()
plt.title('Absolute error')
plt.xlabel('Omega')
plt.ylabel('Energy [Hartree]')

plt.show()
