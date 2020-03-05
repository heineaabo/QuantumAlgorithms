import numpy as np
import matplotlib.pyplot as plt

hf   = np.load('hf.npy') 
cisd = np.load('cisd.npy') 
ccsd = np.load('ccsd.npy')  
fci  = np.load('fci.npy')  

powell = np.load('Powell.npy') 
nelder = np.load('Nelder-Mead.npy')
cobyla = np.load('Cobyla.npy')
spsa   = np.load('SPSA.npy') 

x = np.linspace(0.5,1.5,11)

#plt.plot(x,hf-fci,label='Hartree-Fock')
#plt.plot(x,cisd,label='CISD')
plt.plot(x,ccsd-fci,label='CCSD')
#plt.plot(x,fci,label='FCI')
#plt.plot(x,powell,label='Powell')
plt.plot(x,nelder-fci,label='Nelder-Mead')
plt.plot(x,cobyla-fci,label='Cobyla')
plt.plot(x,spsa-fci,label='SPSA')

plt.legend()
plt.xlabel('Bond length [Å]')
plt.ylabel('Energy Difference [Hartree]')
plt.show()


