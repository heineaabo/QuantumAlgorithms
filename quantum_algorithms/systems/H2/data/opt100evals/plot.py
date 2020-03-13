import numpy as np
import matplotlib.pyplot as plt

#nlm= np.load('Nelder-Mead500.npy')
#nlm_avg= np.load('Nelder-Mead500_avg10.npy')
#nlm_std= np.load('Nelder-Mead500_std10.npy')
#pwl = np.load('Powell500.npy')
#pwl_avg = np.load('Powell500_avg10.npy')
#pwl_std = np.load('Powell500_std10.npy')
cob = np.load('Cobyla500.npy')
cob_avg = np.load('Cobyla500_avg10.npy')
cob_std = np.load('Cobyla500_std10.npy')
spsa = np.load('SPSA500.npy')
spsa_avg = np.load('SPSA500_avg10.npy')
spsa_std = np.load('SPSA500_std10.npy')

hf   = np.load('hf.npy') 
ccsd = np.load('ccsd.npy')  
fci  = np.load('fci.npy')  

#plt.plot(nlm,'b-',label='Nelder-Mead')
#plt.plot(nlm_avg,'b--')
#plt.plot(pwl,'y-',label='Powell')
#plt.plot(pwl_avg,'y--')
plt.plot(cob,'g-',label='Cobyla')
plt.plot(cob_avg,'g--')
plt.plot(spsa,'r-',label='SPSA')
plt.plot(spsa_avg,'r--')

x = np.arange(max([len(i) for i in [cob,spsa]]))

plt.plot(hf,label='Hartree-Fock')
plt.plot(ccsd,label='CCSD')
plt.plot(fci,label='FCI')

plt.legend()
plt.show()

