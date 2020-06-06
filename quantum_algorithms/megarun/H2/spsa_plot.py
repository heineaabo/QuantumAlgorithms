import numpy as np
import matplotlib.pyplot as plt
from tikzplotlib import save

uccsd_i = np.load('SPSA/E_ideal_UCCSD.npy') 
#uccd_n = np.load('SPSA/E_noisy_UCCD.npy') 
uccsdr_i = np.load('SPSA/E_ideal_UCCSDr.npy') 
#uccsd_n = np.load('SPSA/E_noisy_UCCSD.npy') 
ryrz = np.load('SPSA/E_ideal_RYRZ.npy') 

x = np.load('bonds.npy')
fci = np.load('cc.npy')
hf = np.load('hf.npy')


plt.figure()
plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
#plt.plot(x,uccsd_i,label='UCCSD')
#plt.plot(x,uccsdr_i,label='UCCSDr')
plt.plot(x,ryrz,label='RYRZ')
plt.legend()
plt.title('Ideal simulation')
plt.xlabel('Bond length [Å]')
plt.ylabel('Energy [Hartree]')

save('tex/h2_bond_ideal_5.tex')

#plt.figure()
#plt.plot(x,fci,label='FCI')
#plt.plot(x,hf,label='HF')
#plt.plot(x,uccd_n,label='UCCD full')
#plt.plot(x,uccsd_n,label='UCCSD')
#plt.legend()
#plt.title('Noisy simulation')
#plt.xlabel('Bond length [Å]')
#plt.ylabel('Energy [Hartree]')
#
#save('tex/h2_bond_noisy_5.tex')

plt.figure()
plt.fill_between(x,0,0.0016, alpha=0.2,label='Chemical accuracy')
#plt.plot(x,np.abs(uccsd_i-fci),label='UCCSD')
#plt.plot(x,np.abs(uccsdr_i-fci),label='UCCSDr')
plt.plot(x,np.abs(ryrz-fci),label='RYRZ')
plt.legend()
plt.title('Absolute error')
plt.xlabel('Bond length [Å]')
plt.ylabel('Energy [Hartree]')

save('tex/h2_error_5.tex')

plt.show()
