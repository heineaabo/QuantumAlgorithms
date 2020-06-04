import numpy as np
import matplotlib.pyplot as plt
from tikzplotlib import save
import sys

save_plot = False
if len(sys.argv) > 1:
    if sys.argv[1].lower() == 'save':
        save_plot = True

ry_i  = np.load('cobyla/E_ideal_RY.npy') 
ry_n  = np.load('cobyla/E_noisy_RY.npy') 
ucc_i = np.load('cobyla/E_ideal_UCC.npy') 
ucc_n = np.load('cobyla/E_noisy_UCC.npy') 
uccd_i = np.load('cobyla/E_ideal_UCC_full.npy') 
uccd_n = np.load('cobyla/E_noisy_UCC_full.npy') 
uccsd_i = np.load('cobyla/E_ideal_UCCSD.npy') 
uccsd_n = np.load('cobyla/E_noisy_UCCSD.npy') 

x = np.load('omegas.npy')
#fci = np.load('cc.npy')
fci = np.load('fci.npy')
hf = np.load('hf.npy')


plt.figure()
plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.plot(x,ry_i,label='RY')
#plt.plot(x,ucc_i,label='UCCD')
plt.plot(x,uccd_i,label='UCCD full')
plt.plot(x,uccsd_i,label='UCCSD')
plt.legend()
plt.title('Ideal simulation')
plt.xlabel('Omega')
plt.ylabel('Energy [Hartree]')

if save_plot:
    save('tex/ideal.tex')

plt.figure()
plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.plot(x,ry_n,label='RY')
#plt.plot(x,ucc_n,label='UCCD')
plt.plot(x,uccd_n,label='UCCD full')
plt.plot(x,uccsd_n,label='UCCSD')
plt.legend()
plt.title('Noisy simulation')
plt.xlabel('Omega')
plt.ylabel('Energy [Hartree]')

if save_plot:
    save('tex/noisy.tex')

plt.figure()
plt.fill_between(x,0,0.0016, alpha=0.2,label='Chemical accuracy')
plt.plot(x,np.abs(fci-ry_i),label='RY')
#plt.plot(x,np.abs(fci-ucc_i),label='UCCD')
plt.plot(x,np.abs(fci-uccd_i),label='UCCD full')
plt.plot(x,np.abs(fci-uccsd_i),label='UCCSD')
plt.legend()
plt.title('Absolute error')
plt.xlabel('Omega')
plt.ylabel('Energy [Hartree]')

if save_plot:
    save('tex/error.tex')

plt.show()
