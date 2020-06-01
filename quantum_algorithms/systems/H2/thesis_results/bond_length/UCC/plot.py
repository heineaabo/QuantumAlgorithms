import numpy as np
import matplotlib.pyplot as plt

x   = np.load('bonds.npy')
fci = np.load('fci.npy')
hf  = np.load('hf.npy')
cc  = np.load('cc.npy')
cobE   = np.load('cobyla/E.npy')
cobVar = np.load('cobyla/var.npy')
powE   = np.load('powell/E.npy')
powVar = np.load('powell/var.npy')
nelE   = np.load('nelder/E.npy')
nelVar = np.load('nelder/var.npy')


plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.plot(x,cc,label='CCSD')
plt.errorbar(x,cobE,yerr=cobVar,label='Cobyla')
plt.errorbar(x,powE,yerr=powVar,label='Powell')
plt.errorbar(x,nelE,yerr=nelVar,label='Nelder-Mead')
plt.legend()
plt.show()
