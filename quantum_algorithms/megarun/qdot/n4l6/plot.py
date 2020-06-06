import numpy as np
import matplotlib.pyplot as plt

x = np.load('omegas.npy')
fci = np.load('fci.npy')
hf  = np.load('hf.npy')
cc  = np.load('cc.npy')

ryrz = np.load('E_ideal_RYRZ.npy')
ucc = np.load('E_ideal_UCCSD.npy')

plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.plot(x,cc,label='CCSD')
plt.plot(x,ryrz,label='RYRZ')
plt.plot(x,ucc,label='UCCSD')
plt.legend()

plt.show()
