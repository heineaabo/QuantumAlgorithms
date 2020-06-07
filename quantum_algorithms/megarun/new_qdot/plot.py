import numpy as np
import matplotlib.pyplot as plt

x = np.load('omegas.npy')
fci = np.load('fci.npy')
hf = np.load('hf.npy')

ryrz_i  = np.load('RYRZ_E_i.npy')
ryrz_n  = np.load('RYRZ_E_n.npy')
uccdr_i = np.load('UCCDr_E_i.npy')
uccdr_n = np.load('UCCDr_E_n.npy')
uccsd_i = np.load('UCCSD_E_i.npy')
uccsd_n = np.load('UCCSD_E_n.npy')


plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.plot(x,ryrz_i,'.',label='RYRZ')
plt.plot(x,uccdr_i,'+',label='UCCDr')
plt.plot(x,uccsd_i,'*',label='UCCSD')
plt.legend()

plt.show()
