import numpy as np
import matplotlib.pyplot as plt

x = np.arange(-2.0,2.0,0.4)

fci = np.load('pair_fci_n4l8.npy')
hf  = np.load('pair_hf_n4l8.npy')
cc = np.load('pair_ccd_n4l8.npy')

ryrz = np.load('8000shots/E_RYRZ_withcx_8000shots.npy')
#ryrz = np.load('8000shots/E_RYRZ_8000shots.npy')
uccsd = np.load('8000shots/E_UCCD_8000shots.npy')


plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.plot(x,cc,label='CCD')
plt.plot(x,ryrz,label='RYRZ')
plt.plot(x,uccsd,label='UCCD')
plt.legend()
plt.show()
