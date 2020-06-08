import numpy as np
import matplotlib.pyplot as plt

x   = np.load('omegas.npy')
fci = np.load('fci.npy')
hf  = np.load('hf.npy')
cc  = np.load('cc.npy')

ucc = np.load('UCCSDr_E_i.npy')

plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.plot(x,ucc,label='UCCSDr')
plt.legend()

plt.show()

