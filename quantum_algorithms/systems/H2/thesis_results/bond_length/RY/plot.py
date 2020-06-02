import numpy as np
import matplotlib.pyplot as plt

E_i = np.load('cobyla/E_ideal.npy')
E_n = np.load('cobyla/E_noisy.npy')

fci = np.load('fci.npy')
hf = np.load('hf.npy')

x = np.load('bonds.npy')

plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.plot(x,E_i,label='VQE ideal')
plt.plot(x,E_n,label='VQE noisy')
plt.legend()
plt.show()
