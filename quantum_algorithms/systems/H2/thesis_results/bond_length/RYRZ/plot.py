import numpy as np
import matplotlib.pyplot as plt

x = np.load('bonds.npy')
fci = np.load('fci.npy')
hf = np.load('hf.npy')
E = np.load('E.npy')
var = np.load('var.npy')
print(var)


plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
plt.errorbar(x,E,yerr=var,fmt='-o',label='VQE')
plt.show()
