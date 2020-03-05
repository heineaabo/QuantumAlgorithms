import numpy as np
import matplotlib.pyplot as plt


methods = ['Powell','Nelder-Mead','Cobyla','L_BFGS_B',]

powell = np.load('data/Powell.npy')
cobyla = np.load('data/Cobyla.npy')
L_BFGS_B = np.load('data/L-BFGS-B.npy')
nelder = np.load('data/Nelder-Mead.npy')
spsa = np.load('data/SPSA.npy')

x_powell = np.arange(len(powell))
x_cobyla = np.arange(len(cobyla))
x_lbfgsb = np.arange(len(L_BFGS_B))
x_nelder= np.arange(len(nelder))
x_spsa= np.arange(len(spsa))

plt.plot(x_powell,powell,label='Powell')
plt.plot(x_cobyla,cobyla,label='Cobyla')
plt.plot(x_lbfgsb,L_BFGS_B,label='L-BFGS-B')
plt.plot(x_nelder,nelder,label='Nelder-Mead')
plt.plot(x_spsa,spsa,label='SPSA')
plt.xlabel('Function evaluations')
plt.ylabel('Energy [Hartrees]')
plt.legend()
plt.show()
