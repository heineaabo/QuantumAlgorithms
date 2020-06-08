import numpy as np
import matplotlib.pyplot as plt

E = np.load('ibmq_energies.npy')
t = np.load('ibmq_theta.npy')

plt.plot(t,E)
plt.show()
