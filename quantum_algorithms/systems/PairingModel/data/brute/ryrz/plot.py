import numpy as np
import matplotlib.pyplot as plt

x = np.load('grid.npy')
params = np.load('parameters.npy')
Es = np.load('values10000.npy')

print(np.min(Es))
Es = Es.reshape(50,50)
plt.imshow(Es)
plt.show()
