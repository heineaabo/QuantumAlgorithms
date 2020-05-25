import numpy as np
import matplotlib.pyplot as plt

x = np.load('ry_grid.npy')
params = np.load('ry_parameters.npy')
Es = np.load('ry_values1000.npy')

print(np.min(Es))
plt.plot(x,Es)
plt.show()
