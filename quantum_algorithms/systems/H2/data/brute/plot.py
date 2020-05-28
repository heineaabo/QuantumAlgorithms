import numpy as np
import matplotlib.pyplot as plt

grid = np.load('grid.npy')
x = np.load('parameters.npy')
y = np.load('values1000.npy')
print('Min: x =',min(y))

plt.imshow(y.reshape((len(grid),len(grid))),extent=[grid[0],grid[-1],grid[0],grid[-1]])
plt.show()
