import numpy as np
import matplotlib.pyplot as plt

x = np.load('grid.npy')
params1000 = np.load('parameters.npy')
vals1000 = np.load('values1000.npy')


minx1000 = np.argmin(vals1000)
minval1000 = vals1000[minx1000]
print('1000 function evaluations:')
print('Minimum value at theta = {}\n<E> = {}'.format(x[minx1000],minval1000))

im = vals1000.reshape((len(x),len(x)))

plt.imshow(im,extent=(0,2*np.pi,0,2*np.pi))
plt.show()
