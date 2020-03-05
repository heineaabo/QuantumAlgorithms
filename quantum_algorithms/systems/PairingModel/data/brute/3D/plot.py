import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

x = np.load('grid.npy')
params = np.load('parameters.npy')
vals1000 = np.load('values1000.npy')

minx1000 = np.argmin(vals1000)
minval1000 = vals1000[minx1000]
print('1000 function evaluations:')
print('Minimum value at theta = {}\n<E> = {}'.format(x[minx1000%len(x)],minval1000))

n = len(x)
ims = vals1000.reshape(n,n,n)

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = plt.plot([], [], 'ro')


im = plt.imshow(ims[:,:,0],extent=(0,2*np.pi,0,2*np.pi))
plt.title(str(0))
cbar = fig.colorbar(im)
#def init():
#    im = ims[:,:,0]
#    return [im]

def update(i):
    im.set_array(ims[:,:,i])
    plt.title(str(x[i]))
    return [im]

ani = FuncAnimation(fig, update, frames=np.arange(n)) #,blit=True
plt.show()


