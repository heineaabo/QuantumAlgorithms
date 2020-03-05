import numpy as np
import matplotlib.pyplot as plt

x = np.load('shots.npy')
values = np.load('values.npy')
legal = np.load('legal.npy')
illegal = np.load('illegal.npy')

total = legal+illegal
legal /= total
illegal /= total

#plt.plot(x,values)
plt.plot(x,legal,label='Legal')
plt.plot(x,illegal,label='Illegal')
plt.legend()
plt.show()
