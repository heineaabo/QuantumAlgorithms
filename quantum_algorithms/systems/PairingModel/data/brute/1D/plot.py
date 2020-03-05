import numpy as np
import matplotlib.pyplot as plt


x = np.load('x.npy')
vals100k = np.load('values100k.npy')
vals10k = np.load('values10k.npy')
vals1000 = np.load('values1000.npy')
#vals500 = np.load('values500.npy')

legal1000 = np.load('legal1000.npy')
illegal1000 = np.load('illegal1000.npy')
tot = legal1000+illegal1000
legal1000 /= tot
illegal1000 /= tot

minx1000 = np.argmin(vals1000)
minval1000 = vals1000[minx1000]
#minx500 = np.argmin(vals500)
#minval500 = vals500[minx500]
minx10k = np.argmin(vals10k)
minval10k = vals10k[minx10k]
minx100k = np.argmin(vals100k)
minval100k = vals100k[minx100k]

print('100000 function evaluations:')
print('Minimum value at theta = {}\n<E> = {}'.format(x[minx100k],minval100k))
print('10000 function evaluations:')
print('Minimum value at theta = {}\n<E> = {}'.format(x[minx10k],minval10k))
print('1000 function evaluations:')
print('Minimum value at theta = {}\n<E> = {}'.format(x[minx1000],minval1000))
#print('500 function evaluations:')
#print('Minimum value at theta = {}\n<E> = {}'.format(x[minx500],minval500))

plt.figure()
plt.plot(x,vals100k,'-',label='100000 evaluations')
plt.plot(x,vals10k,'-',label='10000 evaluations')
plt.plot(x,vals1000,'-',label='1000 evaluations')
#plt.plot(x,vals500,label='500 evaluations')
plt.legend()


plt.figure()
plt.plot(x,legal1000,label='Legal')
plt.plot(x,illegal1000,label='Ilegal')
plt.legend()

plt.show()
