import numpy as np
import matplotlib.pyplot as plt

x = np.arange(-2,2,0.05)

uccE = np.load('london_UCCD_1000shots.npy')
ryE = np.load('london_RYPAIRING_1000shots.npy')
ry2E = np.load('london_RYPAIRING_10000shots.npy')
fci = np.load('fci.npy')

uccE1 = uccE[:,0]
uccE2 = uccE[:,1]
uccE3 = uccE[:,2]

ryE1 = ryE[:,0]
ryE2 = ryE[:,1]
ryE3 = ryE[:,2]

ry2E1 = ry2E[:,0]
ry2E2 = ry2E[:,1]
ry2E3 = ry2E[:,2]

plt.figure()
plt.plot(x,fci,label='FCI')
plt.plot(x,uccE1,label='UCCD ideal')
plt.plot(x,uccE2,label='UCCD noisy_n')
plt.plot(x,uccE3,label='UCCD noisy_w')
plt.title('UCCD N=1000')
plt.legend()

plt.figure()
plt.plot(x,fci,label='FCI')
plt.plot(x,ryE1,label='RY ideal')
plt.plot(x,ryE2,label='RY noisy_n')
plt.plot(x,ryE3,label='RY noisy_w')
plt.title('RY N=1000')
plt.legend()

plt.figure()
plt.plot(x,fci,label='FCI')
plt.plot(x,ry2E1,label='RY ideal')
plt.plot(x,ry2E2,label='RY noisy_n')
plt.plot(x,ry2E3,label='RY noisy_w')
plt.title('RY N=10000')
plt.legend()

plt.figure()
plt.plot(x,np.abs(ryE1-fci),label='RY N=1e03')
plt.plot(x,np.abs(ry2E1-fci),label='RY N=1e04')
plt.plot(x,np.abs(uccE1-fci),label='UCCD')
plt.legend()

plt.show()
