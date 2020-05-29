import numpy as np
import matplotlib.pyplot as plt

normal = np.load('ESSEX_GC_normal.npy')
noise_no = np.load('ESSEX_GC_noisy_n.npy')
noise_w = np.load('ESSEX_GC_noisy_w.npy')

x = np.arange(1,2*np.pi+1,2*np.pi/500)
plt.figure(0)
plt.axhline(y=-0.618, color='r', linestyle='--',label='FCI')
plt.plot(x,normal,label='ideal')
plt.plot(x,noise_no,label='noisy (no meas_fit)')
plt.plot(x,noise_w,label='noisy (with meas_fit)')
plt.legend()
plt.xlabel('Theta')
plt.ylabel('Energy [Hartree]')
plt.title('Pairing model n=2, l=4 (N=1000) LONDON')

plt.show()
