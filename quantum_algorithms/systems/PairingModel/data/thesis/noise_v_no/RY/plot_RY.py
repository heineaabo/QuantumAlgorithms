import numpy as np
import matplotlib.pyplot as plt
import tikzplotlib

normal = np.load('normal.npy')
noise_no = np.load('noisy_no_meas_fit.npy')
noise_w = np.load('noisy_with_meas_fit.npy')

normal10 = np.load('normal_10k.npy')
noise_no10 = np.load('noisy_no_meas_fit_10k.npy')
noise_w10 = np.load('noisy_with_meas_fit_10k.npy')

x = np.arange(0,2*np.pi,2*np.pi/500)
plt.figure(0)
plt.axhline(y=-0.618, color='r', linestyle='--',label='FCI')
plt.plot(x,normal,label='N=1000')
#plt.plot(x,normal,label='ideal')
#plt.plot(x,noise_no,label='noisy (no meas_fit)')
#plt.plot(x,noise_w,label='noisy (with meas_fit)')
#plt.legend()
#plt.xlabel('Theta')
#plt.ylabel('Energy [Hartree]')
#plt.title('Pairing model n=2, l=4 (N=1000) LONDON')
#
#x = np.arange(1,2*np.pi+1,2*np.pi/500)
#plt.figure(1)
#plt.axhline(y=-0.618, color='r', linestyle='--',label='FCI')
plt.plot(x,normal10,label='N=10000')
#plt.plot(x,normal10,label='ideal')
#plt.plot(x,noise_no10,label='noisy (no meas_fit)')
#plt.plot(x,noise_w10,label='noisy (with meas_fit)')
plt.legend()
plt.xlabel('Theta')
plt.ylabel('Energy [Hartree]')
plt.title('Pairing model with n=2 and l=4')
#plt.title('Pairing model n=2, l=4 (N=10000) LONDON')
plt.xlim([3.6,6.1])
plt.ylim([-0.8,0.4])

tikzplotlib.save('shots.tex')

plt.show()
