import numpy as np
import matplotlib.pyplot as plt

london_n = np.load('noisy_no_meas_fit.npy')
london_w = np.load('noisy_with_meas_fit.npy')

rome_n= np.load('ROME_noisy_no_meas_fit.npy')
rome_w = np.load('ROME_noisy_with_meas_fit.npy')

essex_n = np.load('ESSEX_noisy_no_meas_fit.npy')
essex_w = np.load('ESSEX_noisy_with_meas_fit.npy')

x = np.arange(0,2*np.pi,2*np.pi/500)
plt.figure(0)
plt.axhline(y=-0.618, color='r', linestyle='--',label='FCI')
plt.plot(x,london_n,label='London')
plt.plot(x,rome_n,label='Rome')
plt.plot(x,essex_n,label='Essex')
plt.legend()
plt.xlabel('Theta')
plt.ylabel('Energy [Hartree]')
plt.title('Without meas_fitter')

x = np.arange(1,2*np.pi+1,2*np.pi/500)
plt.figure(1)
plt.axhline(y=-0.618, color='r', linestyle='--',label='FCI')
plt.plot(x,london_w,label='London')
plt.plot(x,rome_w,label='Rome')
plt.plot(x,essex_w,label='Essex')
plt.legend()
plt.xlabel('Theta')
plt.ylabel('Energy [Hartree]')
plt.title('With meas_fitter')

plt.show()
