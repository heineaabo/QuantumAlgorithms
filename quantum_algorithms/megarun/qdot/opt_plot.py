import numpy as np
import matplotlib.pyplot as plt
from tikzplotlib import save


# RY ansatz
cob = np.load('cobyla/E_conv_w1_RYPAIRING.npy')
powl= np.load('powell/E_conv_w1_RYPAIRING.npy')
nel = np.load('nelder-mead/E_conv_w1_RYPAIRING.npy')
spsa= np.load('SPSA/E_conv_w1_RYPAIRING.npy')
plt.figure()
plt.plot(cob,label='Cobyla')
plt.plot(powl,label='Powell')
plt.plot(nel,label='Nelder-Mead')
plt.plot(spsa,label='SPSA')
plt.xlabel('Function evaluations')
plt.ylabel('Energy [Hartree]')
plt.title('UCCD ansatz')
plt.legend()

# UCCD ansatz
cob = np.load('cobyla/E_conv_w1_UCCD.npy')
powl= np.load('powell/E_conv_w1_UCCD.npy')
nel = np.load('nelder-mead/E_conv_w1_UCCD.npy')
spsa= np.load('SPSA/E_conv_w1_UCCD.npy')
plt.figure()
plt.plot(cob,label='Cobyla')
plt.plot(powl,label='Powell')
plt.plot(nel,label='Nelder-Mead')
plt.plot(spsa,label='SPSA')
plt.xlabel('Function evaluations')
plt.ylabel('Energy [Hartree]')
plt.title('UCCD ansatz')
plt.legend()

# UCCDr ansatz
cob = np.load('cobyla/E_conv_w1_UCCDr.npy')
powl= np.load('powell/E_conv_w1_UCCDr.npy')
nel = np.load('nelder-mead/E_conv_w1_UCCDr.npy')
spsa= np.load('SPSA/E_conv_w1_UCCDr.npy')
plt.figure()
plt.plot(cob,label='Cobyla')
plt.plot(powl,label='Powell')
plt.plot(nel,label='Nelder-Mead')
plt.plot(spsa,label='SPSA')
plt.xlabel('Function evaluations')
plt.ylabel('Energy [Hartree]')
plt.title('UCCDr ansatz')
plt.legend()

# UCCSD ansatz
cob = np.load('cobyla/E_conv_w1_UCCSD.npy')
powl= np.load('powell/E_conv_w1_UCCSD.npy')
nel = np.load('nelder-mead/E_conv_w1_UCCSD.npy')
spsa= np.load('SPSA/E_conv_w1_UCCSD.npy')
plt.figure()
plt.plot(cob,label='Cobyla')
plt.plot(powl,label='Powell')
plt.plot(nel,label='Nelder-Mead')
plt.plot(spsa,label='SPSA')
plt.xlabel('Function evaluations')
plt.ylabel('Energy [Hartree]')
plt.title('UCCSD ansatz')
plt.legend()

# UCCSDr ansatz
cob = np.load('cobyla/E_conv_w1_UCCSDr.npy')
powl= np.load('powell/E_conv_w1_UCCSDr.npy')
nel = np.load('nelder-mead/E_conv_w1_UCCSDr.npy')
spsa= np.load('SPSA/E_conv_w1_UCCSDr.npy')
plt.figure()
plt.plot(cob,label='Cobyla')
plt.plot(powl,label='Powell')
plt.plot(nel,label='Nelder-Mead')
plt.plot(spsa,label='SPSA')
plt.xlabel('Function evaluations')
plt.ylabel('Energy [Hartree]')
plt.title('UCCSDr ansatz')
plt.legend()


plt.show()
