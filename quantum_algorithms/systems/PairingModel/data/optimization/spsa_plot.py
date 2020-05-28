import numpy as np
import matplotlib.pyplot as plt

spsa   = np.load('spsa.npy') 
spsa_n = np.load('spsa_noisy.npy') 

avgs = [1,2,5,10]

#for i,avg in enumerate(avgs):
#    plt.figure(i)
#    plt.plot(spsa[i],label='spsa')
#    plt.plot(spsa_n[i],label='spsa_n')
#    plt.legend()
#    plt.title('Avgs: '+str(avg))
plt.figure(0)
for i,avg in enumerate(avgs):
    plt.plot(spsa[i],label=str(avg))
    plt.legend()
    plt.title('SPSA')
plt.figure(1)
for i,avg in enumerate(avgs):
    plt.plot(spsa_n[i],label=str(avg))
    plt.legend()
    plt.title('Noisy')
plt.show()
