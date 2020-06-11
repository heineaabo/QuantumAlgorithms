import numpy as np
import matplotlib.pyplot as plt
from tikzplotlib import save

strs = ['030','040','050','060','074','100','130','160','190','220','250','280']
Rs = [0.3,0.4,0.5,0.6,0.74,1.0,1.3,1.6,1.9,2.2,2.5,2.8]

hf = np.load('hf.npy')
fci = np.load('fci.npy')


means = np.zeros_like(Rs)

for i,s in enumerate(strs):
    E = np.load('E_{}.npy'.format(s))
    theta = np.load('t_{}.npy'.format(s))
    #print('Energies    Theta')
    #for e,t in zip(E,theta):
    #    print('{:<.7f}  {:<.7f}'.format(e,t))
    #print('R = {}:'.format(Rs[i]))
    #print('FCI         HF          |Mean        Var')
    #print('{:<.7f}  {:<.7f}  |{:<.7f}  {:<.7f}'.format(fci[i],hf[i],np.mean(E[-5:]),np.var(E[-5:])))
    means[i] = np.mean(E[-5:])
    if s == '130':
        means[i] = np.mean(E[15:20])
    if s == '250':
        means[i] = np.mean(E[10:17])

fci_real = np.load('../../../../../megarun/H2/cc.npy')
hf_real = np.load('../../../../../megarun/H2/hf.npy')
x_real = np.load('../../../../../megarun/H2/bonds.npy')
x2 = [0,5,10,15,22,35,50,65,80,95,110,125]
fci2 = fci_real[x2]
plt.figure()
plt.plot(x_real,fci_real,label='FCI')
plt.plot(x_real,hf_real,label='HF')
plt.plot(Rs,means,'.',label='VQE')
plt.legend()
save('h2_london.tex')

plt.figure()
plt.plot(Rs,np.abs(fci2-means),'.',label='VQE')
plt.legend()


plt.show()




