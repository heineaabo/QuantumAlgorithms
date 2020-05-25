import numpy as np
import matplotlib.pyplot as plt
import sys
#plt.cm.get_cmap('tab20')

if len(sys.argv) > 1:
    prefer = [sys.argv[i] for i in range(1,len(sys.argv))]
else:
    prefer = None

data = np.load('QGD.npy')
data_noisy = np.load('QGD_noisy.npy')
methods = ['GD','Adam','Adagrad','RMSprop']
powell = np.load('Powell.npy')
nelder = np.load('Nelder-Mead.npy')
powell_n = np.load('Powell_n.npy')
nelder_n = np.load('Nelder-Mead_n.npy')

for i,method in enumerate(methods):
    if prefer == None or method in prefer:
        plt.plot(np.arange(len(data[i])*2),np.repeat(data[i],2),label=method)
        plt.plot(np.arange(len(data_noisy[i])*2),np.repeat(data_noisy[i],2),'--',label=method+'_n')
plt.plot(np.arange(200 if len(nelder) > 200 else len(nelder)),nelder[:200],label='Nelder-Mead')
plt.plot(np.arange(200 if len(powell) > 200 else len(powell)),powell[:200],label='Powell')
plt.plot(np.arange(200 if len(nelder_n) > 200 else len(nelder_n)),nelder_n[:200],'--',label='Nelder-Mead_n')
plt.plot(np.arange(200 if len(powell_n) > 200 else len(powell_n)),powell_n[:200],'--',label='Powell_n')
plt.legend()
plt.xlabel('Iterations (2 function evaluations)')
plt.show()


