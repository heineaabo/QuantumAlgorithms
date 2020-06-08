import numpy as np
import matplotlib.pyplot as plt

x   = np.load('omegas.npy')
x2  = [0.5,1.0,1.5,2.0]
fci = np.load('fci.npy')
hf  = np.load('hf.npy')
cc  = np.load('cc.npy')

#ucc = np.load('UCCSDr_E_i.npy')
ucc_new = np.load('UCCSDr_E_new.npy')
print(len(ucc_new))

plt.plot(x,fci,label='FCI')
plt.plot(x,hf,label='HF')
#plt.plot(x,ucc,'.',label='UCCSDr')
plt.plot(x2,ucc_new,'*',label='UCCSDr_new')
plt.legend()

plt.show()

