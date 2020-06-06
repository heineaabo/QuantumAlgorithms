import numpy as np
import matplotlib.pyplot as plt
from tikzplotlib import save

fci = np.load('../cc.npy')
hf = np.load('../hf.npy')
x1 = np.load('../bonds.npy')
x2 = x1[[i for i in range(0,len(x1),5)]]
fci2 = fci[[i for i in range(0,len(x1),5)]]

ryrz = np.load('RYRZ_good_E.npy')
ryrz_n = np.load('RYRZ_noisy_E.npy')
uccsd = np.load('UCCSD_good_E.npy')
uccsd_n = np.load('UCCSD_noisy_E.npy')
uccsdr = np.load('E_ideal_UCCSDr.npy')

plt.figure()
plt.plot(x1,fci,label='FCI')
plt.plot(x1,hf,label='HF')
plt.plot(x2,ryrz,'.',label='RYRZ')
plt.plot(x2,uccsd,'+',label='UCCSD')
plt.plot(x2,uccsdr,'+',label='UCCSDr')
plt.title('Ideal simulation')
plt.legend()

save('../tex/h2_bond_ideal_5.tex')

plt.figure()
plt.plot(x1,fci,label='FCI')
plt.plot(x1,hf,label='HF')
plt.plot(x2,ryrz_n,'.',label='RYRZ')
plt.plot(x2,uccsd_n,'+',label='UCCSD')
plt.title('Noisy simulation')
plt.legend()

save('../tex/h2_bond_noisy_5.tex')

plt.figure()
plt.fill_between(x2,0,0.0016, alpha=0.2,label='Chemical accuracy')
plt.plot(x2,np.abs(fci2-ryrz),'.',label='RYRZ')
plt.plot(x2,np.abs(fci2-uccsd),'+',label='UCCSD')
plt.plot(x2,np.abs(fci2-uccsdr),'*',label='UCCSDr')
plt.legend()

save('../tex/h2_error_5.tex')

plt.show()
