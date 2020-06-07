import numpy as np
import matplotlib.pyplot as plt
from tikzplotlib import save
from sys import argv

x = np.load('../bonds.npy')
fci = np.load('../../H2/cc.npy')
fci = fci[[i for i in range(0,len(fci),5)]]
x_full = np.load('../../H2/bonds.npy')
fci_full = np.load('../../H2/cc.npy')
hf = np.load('../../H2/hf.npy')

ry = np.load('single_param/E_RY.npy')
uccd = np.load('single_param/E_UCCD.npy')
uccdr = np.load('single_param/E_UCCDr.npy')

ryrz = np.load('five_param/E_RYRZ_s100.npy')
uccsd = np.load('five_param/E_UCCSD.npy')
uccsdr = np.load('five_param/E_UCCSDr_s100.npy')



plt.figure()
plt.title('Noisy simulation single param')
plt.plot(x_full,fci_full,label='FCI')
plt.plot(x_full,hf,label='HF')
plt.plot(x,ry,'.',label='RY')
plt.plot(x,uccd,'+',label='UCCD')
plt.plot(x,uccdr,'*',label='UCCDr')
plt.legend()

save('tex/h2_n_single_param.tex')

plt.figure()
plt.title('Noisy simulation single param')
plt.plot(x,np.abs(fci-ry),'.',label='RY')
plt.plot(x,np.abs(fci-uccd),'+',label='UCCD')
plt.plot(x,np.abs(fci-uccdr),'*',label='UCCDr')
plt.legend()

save('tex/h2_n_abs_single_param.tex')

plt.figure()
plt.title('Noisy simulation multi param')
plt.plot(x_full,fci_full,label='FCI')
plt.plot(x_full,hf,label='HF')
plt.plot(x,ryrz,'.',label='RYRZ')
plt.plot(x,uccsd,'+',label='UCCSD')
plt.plot(x,uccsdr,'*',label='UCCSDr')
plt.legend()

save('tex/h2_n_five_param.tex')

plt.figure()
plt.title('Noisy simulation multi param')
plt.plot(x,np.abs(fci-ryrz),'.',label='RYRZ')
plt.plot(x,np.abs(fci-uccsd),'+',label='UCCSD')
plt.plot(x,np.abs(fci-uccsdr),'*',label='UCCSDr')
plt.legend()

save('tex/h2_n_abs_five_param.tex')

plt.show()
