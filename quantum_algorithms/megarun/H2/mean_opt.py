import numpy as np
import matplotlib.pyplot as plt
from tikzplotlib import save
from pandas import DataFrame

print('FCI   :  -1.1372838344855132')


# RY ansatz
cob = np.load('cobyla/E_R074_final_RYPAIRING.npy')
powl= np.load('powell/E_R074_final_RYPAIRING.npy')
nel = np.load('nelder-mead/E_R074_final_RYPAIRING.npy')
spsa= np.load('SPSA/E_R074_final_RYPAIRING.npy')
cobA = np.load('cobyla/E_conv_R074_RYPAIRING.npy')
powlA= np.load('powell/E_conv_R074_RYPAIRING.npy')
nelA = np.load('nelder-mead/E_conv_R074_RYPAIRING.npy')
spsaA= np.load('SPSA/E_conv_R074_RYPAIRING.npy')
print('RY ansatz')
print('Cobyla: ',cob,len(cobA))
print('Powell: ',powl,len(powlA))
print('Nelder: ',nel,len(nelA))
print('SPSA  : ',spsa,len(spsaA))
print('')
RY_E = [cob,powl,nel,spsa]
RY_N = [len(cobA),len(powlA),len(nelA),len(spsaA)]

# UCCD ansatz
cob = np.load('cobyla/E_R074_final_UCCD.npy')
powl= np.load('powell/E_R074_final_UCCD.npy')
nel = np.load('nelder-mead/E_R074_final_UCCD.npy')
spsa= np.load('SPSA/E_R074_final_UCCD.npy')
cobA = np.load('cobyla/E_conv_R074_UCCD.npy')
powlA= np.load('powell/E_conv_R074_UCCD.npy')
nelA = np.load('nelder-mead/E_conv_R074_UCCD.npy')
spsaA= np.load('SPSA/E_conv_R074_UCCD.npy')
print('UCCD ansatz')
print('Cobyla: ',cob,len(cobA))
print('Powell: ',powl,len(powlA))
print('Nelder: ',nel,len(nelA))
print('SPSA  : ',spsa,len(spsaA))
print('')

# UCCDr ansatz
cob = np.load('cobyla/E_R074_final_UCCDr.npy')
powl= np.load('powell/E_R074_final_UCCDr.npy')
nel = np.load('nelder-mead/E_R074_final_UCCDr.npy')
spsa= np.load('SPSA/E_R074_final_UCCDr.npy')
cobA = np.load('cobyla/E_conv_R074_UCCDr.npy')
powlA= np.load('powell/E_conv_R074_UCCDr.npy')
nelA = np.load('nelder-mead/E_conv_R074_UCCDr.npy')
spsaA= np.load('SPSA/E_conv_R074_UCCDr.npy')
print('UCCDr ansatz')
print('Cobyla: ',cob,len(cobA))
print('Powell: ',powl,len(powlA))
print('Nelder: ',nel,len(nelA))
print('SPSA  : ',spsa,len(spsaA))
print('')
UCCD_E = [cob,powl,nel,spsa]
UCCD_N = [len(cobA),len(powlA),len(nelA),len(spsaA)]

# UCCSD ansatz
cob = np.load('cobyla/E_R074_final_UCCSD.npy')
powl= np.load('powell/E_R074_final_UCCSD.npy')
nel = np.load('nelder-mead/E_R074_final_UCCSD.npy')
spsa= np.load('SPSA/E_R074_final_UCCSD.npy')
cobA = np.load('cobyla/E_conv_R074_UCCSD.npy')
powlA= np.load('powell/E_conv_R074_UCCSD.npy')
nelA = np.load('nelder-mead/E_conv_R074_UCCSD.npy')
spsaA= np.load('SPSA/E_conv_R074_UCCSD.npy')
print('UCCSD ansatz')
print('Cobyla: ',cob,len(cobA))
print('Powell: ',powl,len(powlA))
print('Nelder: ',nel,len(nelA))
print('SPSA  : ',spsa,len(spsaA))
print('')

# UCCSDr ansatz
cob = np.load('cobyla/E_R074_final_UCCSDr.npy')
powl= np.load('powell/E_R074_final_UCCSDr.npy')
nel = np.load('nelder-mead/E_R074_final_UCCSDr.npy')
spsa= np.load('SPSA/E_R074_final_UCCSDr.npy')
cobA = np.load('cobyla/E_conv_R074_UCCSDr.npy')
powlA= np.load('powell/E_conv_R074_UCCSDr.npy')
nelA = np.load('nelder-mead/E_conv_R074_UCCSDr.npy')
spsaA= np.load('SPSA/E_conv_R074_UCCSDr.npy')
print('UCCSDr ansatz')
print('Cobyla: ',cob,len(cobA))
print('Powell: ',powl,len(powlA))
print('Nelder: ',nel,len(nelA))
print('SPSA  : ',spsa,len(spsaA))
UCCSD_E = [cob,powl,nel,spsa]
UCCSD_N = [len(cobA),len(powlA),len(nelA),len(spsaA)]

rnd = 5
index = ['Cobyla','Powell','Nelder','SPSA']
columns = ['RY_E','RY_N','UCCD_E','UCCD_N','UCCSD_E','UCCSD_N']
print()
data = {'RY_E':np.round(np.asarray(RY_E),rnd), 'RY_N':RY_N,
        'UCCD_E':np.round(np.asarray(UCCD_E),rnd), 'UCCD_N':UCCD_N,
        'UCCSD_E':np.round(np.asarray(UCCSD_E),rnd), 'UCCSD_N':UCCSD_N}

df = DataFrame(data,columns=columns,index=index)
print(df.to_latex())
