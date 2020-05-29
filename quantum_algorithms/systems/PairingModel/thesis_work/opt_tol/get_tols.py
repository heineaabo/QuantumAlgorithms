import sys
import numpy as np
from pandas import DataFrame

data = np.load(sys.argv[1])
device = sys.argv[1].split('_')
device = device[1]

col = ['normal','noisy_n','noisy_w']
index = ['1e03','1e04']

n,m,z = data.shape

df_data = {}

for j in range(m):
    varz = []
    for i in range(n):
        v = np.ceil(np.log10(np.var(data[i,j,:])))
        varz.append(np.power(10,v))
    df_data[col[j]] = varz

dataframe = DataFrame(df_data,columns=col,index=index)

print('On Quantum computer {}:'.format(device))
print('')
print(dataframe)

