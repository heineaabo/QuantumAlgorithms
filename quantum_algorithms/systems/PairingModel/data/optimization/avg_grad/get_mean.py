import numpy as np
import sys

avg_nums = [0,2,5,10]

method = sys.argv[1]
val = int(sys.argv[2])

data = np.load(method+'.npy')
print('Mean and variance for',method+':')
print('avg_num   Mean        Variance')
for i,a in enumerate(avg_nums):
    mean = np.mean(data[i][-val:])
    var = np.var(data[i][-val:])
    print('{:<10}{:<12.8f}{:.8f}'.format(a,mean,var))
