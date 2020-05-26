import numpy as np
import matplotlib.pyplot as plt
import sys

avg_nums = [0,2,5,10]

method = sys.argv[1]

data = np.load(method+'.npy')
x = np.arange(len(data[0]))
for i,a in enumerate(avg_nums):
    plt.plot(x,data[i],label='avg_num={}'.format(a))
plt.legend()
plt.title('Average gradient of {}'.format(method))
plt.xlabel('Iterations (2 function evaluations)')
plt.ylabel('Energy [Hartree]')
plt.show()
