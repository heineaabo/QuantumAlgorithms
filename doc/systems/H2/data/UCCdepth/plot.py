import numpy as np
import matplotlib.pyplot as plt

x = np.load('depths.npy')
legal = np.load('legal.npy')
illegal = np.load('illegal.npy')
legal_unopt = np.load('legal_unopt.npy')
illegal_unopt = np.load('illegal_unopt.npy')

total = legal + illegal
legal /= total
illegal /= total

total_unopt = legal_unopt + illegal_unopt
legal_unopt /= total_unopt
illegal_unopt /= total_unopt

plt.plot(x,legal,label='Legal states (opt)')
plt.plot(x,illegal,label='Illegal states (opt)')
plt.plot(x,legal_unopt,label='Legal states (not opt)')
plt.plot(x,illegal_unopt,label='Illegal states (not opt)')
plt.legend()
plt.show()
