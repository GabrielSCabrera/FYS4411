import matplotlib.pyplot as plt
import numpy as np

file = np.sort(np.loadtxt('E_5.dat'))[::-1]
plt.plot(np.arange(1,len(file)+1), file)
plt.show()
plt.semilogy(np.arange(1,len(file)+1), file)
plt.show()
