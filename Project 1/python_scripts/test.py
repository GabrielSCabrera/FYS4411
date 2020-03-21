import numpy as np
import matplotlib.pyplot as plt

h_bar = 6.582e-16  # eVs
# h_bar = 10
mu = 939.6e6/(2*299_792_458**2) # eV / (m/s)^2
# mu = 939.6e6/(2*299_792_458**2) # eV / (m/s)^2
R = 2.1e-15 # m
V0 = 35e6 # eV

E = np.linspace(-35e6, 0, 100000)

def left(E):
    k = np.sqrt(2*mu*(E+V0)) / h_bar
    return k*np.cos(k*R)/np.sin(k*R)

def right(E):
    kappa = np.sqrt(-2*mu*E) / h_bar
    return -kappa


plt.plot(E, left(E))
plt.plot(E, right(E))
plt.ylim(-1E-15, 0.5E-15)
plt.xlabel("Energi (eV)")
plt.show()
