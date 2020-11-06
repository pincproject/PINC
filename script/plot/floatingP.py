import numpy as np
from scipy.special import erfc
import matplotlib.pyplot as plt
import scipy.constants as c

ne = 1e8
k = c.Boltzmann
me = c.electron_mass
q = c.elementary_charge
Te = 1.16e6
r = 0.9
l = 0.9


V = np.linspace(0,200,4000)
y = np.zeros_like(V)

y.fill(1.55e-13)

#cyl = (2/(np.sqrt(np.pi))) * np.sqrt(V) + np.exp(V)*erfc(np.sqrt(V))

cyl_alt = ne*q*np.sqrt((k*Te)/(2*np.pi*me))*2*np.pi*r*l * 2/np.sqrt(np.pi) * np.sqrt(1 + ((q*V)/(k*Te)))

plt.plot(V,cyl_alt)
plt.plot(V,y)
#plt.plot(x,cyl)
plt.show()