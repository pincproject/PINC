#Some rudimentary plotting of distributions

import h5py
import numpy as np
import pylab as plt

fileRho = h5py.File('../test_rho.grid.h5','r')
rho = fileRho['/n=0.0']

filePhi = h5py.File('../test_phi.grid.h5','r')
phi = filePhi["/n=0.0"]


x = np.arange(rho.shape[1])
y = np.arange(rho.shape[0])

X,Y = np.meshgrid(x,y)


#Plotting
plt.figure()
contourRho = plt.contourf(X,Y,rho)
plt.colorbar(contourRho)

plt.figure()
contourPhi = plt.contourf(X,Y,phi)
plt.colorbar(contourPhi)

plt.show()
