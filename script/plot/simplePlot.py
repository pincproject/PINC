# Simple plot of a 2D slice through a single field quantity

import h5py
import numpy as np
import pylab as plt


fileRho = h5py.File('../../data/phi.grid.h5','r')
rho = fileRho['/n=500.0']
rho = np.squeeze(rho)
rho = np.transpose(rho)
rho_denorm = fileRho.attrs.__getitem__("Quantity denormalization factor")
print(rho_denorm)

rho = rho_denorm * rho
print(rho.shape)
x = np.arange(rho.shape[0])
y = np.arange(rho.shape[1])

X,Y = np.meshgrid(x,y,indexing='ij')

fig, ax = plt.subplots(1)
im = ax.contourf(X,Y,rho[:,:,32], 50)

fig.subplots_adjust(bottom = 0.25)
cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")
plt.show()
