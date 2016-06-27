# To test that a (hopefully) sinusoidal small perturbation of a uniform
# positions makes indeed a sinusoidal charge density.

import h5py
import numpy as np
import pylab as plt




fileRho = h5py.File('../../test_rho.grid.h5','r')
rho = fileRho['/n=0.0']
rho = np.transpose(rho,(3,2,1,0))
rho = np.squeeze(rho)
print(rho.shape)

x = np.arange(rho.shape[0])
y = np.arange(rho.shape[1])

X,Y = np.meshgrid(x,y,indexing='ij')

fig, ax = plt.subplots(1)
im = ax.contourf(X,Y,rho[:,:,3], 50)

fig.subplots_adjust(bottom = 0.25)
cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")

plt.show()

