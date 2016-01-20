#Some rudimentary plotting of distributions

import h5py
import numpy as np
import pylab as plt

fileRho = h5py.File('../test_rho.grid.h5','r')
rho = fileRho['/n=0.0']
rho = rho[:,:,0]

filePhi = h5py.File('../test_phi.grid.h5','r')

phi0 = filePhi["/n=0.0"]
phi0 = phi0[:,:,0]
phi1 = filePhi["/n=1.0"]
phi1 = phi1[:,:,0]
phi2 = filePhi["/n=2.0"]
phi2 = phi2[:,:,0]
phi3 = filePhi["/n=3.0"]
phi3 = phi3[:,:,0]
phi4 = filePhi["/n=4.0"]
phi4 = phi4[:,:,0]
phi5 = filePhi["/n=5.0"]
phi5 = phi5[:,:,0]

x = np.arange(rho.shape[1])
y = np.arange(rho.shape[0])

X,Y = np.meshgrid(x,y)


vmin = np.min(phi5)
vmax = np.max(phi5)
levels=np.arange(vmin,vmax, 2./50.)

fig, ax = plt.subplots(1)
im = ax.contourf(X,Y,rho)

fig.subplots_adjust(bottom = 0.25)
cbar_ax = fig.add_axes([0.10, 0.05, 0.8, 0.10])
fig.colorbar(im, cax=cbar_ax, orientation = "horizontal")


fig, ax = plt.subplots(2,3, sharex= True , sharey = True)
ax[0,0].contourf(X,Y,phi0, levels = levels)
ax[0,1].contourf(X,Y,phi1, levels = levels)
ax[0,2].contourf(X,Y,phi2, levels = levels)
ax[1,0].contourf(X,Y,phi3, levels = levels)
ax[1,1].contourf(X,Y,phi4, levels = levels)
ax[1,2].contourf(X,Y,phi5, levels = levels)

fig.subplots_adjust(bottom = 0.25)

cbar_ax = fig.add_axes([0.10, 0.05, 0.8, 0.10])
fig.colorbar(im, cax=cbar_ax, orientation = "horizontal")

print np.max(phi4)
print np.max(phi5)

plt.show()
