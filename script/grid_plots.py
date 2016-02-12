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
phi6 = filePhi["/n=6.0"]
phi6 = phi6[:,:,0]
phi7 = filePhi["/n=7.0"]
phi7 = phi7[:,:,0]
phi8 = filePhi["/n=8.0"]
phi8 = phi8[:,:,0]
phi9 = filePhi["/n=9.0"]
phi9 = phi9[:,:,0]
phi10 = filePhi["/n=10.0"]
phi10 = phi10[:,:,0]
phi11 = filePhi["/n=11.0"]
phi11 = phi11[:,:,0]


x = np.arange(rho.shape[1])
y = np.arange(rho.shape[0])

X,Y = np.meshgrid(x,y)


vmin = np.min(phi11)-1
vmax = np.max(phi11)+1
# vmin = -100
# vmax = 200
# vmin = -1
# vmax = 4

levels=np.arange(vmin,vmax, np.abs(vmax - vmin)/50.)

fig, ax = plt.subplots(1)
im = ax.contourf(X,Y,rho)

fig.subplots_adjust(bottom = 0.25)
cbar_ax = fig.add_axes([0.10, 0.05, 0.8, 0.10])
fig.colorbar(im, cax=cbar_ax, orientation = "horizontal")


fig, ax = plt.subplots(4,3, sharex= True , sharey = True)
# im = ax[0,0].contourf(X,Y,phi0, levels = levels)
# ax[0,1].contourf(X,Y,phi1, levels = levels)
# ax[0,2].contourf(X,Y,phi2, levels = levels)
# ax[1,0].contourf(X,Y,phi3, levels = levels)
# ax[1,1].contourf(X,Y,phi4, levels = levels)
# ax[1,2].contourf(X,Y,phi5, levels = levels)
# ax[2,0].contourf(X,Y,phi6, levels = levels)
# ax[2,1].contourf(X,Y,phi7, levels = levels)
# ax[2,2].contourf(X,Y,phi8, levels = levels)
# ax[3,0].contourf(X,Y,phi9, levels = levels)
# ax[3,1].contourf(X,Y,phi10, levels = levels)
# ax[3,2].contourf(X,Y,phi11, levels = levels)

im = ax[0,0].contourf(X,Y,phi0, 50)
ax[0,1].contourf(X,Y,phi1, 50)
ax[0,2].contourf(X,Y,phi2, 50)
ax[1,0].contourf(X,Y,phi3, 50)
ax[1,1].contourf(X,Y,phi4, 50)
ax[1,2].contourf(X,Y,phi5, 50)
ax[2,0].contourf(X,Y,phi6, 50)
ax[2,1].contourf(X,Y,phi7, 50)
ax[2,2].contourf(X,Y,phi8, 50)
ax[3,0].contourf(X,Y,phi9, 50)
ax[3,1].contourf(X,Y,phi10, 50)
ax[3,2].contourf(X,Y,phi11, 50)


fig.subplots_adjust(bottom = 0.25)

cbar_ax = fig.add_axes([0.10, 0.05, 0.8, 0.10])
fig.colorbar(im, cax=cbar_ax, orientation = "horizontal")


plt.show()
