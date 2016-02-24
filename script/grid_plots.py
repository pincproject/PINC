#Some rudimentary plotting of distributions

import h5py
import numpy as np
import pylab as plt

fileRho = h5py.File('../test_rho.grid.h5','r')
fileRes = h5py.File('../test_res.grid.h5','r')
filePhi = h5py.File('../test_phi.grid.h5','r')

#Time steps
rho = fileRho['/n=0.0']
res = fileRes['/n=0.0']
phi0 = filePhi["/n=0.0"]
# phi1 = filePhi["/n=1.0"]
# phi2 = filePhi["/n=2.0"]
# phi3 = filePhi["/n=3.0"]

# #Extracting wanted layers 3D
# rho = rho[0,:,:,0]
#
# phi0 = phi0[0,:,:,0]
# phi1 = phi1[0,:,:,0]
# phi2 = phi2[0,:,:,0]
# phi3 = phi3[0,:,:,0]

#2D
rho = rho[:,:,0]
res = res[:,:,0]

# phi0 = phi0[:,:,0]
# phi1 = phi1[:,:,0]
# phi2 = phi2[:,:,0]
# phi3 = phi3[:,:,0]



x = np.arange(rho.shape[1])
y = np.arange(rho.shape[0])

X,Y = np.meshgrid(x,y)
#
#
# vmin = np.min(phi3)-1
# vmax = np.max(phi3)+1
# # # vmin = -100
# # # vmax = 200
# # # vmin = -1
# # # vmax = 4
# #
# levels=np.arange(vmin,vmax, np.abs(vmax)/50.)
# # # levels=np.arange(0, vmax)

fig, ax = plt.subplots(1)
im = ax.contourf(X,Y,rho, 50)

fig.subplots_adjust(bottom = 0.25)
cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")

fig, ax = plt.subplots(1)
im = ax.contourf(X,Y,res, 50)

fig.subplots_adjust(bottom = 0.25)
cbar_res = fig.add_axes([0.10, 0.05, 0.8, 0.10])
fig.colorbar(im, cax=cbar_res, orientation = "horizontal")



# fig, ax = plt.subplots(2,2, sharex= True , sharey = True)
# # im = ax[0,0].contourf(X,Y,phi0, levels = levels)
# # ax[0,1].contourf(X,Y,phi1, levels = levels)
# # ax[1,0].contourf(X,Y,phi2, levels = levels)
# # ax[1,1].contourf(X,Y,phi3, levels = levels)
#
# im = ax[0,0].contourf(X,Y,phi0, 50)
# ax[0,1].contourf(X,Y,phi1, 50)
# ax[1,0].contourf(X,Y,phi2, 50)
# ax[1,1].contourf(X,Y,phi3, 50)
#
# fig.subplots_adjust(bottom = 0.25)
#
# cbar_ax = fig.add_axes([0.10, 0.05, 0.8, 0.10])
# fig.colorbar(im, cax=cbar_ax, orientation = "horizontal")
# #
#
plt.show()
