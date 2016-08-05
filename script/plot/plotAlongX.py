import h5py
import numpy as np
import pylab as plt

rhoFile = h5py.File('../../test_rho.grid.h5','r')
phiFile = h5py.File('../../test_phi.grid.h5','r')
EFile = h5py.File('../../test_E.grid.h5','r')

dataset = '/n=1.0';

rho = rhoFile[dataset]
phi = phiFile[dataset]
E = EFile[dataset]

# Change order of axes to get (d,x,y,z) like in PINC, d is dimension (0, 1 or 2)
rho = np.transpose(rho,(3,2,1,0))
phi = np.transpose(phi,(3,2,1,0))
E = np.transpose(E,(3,2,1,0))

# Squeeze scalar quantities, and extract x-component of field. Result is 3D arrays.
rho = np.squeeze(rho)
phi = np.squeeze(phi)
E = E[0,:,:,:]

# Average away y-axis, then z-axis. Result is a 1D (averaged) array.
rho = np.average(rho,axis=1)
rho = np.average(rho,axis=1)
phi = np.average(phi,axis=1)
phi = np.average(phi,axis=1)
E = np.average(E,axis=1);
E = np.average(E,axis=1);

plt.plot(rho,label='rho')
plt.plot(phi,label='phi')
plt.plot(E,label='E')
plt.legend(loc='upper right')
plt.show()
