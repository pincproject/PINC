# coding: utf-8
import h5py
import numpy as np
fileRho = h5py.File('../../data/phi.grid.h5', 'r')
rho = fileRho['/n=1000.0']
rho_denorm = fileRho.attrs.__getitem__("Quantity denormalization factor")
rho = np.transpose(rho,(3,2,1,0))
rho = rho.squeeze()
rho = rho_denorm * rho
print(rho.max())
print(rho.min())
print(rho.mean())