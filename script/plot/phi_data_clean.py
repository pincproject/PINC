# coding: utf-8
import h5py
import numpy as np
import matplotlib.pyplot as plt
fileRho = h5py.File('../../data/phi.grid.h5', 'r')
rhos = []


for i in range(1,10,1):
    rho = fileRho['/n=' + str(i) + '.0']
    rho_denorm = fileRho.attrs.__getitem__("Quantity denormalization factor")
    rho = np.transpose(rho,(3,2,1,0))
    rho = rho.squeeze()
    rho = rho_denorm * rho
    #print(rho.shape)
    rhos.append(rho.mean())

x = np.arange(1,10,1)
y = np.asarray(rhos)


fig, ax = plt.subplots()
ax.plot(x, y)
ax.set(xlabel='timestep', ylabel='Max Voltage (V)', title='Object voltage trend over timesteps')
ax.grid()
#fig.savefig("trend.png")
plt.show()
#print(rhos.max)