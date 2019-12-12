# Simple plot of a 2D slice through a single field quantity

import h5py
import numpy as np
import pylab as plt
plt.matplotlib.use("TkAgg")


h5 = h5py.File('../../data/rho.grid.h5','r')

start = 50501
for i in range(start,200000,1):
	dataset = h5["/n=%.1f"%i]
	data = np.squeeze(dataset)
	#data = np.transpose(data,(2,1,0))
	data = data[:,:,32]#np.average(data,axis=0)
	if i==start:
		p = plt.imshow(data, vmin=-20, vmax=20)
		fig = plt.gcf()
		plt.clim()
		plt.title("Charge density, t=%i"%i);
		plt.colorbar(orientation='horizontal')
	else:
		p.set_data(data)
		plt.title("Charge density, t=%i"%i);

	plt.pause(0.1)
