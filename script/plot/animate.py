# Simple plot of a 2D slice through a single field quantity

import h5py
import numpy as np
import pylab as plt



h5 = h5py.File('../../data/phi.grid.h5','r')

start = 1
for i in range(start,200000,1):
	dataset = h5["/n=%.1f"%i]
	data = np.squeeze(dataset)
	data = data[3,:,:]#np.average(data,axis=0)
	if i==start:
		p = plt.imshow(data)
		fig = plt.gcf()
		plt.clim()
		plt.title("Charge density, t=%i"%i);
		plt.colorbar(orientation='horizontal')
	else:
		p.set_data(data)
		plt.title("Charge density, t=%i"%i);

	plt.pause(0.1)
