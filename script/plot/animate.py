# Simple plot of a 2D slice through a single field quantity

import h5py
import numpy as np
import pylab as plt


h5 = h5py.File('../../data/phi.grid.h5','r')

for i in range(89001,90000):
	dataset = h5["/n=%.1f"%i]
	data = np.squeeze(dataset)
	#data = np.average(data,axis=1)		#animate average on axis
	data = data[:][16][:]			#animate slice
	if i==89001:
		p = plt.imshow(data)
		fig = plt.gcf()
		plt.clim()
		plt.title("Charge density, t=%i"%i);
		plt.colorbar(orientation='horizontal')
	else:
		p.set_data(data)
		plt.title("Charge density, t=%i"%i);
	#plt.savefig("animation/dt%i"%i) #save to file
	plt.pause(0.0002)
