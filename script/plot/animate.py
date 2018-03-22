# Simple plot of a 2D slice through a single field quantity

import h5py
import numpy as np
import pylab as plt


h5 = h5py.File('../../data/rho.grid.h5','r')

#Nt = h5["/n=90000.0"].shape[1]	# Number of timesteps
#print(Nt)
startindex =10000
stopindex = 20000

"""
for i in range(startindex,stopindex,1000):#start and stop timestep
	dataset = h5["/n=%.1f"%i]
	data = np.squeeze(dataset)
	data = np.average(data,axis=2)		#animate average on axis
	#data = data[:][:][8]			#animate slice
	if i==startindex:
		p = plt.imshow(data)
		fig = plt.gcf()
		#plt.clim(vmin=-0.0016, vmax=0.0016)
		plt.clim()    # auto range color
		plt.title("Potential, t=%i"%i);
		plt.colorbar(orientation='horizontal')
		plt.savefig("animation/dt%i"%i) #save to file
	else:
		p.set_data(data)
		#plt.clim()
		plt.title("Potential, t=%i"%i)
		plt.savefig("animation/dt%i"%i) #save to file
	
	plt.pause(0.1)

#alternate saves .png to make animations with e.g ffmpg
"""
count = 0
for i in range(startindex,stopindex,1000):#start and stop timestep
	dataset = h5["/n=%.1f"%i]
	data = np.transpose(dataset,(3,2,1,0))
	data = np.squeeze(data)
	#print(data.shape)

	x = np.arange(data.shape[0])
	y = np.arange(data.shape[1])

	X,Y = np.meshgrid(x,y,indexing='ij')

	fig, ax = plt.subplots(1)
	im = ax.contourf(X,Y,data[:,:,32], 100)

	fig.subplots_adjust(bottom = 0.25)
	cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
	fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")
	plt.title("Ion density (rho_e) perpendicular to B_0, t=%i"%i);
	
	plt.savefig("animation/dt%i"%count) #save to file
	count +=1
	#plt.show()




