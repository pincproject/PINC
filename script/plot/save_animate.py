# Simple plot of a 2D slice through a single field quantity

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


h5 = h5py.File('../../data/phi.grid.h5','r')

#for i in range(1,100):
#	dataset = h5["/n=%.1f"%i]
#	data = np.squeeze(dataset)
#	#data = np.average(data,axis=1)		#animate average on axis
#	data = data[:][16][:]			#animate slice
#	if i==1:
#		p = plt.imshow(data)
#		fig = plt.gcf()
#		plt.clim()
#		plt.title("Charge density, t=%i"%i);
#		plt.colorbar(orientation='horizontal')
#	else:
#		p.set_data(data)
#		plt.title("Charge density, t=%i"%i);
#
#	plt.pause(0.0002)



def data_gen(t=0):
    for i in range(89001,90000):
	dataset = h5["/n=%.1f"%i]
	data = np.squeeze(dataset)
	#data = np.average(data,axis=1)		#animate average on axis
	data = data[:][16][:]			#animate slice
	yield data


dataset = h5["/n=89001.0"]
data = np.squeeze(dataset)
#data = np.average(data,axis=1)		#animate average on axis
data = data[:][16][:]
p = plt.imshow(data)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.clim()
plt.title("Charge density, t=89001.0");
plt.colorbar(orientation='horizontal')



def run(data):
    # update the data
    p.set_data(data)
    plt.title("Charge density");

ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=10,
                              repeat=False)
writer = animation.writers['ffmpeg'](fps=30)
ani.save('test.mp4',writer=writer,dpi=600)
plt.show()
