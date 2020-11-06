
import h5py
#import numpy as np
#import pylab as plt


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

## Setup Params: #######

file_type = "phi"
file_names = ["PINC-daed1-polar-mag/data/","PINC-daed2-polar-mag/data/","PINC-daed3-polar-mag/data/","PINC-daed4-polar-mag/data/"]#"rhoNeutral" #"P"

ppc = 32 # particle per cell (for rho plots)

# timesteps:
use_last_timesteps = 500#50600 #4950#50713#45715 # Must exist in dataset
#step = 1

# Plot:
levels = 500 ## granularity of contourf

#Restrict data values (can be values from 0-1):
restr_max = 1 # (0.5 = half of positive values)
restr_min = 1 #(1 = all of negative values)

cmap = 'jet'
line = 'Z' # X, Y, Z
offset1 = -20
offset2 = 0


save_figs = False#True

########################



## Colormaps:
# : 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu','RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic', 'YlGnBu'
# : 'viridis', 'plasma', 'inferno', 'magma', 'cividis'



def plot_line_in_file(file_name):
	h5 = h5py.File(file_name+file_type+'.grid.h5','r')

	label_name = file_name.split("/")[0].split("-")[1]
	dimen = h5.attrs["Axis denormalization factor"][0]
	denorm = h5.attrs["Quantity denormalization factor"][0]

	timesteps=[]
	for item in h5:
		timesteps.append(item.split("=")[-1])
	timesteps = np.array(timesteps,dtype =float)
	timesteps = np.sort(timesteps)


	#for i in range(len(timesteps)):
	    #if (timestep==timesteps[i]):
	        #start_index=i
	timesteps=timesteps[-use_last_timesteps:]

	#timesteps = timesteps[start_index-1:] # should include an offset +- for averaging over n timesteps
	DATA = []

	for q in range(len(timesteps)):
		i = timesteps[q]
		dataset = h5["/n=%.1f"%i]
		#data = np.transpose(dataset,(3,2,1,0))
		data = np.squeeze(dataset)

		data = np.transpose(data)
		if (line == 'X'):
			data = np.transpose((data[:,int(len(data[0,:,0])/2)+offset1,int(len(data[0,0,:])/2)+offset2 ]))*denorm #int(len(data[0,0,:])/2)
		if (line == 'Y'):
			data = np.transpose((data[int(len(data[:,0,0])/2)+offset1,:,int(len(data[0,0,:])/2)+offset2]))*denorm #int(len(data[0,0,:])/2)
		if (line == 'Z'):
			data = np.transpose((data[int(len(data[:,0,0])/2)+offset1,int(len(data[0,:,0])/2)+offset2,:]))*denorm #int(len(data[0,0,:])/2)
		if("rho" in file_type ):
			if("rho_e" in file_type  ):
				data = -1*data
			data = (data/denorm)/ppc # Number density
		#data = np.transpose(np.average(data,axis=2))*denorm

		DATA.append(data)
	DATA = np.array(DATA)
	data = DATA[0]
	for i in range(1,len(DATA)):
		data+=DATA[i]
	data /= use_last_timesteps



	x = np.linspace(0,len(data)-1,len(data))*dimen

	plt.plot(x,data,label=label_name)
	return timesteps,denorm,dimen


for file_name in file_names:
	timesteps,denorm,dimen = plot_line_in_file(file_name)


if (line == 'X'):
	plt.xlabel("X (m)")
	plt.ylabel(file_type )
if (line == 'Y'):
	plt.xlabel("Y (m)")
	plt.ylabel(file_type )
if (line == 'Z'):
	plt.xlabel("Z (m)")
	plt.ylabel(file_type )


plt.legend()
plt.title("offset (1,2) = center+(%.1f,%.1f) m, averaged over %03d"%((offset1*dimen),(offset2*dimen),use_last_timesteps))
if(save_figs == True):
            plt.savefig(file_type +"_timestep_%03d"%(timesteps[1])+".png")
plt.show()
