
import h5py
#import numpy as np
#import pylab as plt


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

## Setup Params: #######

file_name = "phi"#"rhoNeutral" #"P"

ppc = 32 # particle per cell (for rho plots)

# timesteps:
timestep = 10#100000#29000#50600 #4950#50713#45715 # Must exist in dataset
#step = 1

# Plot:
levels = 500 ## granularity of contourf

#Restrict data values (can be values from 0-1):
restr_max = 1 # (0.5 = half of positive values)
restr_min = 1 #(1 = all of negative values)

cmap = 'jet'
line = 'X' # Y, Z


save_figs = True

########################



## Colormaps:
# : 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu','RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic', 'YlGnBu'
# : 'viridis', 'plasma', 'inferno', 'magma', 'cividis'




h5 = h5py.File('../../data/'+file_name+'.grid.h5','r')

dimen = h5.attrs["Axis denormalization factor"][0]
denorm = h5.attrs["Quantity denormalization factor"][0]

timesteps=[]
for item in h5:
	timesteps.append(item.split("=")[-1])
timesteps = np.array(timesteps,dtype =float)
timesteps = np.sort(timesteps)
#print(timesteps)

for i in range(len(timesteps)):
    if (timestep==timesteps[i]):
        start_index=i+1

timesteps = timesteps[start_index-1:] # should include an offset +- for averaging over n timesteps
DATA = []
print(timesteps)
for q in range(len(timesteps)):
	i = timesteps[q]
	dataset = h5["/n=%.1f"%i]
	#data = np.transpose(dataset,(3,2,1,0))	
	data = np.squeeze(dataset)
	#print(" ")
	#print(data.shape)
	data = np.transpose(data)
	if (line == 'X'):
		print("using pos %f in Y"%( int(len(data[0,0,:])/2)*dimen ))
		print("using pos %f in Z"%( int(len(data[0,0,:])/2)*dimen ))
		data = np.transpose((data[:,int(len(data[0,0,:])/2),int(len(data[0,0,:])/2) ]))*denorm #int(len(data[0,0,:])/2)

	if (line == 'Y'):
		print("using pos %f in X"%( int(len(data[0,0,:])/2)*dimen ))
		print("using pos %f in Z"%( int(len(data[0,0,:])/2)*dimen ))
		data = np.transpose((data[int(len(data[0,0,:])/2),:,int(len(data[0,0,:])/2)]))*denorm #int(len(data[0,0,:])/2)

	if (line == 'Z'):
		print("using pos %f in X"%( int(len(data[0,0,:])/2)*dimen ))
		print("using pos %f in Y"%( int(len(data[0,0,:])/2)*dimen ))
		data = np.transpose((data[int(len(data[0,0,:])/2),int(len(data[0,0,:])/2),:]))*denorm #int(len(data[0,0,:])/2)

	if("rho" in file_name ):
		if("rho_e" in file_name ):
			data = -1*data
		data = (data/denorm)/ppc # Number density
	#data = np.transpose(np.average(data,axis=2))*denorm 
	#print(data[0,0])
	DATA.append(data)


DATA = np.array(DATA)
print(DATA.shape)
avgData = DATA[0]
for i in range(1,len(DATA)):
	avgData += DATA[i]
avgData /= len(DATA)


print(avgData.shape)
plt.plot(range(len(avgData))*dimen, avgData)
plt.title("avg_last_%03d_dt_dim-"%(timesteps[-1]-timesteps[0])+line )
if (line == 'X'):
	plt.xlabel("X (m)")
	plt.ylabel(file_name)
if (line == 'Y'):
	plt.xlabel("Y (m)")
	plt.ylabel(file_name)
if (line == 'Z'):
	plt.xlabel("X (m)")
	plt.ylabel(file_name)
if(save_figs == True): 
            plt.savefig("anim_output/"+file_name+"_avg_last_%03d_dt_dim-"%(timesteps[-1]-timesteps[0])+line+".png")
plt.show()



