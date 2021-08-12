
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
start = 280#59500# # Must exist in dataset
#step = 1

# Plot:
levels = 500 ## granularity of contourf
interval = 0.1#in seconds

#Restrict data values (can be values from 0-1):
restr_max = 1 # (0.5 = half of positive values)
restr_min = 1 #(1 = all of negative values)

cmap = 'jet'
plane = 'XZ' # XY, XZ, YZ

show_plot = True 

save_figs = False#True


## Needs ffmpeg codec
save_anim = False ## Bool (if false anim is only shown on screen)
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
end = timesteps[-1]

for i in range(len(timesteps)):
	if timesteps[i] == start:
		timesteps = timesteps[i:]
		break


DATA = []
for i in timesteps:
	dataset = h5["/n=%.1f"%i]
	#data = np.transpose(dataset,(3,2,1,0))	
	data = np.squeeze(dataset)
	#print(" ")
	#print(data.shape)
	data = np.transpose(data)
	if (plane == 'XY'):
		data = np.transpose((data[:,:,int(len(data[0,0,:])/2) ]))*denorm #int(len(data[0,0,:])/2)
	if (plane == 'YZ'):
		data = np.transpose((data[int(len(data[0,0,:])/2),:,:]))*denorm #int(len(data[0,0,:])/2)
	if (plane == 'XZ'):
		data = np.transpose((data[:,int(len(data[0,0,:])/2),:]))*denorm #int(len(data[0,0,:])/2)
	if("rho" in file_name ):
		if("rho_e" in file_name ):
			data = -1*data
		if("rho" == file_name ):
			data = (data/denorm)/(ppc) # Number density
		else: data = (data/denorm)/(ppc)
	#data = np.transpose(np.average(data,axis=2))*denorm 
	#print(data[0,0])
	DATA.append(data)
DATA = np.array(DATA)
avgData = DATA[0]
for i in range(1,len(DATA)):
	avgData += DATA[i]
avgData /= len(DATA)


print("max value = %f"%np.amax(DATA))
print("min value = %f"%np.amin(DATA))


vMin=restr_min*np.amin(DATA)
vMax=restr_max*np.amax(DATA)

print("restricting values to %f, %f"%(vMin,vMax))

fig,ax = plt.subplots()
ax.set_aspect('equal')
from mpl_toolkits.axes_grid1 import make_axes_locatable

div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')


print("dx = %f"%dimen)
print("denorm = %e"%denorm )


x= np.linspace(0,dimen*(len(DATA[0,0,:])-1),len(DATA[0,0,:]))
y= np.linspace(0,dimen*(len(DATA[0,:,0])-1),len(DATA[0,:,0]))




X,Y = np.meshgrid(x,y) # not necessarily actual x, y dimensions

cax.cla()
#ax.clear()
img = ax.contourf(X, Y, avgData ,cmap = cmap , levels=levels, vmin=vMin,vmax=vMax)#,cmap = 'RdYlBu'

#img = ax.imshow(DATA[i,:,:],extent=[0,x[-1],0,y[-1]],cmap = cmap, vmin=vMin,vmax=vMax) 
ax.set_title(file_name+' Averaged over %03d timesteps'%(end-start) )
if (plane == 'XY'):
        ax.set_xlabel("X (m)")
        ax.set_ylabel("Y (m)")
if (plane == 'YZ'):
        ax.set_xlabel("Y (m)")
        ax.set_ylabel("Z (m)")
if (plane == 'XZ'):
        ax.set_xlabel("X (m)")
        ax.set_ylabel("Z (m)")
mesh = cax.pcolormesh(data, cmap = cmap)
mesh.set_clim(vMin,vMax)
if("rho" in file_name ):
    cax.set_title('n/n0' )
if("phi" in file_name ):
    cax.set_title('V' )
fig.colorbar(mesh,cax=cax)


if(save_figs == True): 
    plt.savefig("anim_output/"+file_name+"_Averaged_over_%03d_timesteps"%(end-start)+".png")

if (show_plot == True):
    plt.show()





### old ver


#for i in range(start,200000,1):
#	dataset = h5["/n=%.1f"%i]*denorm
#	data = np.squeeze(dataset)
#	#data = np.transpose(data,(2,1,0))
#	data = data[32,:,:]#np.average(data,axis=0)
#	if i==start:
#		p = plt.imshow(data, interpolation='bilinear', cmap = 'YlGnBu', vmin=-0.8,vmax=0.1) # 
#		
#		fig = plt.gcf()	# 

#				# interpolation= 'none', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos'.
#		plt.clim()
#		plt.title("Charge density, t=%i"%i);
#		plt.colorbar(orientation='horizontal')
#	else:
#		p.set_data(data)
#		plt.title("Charge density, t=%i"%i);

#	plt.pause(0.1)


