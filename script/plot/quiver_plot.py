
import matplotlib.pyplot as plt
import numpy as np
import h5py

file_name = "E"#"gradBulkV"

xSize = 64
ySize = 64

plane = 'YZ' 

timestep = 100 #200

h5 = h5py.File('../../data/'+file_name+'.grid.h5','r')

dimen = h5.attrs["Axis denormalization factor"][0]
denorm = h5.attrs["Quantity denormalization factor"][0]
print(denorm)

#DATA = []

dataset = h5["/n=%.1f"%timestep]
#data = np.transpose(dataset,(3,2,1,0))	
data = np.squeeze(dataset)
#print(" ")

vMin=0.1*np.amin(data)
vMax=0.1*np.amax(data)

data = np.squeeze(data)/vMax

X, Y = np.meshgrid(np.arange(0, xSize, 1), np.arange(0, ySize, 1))

x_shape = X.shape

U = np.zeros(x_shape)
V = np.zeros(x_shape)

#data = np.transpose(data)
print(data[:,:,:,:].shape)
if (plane == 'XY'):
	U = data[int(len (data[0,0,:]) /2),:,:,2] #[Z,Y,X,(z,y,x)]
	V = data[int(len (data[0,0,:]) /2),:,:,1]
	#data = np.transpose((data[:,:,:,int(len (data[0,0,:]) /2)] )) #*denorm #int(len(data[0,0,:])/2)

if (plane == 'YZ'):
	U = data[:,:,int(len (data[0,0,:]) /2),1] #[Z,Y,X,(z,y,x)]
	V = data[:,:,int(len (data[0,0,:]) /2),0]
	#data = np.transpose((  data[:,int(len(data[0,0,:])/2),:,:]   )) #*denorm #int(len(data[0,0,:])/2)
if (plane == 'XZ'):
	U = data[:,int(len (data[0,0,:]) /2),:,0] #[Z,Y,X,(z,y,x)]
	V = data[:,int(len (data[0,0,:]) /2),:,2]
	#data = np.transpose((data[:,:,int(len(data[0,0,:])/2),:])) #*denorm #int(len(data[0,0,:])/2)

#data = np.transpose(np.average(data,axis=2))*denorm 
#print(data[0,0])
#DATA= data
##DATA = np.array(DATA)



print(vMin)
print(vMax)





print(data.shape)

#if (plane == 'XY'):
#	U = DATA[:,:,1]
#	V = DATA[:,:,2]
#if (plane == 'YZ'):
#	U = DATA[:,:,1]
#	V = DATA[:,:,0]
#if (plane == 'XZ'):
#	U = DATA[:,:,0]
#	V = DATA[:,:,2]

#print(U.shape)

fig, ax = plt.subplots()
q = ax.quiver(X, Y, U, V, units='xy' ,scale=1, color='red')

ax.set_aspect('equal')

plt.xlim(0,xSize)
plt.ylim(0,ySize)

plt.title(file_name,fontsize=10)

#plt.savefig('how_to_plot_a_vector_field_in_matplotlib_fig1.png', bbox_inches='tight')
plt.show()
#plt.close()


