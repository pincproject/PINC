#Some rudimentary plotting of distributions

import h5py
import numpy as np
import pylab as plt
# from mayavi import mlab


# def plotPlanesOfGrid(name, grid):
# 	mlab.figure()
# 	im = mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(grid),
# 	                            plane_orientation='x_axes',
# 	                            slice_index=grid.shape[0]/2,
# 	                        )
#
# 	mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(grid),
# 	                            plane_orientation='y_axes',
# 	                            slice_index=grid.shape[1]/2,
# 	                        )
# 	mlab.title(name)
# 	mlab.colorbar(im)
# 	mlab.axes()
#
# 	return

def transformData(dataset):
	grid = dataset['/n=0.0']
	grid = np.transpose(grid, (3,2,1,0))
	grid = np.squeeze(grid)

	return grid


def plot2DSlice(name, grid, saveStr):
	#Format
	x = np.arange(grid.shape[1])
	y = np.arange(grid.shape[0])

	X,Y = np.meshgrid(x,y, indexing= 'ij')

	plt.figure()
	plt.contourf(X,Y,grid, 20)

	plt.colorbar()
	plt.title(name)
	plt.savefig("figures/" + saveStr)

	return


# def plotEField(field):
# 	mlab.figure()
# 	im = mlab.quiver3d(field[:,:,:,0], field[:,:,:,1], field[:,:,:,2])
#
# 	mlab.axes()
#
# 	return


#Loading data/Shaving of last dimension
rho = transformData(h5py.File('../../test_rho.grid.h5','r'))
phi = transformData(h5py.File('../../test_phi.grid.h5','r'))
res = transformData(h5py.File('../../test_res.grid.h5','r'))

#Compute ERROR
# error = np.abs((phi - analytical))

# fileE = h5py.File('../test_E.grid.h5', 'r')
# E = fileE['/n=0.0']

# plotPlanesOfGrid("rho", rho)
# plotPlanesOfGrid("phi",phi)
# plotPlanesOfGrid("analytical", analytical)
# plotPlanesOfGrid("error",error)
# plotPlanesOfGrid("res",res)

slice = int(rho.shape[0]*0.3)

# print rho[5,:,:]

plot2DSlice("$\\rho$", rho[slice,:,:], "rho.pdf")
plot2DSlice("Numerical $\phi$", phi[slice,:,:], "numerical.pdf")
plot2DSlice("Residual", res[slice,:,:], "residual.pdf")
# plot2DSlice("Error $|\phi_{num} - \phi_{ana}|$", error[:,:,30], "error.pdf")
# plot2DSlice("Analytical $\phi$", analytical[:,:,30], "analytical.pdf")


# plotEField(E)

# mlab.show()
plt.show()
