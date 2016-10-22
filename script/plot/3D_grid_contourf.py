#Some rudimentary plotting of distributions

import h5py
import numpy as np
import pylab as plt
# from mayavi import mlab

import sys, os, inspect

# # Workaround to import from different folder
# cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../plot")))
# if cmd_subfolder not in sys.path:
#     sys.path.insert(0, cmd_subfolder)
#
# from utility import *

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

def transformData(dataset, timestep):
	grid = dataset['/n=%.1f'%timestep]
	# grid = np.transpose(grid, (3,2,1,0))
	grid = np.squeeze(grid)
	grid = np.average(grid, axis = 1)
	# print grid.shape

	return grid


def plot2DSlice(name, grid, saveStr):
	#Format
	x = np.arange(grid.shape[1])
	y = np.arange(grid.shape[0])

	X,Y = np.meshgrid(x,y)

	plt.figure()
	plt.contourf(X,Y,grid, 50)

	plt.colorbar()
	plt.title(name)
	plt.savefig("figures/" + saveStr)

	return

path = "../../data/"
#Loading data/Shaving of last dimension
rho = transformData(h5py.File(path +'rho_0.grid.h5','r'),0)
plot2DSlice("$\\rho$", rho, "rho.pdf")
del rho
res = transformData(h5py.File(path +'res_0.grid.h5','r'),0)
plot2DSlice("Residual", res, "residual.pdf")
del res
E = transformData(h5py.File(path +'E_0.grid.h5','r'),0)
plot2DSlice("$E$", E[:,:,0], "E.pdf")
del E
phi = transformData(h5py.File(path +'phi_0.grid.h5','r'),0)
sol = transformData(h5py.File(path +'sol_0.grid.h5','r'),0)
error = sol - phi
plot2DSlice("Numerical $\phi$", phi, "numerical.pdf")
del phi
plot2DSlice("Solution", sol, "analytical.pdf")
del sol
plot2DSlice("Error $|\phi_{num} - \phi_{ana}|$", error, "error.pdf")
del error



# print E.shape
# exit()



# plotEField(E)

# mlab.show()
plt.show()
