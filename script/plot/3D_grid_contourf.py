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

def transformData(dataset, timestep):
	grid = dataset['/n=%.1f'%timestep]
	# grid = np.transpose(grid, (3,2,1,0))
	grid = np.squeeze(grid)
	grid = np.average(grid, axis = 0)
	# print grid.shape

	return grid


def plot2DSlice(name, grid, saveStr):
	#Format
	x = np.arange(grid.shape[1])
	y = np.arange(grid.shape[0])

	X,Y = np.meshgrid(x,y)

	plt.figure()
	plt.contourf(X,Y,grid, 20)

	plt.colorbar()
	plt.title(name)
	plt.savefig("figures/" + saveStr)

	return


#Loading data/Shaving of last dimension
rho = transformData(h5py.File('../../test_rho_0.grid.h5','r'),0)
phi = transformData(h5py.File('../../test_phi_0.grid.h5','r'),0)
res = transformData(h5py.File('../../test_res_0.grid.h5','r'),0)
E = transformData(h5py.File('../../test_E_0.grid.h5','r'),0)
sol = transformData(h5py.File('../../test_sol_0.grid.h5','r'),0)


# print E.shape
# exit()

plot2DSlice("$\\rho$", rho, "rho.pdf")
plot2DSlice("Numerical $\phi$", phi, "numerical.pdf")
plot2DSlice("$E$", E[:,:,0], "E.pdf")
plot2DSlice("Residual", res, "residual.pdf")
# plot2DSlice("Error $|\phi_{num} - \phi_{ana}|$", error[:,:,30], "error.pdf")
plot2DSlice("Solution", sol, "analytical.pdf")


# plotEField(E)

# mlab.show()
plt.show()
