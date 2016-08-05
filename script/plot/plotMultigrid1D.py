#   Plotting a 1D graf of all the grids in a multigrid
#   Used for debugging on a heaviside function

import h5py
import numpy as np
import pylab as plt

def transformData(dataset, timestep):
	grid = dataset['/n=%.1f'%timestep]
	# grid = np.transpose(grid, (3,2,1,0))
	grid = np.squeeze(grid)
	grid = np.average(grid, axis = 0)
	grid = np.average(grid, axis = 1)
	# grid = grid[20,:,:]

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

def plot1DSubgrid(name, grid):
	x = np.arange(grid.shape[0])

	plt.plot(x,grid)

#Plot all grids
for i in range(3):
	path = '../../test_phi_'+ str(i) +'.grid.h5'
	phi = transformData(h5py.File(path,'r'),0)
	plot1DSubgrid("$\\phi$", phi)

plt.show()
