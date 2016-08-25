#   Plotting a 1D graf of all the grids in a multigrid
#   Used for debugging on a heaviside function
#
#	This presupposes that the different grids are stored numbered
#
#	The wanted # of grids to plot is specified by the input.
#
#	python plotMultigrid1D.py 4
#
#	This will plot 4 grid levels
#
#	by Gullik V Killie


import h5py
import numpy as np
import pylab as plt
import sys as sys

def transformData(dataset, timestep):
	grid = dataset['/n=%.1f'%timestep]
	grid = np.squeeze(grid)
	#Accounting for reverse zyx order
	if len(grid.shape)==4:
		grid = grid[:,:,:,0]
	if grid.shape[0] > 1:
		# grid = np.average(grid, axis = 0)
		grid = grid[1,:,:]
	if grid.shape[1] > 1:
		# grid = np.average(grid, axis = 0)
		grid = grid[1,:]
	return grid

def plot1DSubgrid(name, grid, ax):
	length= grid.shape[0]
	x = np.arange(grid.shape[0])
	ax.plot(x,grid)
	# ax.set_title(str(grid.shape[0]), 'right')
	ax.set_xlim([0,length-1])
	ax.text(length/5, 0, " Max = " + str(np.round(np.max(grid))) + "\n Min = " +
		str(np.round(np.min(grid))) + "\n L=" + str(length))
	ax.locator_params(axis='y',nbins=16)
	ax.locator_params(axis='x',nbins=16)
	ax.grid()

	print name + "\t=" + str(np.max(grid))


def plotAllGrids(name, nLevels, totLevels , savePath = 'figures/'):
	#Plot all grids
	f, ax = plt.subplots(nLevels,1)
	if nLevels == 1:
		path = dataPath + 'test_'+name+'_'+ str(0) +'.grid.h5'
		grid = transformData(h5py.File(path,'r'),0)
		plot1DSubgrid(name, grid, ax)

	else:
		for i in range(nLevels):
			path = dataPath + 'test_'+name+'_'+ str(i) +'.grid.h5'
			grid = transformData(h5py.File(path,'r'),0)
			plot1DSubgrid(name, grid, ax[i])
			plt.setp([a.get_xticklabels() for a in f.axes[:]], visible=False)
			plt.setp([a.get_yticklabels() for a in f.axes[:]], visible=False)
			del grid
	f.suptitle(name)
	f.subplots_adjust(hspace=0)
	f.savefig(savePath + name + str(totLevels) + '.eps')


nLevels = int(sys.argv[1])
if len(sys.argv) > 1:
	totLevels = int(sys.argv[2])
else:
	totLevels = 0

dataPath = "../framework/"
if nLevels ==1:
	n = 0
	m = 0
	f, ax = plt.subplots(3,2)
	for name in ("phi", "sol", "E", "rho", "res", "error"):
		path = dataPath+'test_'+name+'_'+ str(0) +'.grid.h5'
		grid = transformData(h5py.File(path,'r'),0)
		plot1DSubgrid(name, grid, ax[n%3,m%2])
		ax[n%3,m%2].set_title(name)
		m+=1
		n+=1-(m%2)
	# rho = h5py.File('../../test_rho_'+ str(0) +'.grid.h5','r')
	# rho	= rho['/n=%.1f'%0]
	# rho = np.squeeze(rho)
	# print rho


else:
	plotAllGrids("phi", nLevels, totLevels, dataPath)
	plotAllGrids("rho", nLevels, totLevels, dataPath)
	plotAllGrids("res", nLevels, totLevels, dataPath)


plt.show()
