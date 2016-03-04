#Some rudimentary plotting of distributions

import h5py
import numpy as np
from mayavi import mlab

def plotPlanesOfGrid(grid):
	mlab.figure()
	im = mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(grid),
	                            plane_orientation='x_axes',
	                            slice_index=grid.shape[0]/2,
	                        )

	mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(grid),
	                            plane_orientation='y_axes',
	                            slice_index=grid.shape[1]/2,
	                        )
	mlab.title("grid")
	mlab.colorbar(im)
	mlab.axes()

	return

def plotEField(field):
	mlab.figure()
	im = mlab.quiver3d(field[:,:,:,0], field[:,:,:,1], field[:,:,:,2])

	mlab.axes()

	return


#Loading data/Shaving of last dimension
filePhi = h5py.File('../test_phi.grid.h5','r')
phi = filePhi['/n=0.0']
phi = phi[:,:,:,0]

fileE = h5py.File('../test_E.grid.h5', 'r')
E = fileE['/n=0.0']


plotPlanesOfGrid(phi)
plotEField(E)

mlab.show()
