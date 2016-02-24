#Some rudimentary plotting of distributions

import h5py
import numpy as np
from mayavi import mlab

#Loading data/Shaving of last dimension
fileRho = h5py.File('../test_rho.grid.h5','r')
rho = fileRho['/n=0.0']
filePhi = h5py.File('../test_phi.grid.h5', 'r')
phi = filePhi['/n=0.0']
fileRes = h5py.File('../test_res.grid.h5', 'r')
res = fileRes['/n=0.0']

rho = rho[:,:,:,0]
phi = phi[:,:,:,0]
res = res[:,:,:,0]


#Plotting
#
# rhoFigure = mlab.figure()
# rhoGrid = mlab.pipeline.scalar_field(rho)
# rhoPlot = mlab.pipeline.volume(rhoGrid)
# mlab.axes()
#
#
# phiFigure = mlab.figure()
# phiGrid = mlab.pipeline.scalar_field(phi)
# phiPlot = mlab.pipeline.volume(phiGrid)
# mlab.axes()
#
# rhoSlice = mlab.figure()
#
im = mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rho),
                            plane_orientation='x_axes',
                            slice_index=rho.shape[0]/2,
                        )

mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rho),
                            plane_orientation='y_axes',
                            slice_index=rho.shape[1]/2,
                        )
mlab.title("rho")
mlab.colorbar(im)
mlab.axes()

phiSlice = mlab.figure()

im = mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(phi),
                            plane_orientation='x_axes',
                            slice_index=phi.shape[0]/2,
                        )

mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(phi),
                            plane_orientation='y_axes',
                            slice_index=phi.shape[1]/2,
                        )
mlab.title("phi")
mlab.colorbar(im)
mlab.axes()

phiSlice = mlab.figure()

im = mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(res),
                            plane_orientation='x_axes',
                            slice_index=res.shape[0]/2,
                        )

mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(res),
                            plane_orientation='y_axes',
                            slice_index=res.shape[1]/2,
                        )
mlab.title("res")
mlab.colorbar(im)
mlab.axes()



mlab.show()
