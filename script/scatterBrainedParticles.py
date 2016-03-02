#Some rudimentary plotting of distributions

import h5py
import numpy as np
from mayavi import mlab

popFile = h5py.File('../test_pop.pop.h5','r')
pop1 = popFile['/pos/specie 0/n=0.0']

nParticles = pop1.shape[0]

size = 1* np.ones(nParticles)

mlab.figure()
mlab.points3d(pop1[:,0], pop1[:,1], pop1[:,2], size, scale_factor=.25)
mlab.axes()
mlab.show()
