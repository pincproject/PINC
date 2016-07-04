# Some rudimentary plotting of distributions

import h5py
import numpy as np
from mayavi import mlab

def scatterPlot(pop, charge):
	nParticles = pop.shape[0]

	size = charge*np.ones(nParticles)

	im = mlab.points3d(pop[:,0], pop[:,1], pop[:,2], size, scale_factor=.25)

	return im


mlab.options.offscreen = True

popFile = h5py.File('../../test_pop.pop.h5','r')

print popFile
#Initial to set axes
# mlab.figure()
pop1 = popFile['/pos/specie 0/n=0.0']

#
# plt = scatterPlot(pop1, 5)
#
# mlab.axes()
# v = mlab.view()
#
# for n in range(0,500):
#
# 	time = str(float(n))
#
# 	pop1 = popFile['/pos/specie 0/n=' +time]
# 	# pop2 = popFile['/pos/specie 0/n=' +time]
#
# 	x1 = pop1[:,0]
# 	y1 = pop1[:,1]
# 	z1 = pop1[:,2]
#
# 	plt.mlab_source.set(x=x1, y = y1, z = z1)
#
# 	mlab.savefig("movie/%04d.jpg" % n)
# # 	# mlab.clf()
#
#
#
# mlab.show()
