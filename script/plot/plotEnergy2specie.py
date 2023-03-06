import h5py
import pylab as plt
import numpy as np


hist = h5py.File('../../instability/two-stream/data/history.xy.h5','r')
pot = hist['/energy/potential/total']
kin = hist['/energy/kinetic/total']


#kin = kin[:,1];		# Extract y-axis
#pot = pot[:,1];	# Extract y-axis and invert
#tot = pot+kin;		# Collect total energy

#specie0
pot0 = hist['/energy/potential/specie 0']
kin0 = hist['/energy/kinetic/specie 0']
kin0 = kin0[:,1];		# Extract y-axis (x is time)

#specie1
pot1 = hist['/energy/potential/specie 1']
kin1 = hist['/energy/kinetic/specie 1']
kin1 = kin1[:,1];		# Extract y-axis (x is time)


#plt.plot(pot,label='potential')
plt.plot(kin0,label='specie0')
plt.plot(kin1,label='specie1')
plt.xlabel("t/dt")
plt.ylabel("Total  Kinetic Energy")
plt.legend(loc='lower left')
plt.show()
