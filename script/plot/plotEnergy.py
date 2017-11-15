import h5py
import pylab as plt
import numpy as np


hist = h5py.File('../../data/history.xy.h5','r')
pot = hist['/energy/potential/total']
kin = hist['/energy/kinetic/total']

kin = kin[:,1];		# Extract y-axis
pot = pot[:,1];	# Extract y-axis and invert
tot = pot+kin;		# Collect total energy

print(len(tot))
avgEn = np.average(tot)
maxEn = np.max(tot)
minEn = np.min(tot)
absError = max(maxEn-avgEn,avgEn-minEn)
relError = absError/avgEn;
print "Relative error: %.2f%%\n"%(relError*100)

#plt.plot(pot,label='potential')
plt.plot(kin,label='kinetic')
plt.plot(tot,label='total')
plt.legend(loc='lower left')
plt.show()
