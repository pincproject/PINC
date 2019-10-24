import h5py
import pylab as plt
import numpy as np

hist = h5py.File('../../data/history.xy.h5','r')
pot = hist['/energy/potential/total']
kin = hist['/energy/kinetic/total']
kin0 = hist['/energy/kinetic/specie 0']
kin1 = hist['/energy/kinetic/specie 1']


kin = kin[:,1];		# Extract y-axis
pot = pot[:,1];	# Extract y-axis and invert
tot = pot+kin;		# Collect total energy

kin0 = kin0[:,1];
kin1 = kin1[:,1];

print(len(tot))
avgEn = np.average(tot)
maxEn = np.max(tot)
minEn = np.min(tot)
absError = max(maxEn-avgEn,avgEn-minEn)
relError = absError/avgEn;
print("Relative error: %.2f%%\n"%(relError*100))

#plt.plot(pot,label='potential')
plt.plot(kin,label='kinetic')
plt.plot(tot,label='total')
plt.legend(loc='lower left')
plt.show()

plt.plot(kin0,label='kinetic specie 0')
plt.plot(kin1,label='kinetic specie 1')
plt.legend(loc='lower left')
plt.show()
