import h5py
import pylab as plt
import numpy as np


hist = h5py.File('../../data/history.xy.h5','r')
pot = hist['/energy/potential/total']
kin = hist['/energy/kinetic/total']

kin = kin[:,1];		# Extract y-axis
pot = pot[:,1];		# Extract y-axis and invert
tot = pot+kin;		# Collect total energy

avgEn = np.average(tot)
maxEn = np.max(tot)
minEn = np.min(tot)
absError = max(maxEn-avgEn,avgEn-minEn)
relError = absError/avgEn;
print "Relative error: %.2f%%\n"%(relError*100)

plt.plot(pot,label='$E_P$')
plt.plot(kin,label='$E_K$')
plt.plot(tot,label='$E_{Tot}$')
plt.ylabel('Energy')
plt.xlabel('$\Delta t$')
plt.legend(loc='lower right')
plt.savefig('figures/energyPlot.pdf')
plt.show()
