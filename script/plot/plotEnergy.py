import h5py
import pylab as plt


hist = h5py.File('../../test_history.xy.h5','r')
pot = hist['/energy/potential/total']
kin = hist['/energy/kinetic/total']

potmul = -1;

plt.plot(potmul*pot[:,1],label='potential')
plt.plot(kin[:,1],label='kinetic')
plt.plot(potmul*pot[:,1]+kin[:,1],label='total')
plt.legend(loc='lower left')
plt.show()
