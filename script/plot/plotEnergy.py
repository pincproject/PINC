import h5py
import pylab as plt


hist = h5py.File('../../test_history.xy.h5','r')
pot = hist['/energy/potential/total']
kin = hist['/energy/kinetic/total']

plt.plot(pot[:,1],label='potential')
plt.plot(kin[:,1],label='kinetic')
plt.legend(loc='upper left')
plt.show()
