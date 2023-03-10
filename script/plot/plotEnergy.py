import h5py
import pylab as plt
import numpy as np
import os

# hist = h5py.File('/mn/fys-server-rp1/simpic/rinkum/two-stream-hb/data/history.xy.h5', 'r')
hist = h5py.File('../../instability/two-stream-hb/data/history.xy.h5', 'r')
# navigate to the directory where the HDF5 file is located
# cd /mn/fys-server-rp1/simpic/rinkum/data/

# open the HDF5 file
# file_path = os.path.join('/mn/fys-server-rp1/simpic', 'rinkum', 'data', 'history.xy.h5')
# print(file_path)
# exit()
# with h5py.File(file_path, 'r') as hist:

pot = hist['/energy/potential/total']
kin = hist['/energy/kinetic/total']
kin0 = hist['/energy/kinetic/specie 0']
kin1 = hist['/energy/kinetic/specie 1']
kin2 = hist['/energy/kinetic/specie 2']
# kin3 = hist['/energy/kinetic/specie 3']
# 

kin = kin[:, 1]		# Extract y-axis
pot = pot[:, 1]	 # Extract y-axis and invert
tot = pot+kin		# Collect total energy

kin0 = kin0[:, 1]
kin1 = kin1[:, 1]
# print(kin1)
# exit()
kin2 = kin2[:, 1]

# kin3 = kin3[:,1]

print(len(tot))
avgEn = np.average(tot)
maxEn = np.max(tot)
minEn = np.min(tot)
absError = max(maxEn-avgEn, avgEn-minEn)
relError = absError/avgEn
print("Relative error: %.2f%%\n" % (relError*100))


# #### FIG SIZE CALC ############

figsize = np.array([200, 200/1.618])  # Figure size in mm
dpi = 300                            # Print resolution
ppi = np.sqrt(1920**2+1200**2)/24    # Screen resolution

# plt.plot(pot,label='potential')
fig, (ax1, ax2,ax3) = plt.subplots(3, 1, figsize=figsize/25.4, constrained_layout=True, dpi=ppi)
ax1.plot(kin, label='kinetic')
ax1.legend(['Kinetic_Energy'])
ax2.plot(pot, label='pot')
ax2.legend(['Potential Energy'])
ax3.plot(tot, label='total')
ax3.legend(['Total Energy'])
# plt.show()
# ax1.legend(loc='lower left')
# ax2.legend(loc='lower left')
# ax3.legend(loc='lower left')


# #### FIG SIZE CALC ############

figsize = np.array([200, 200/1.618])  # Figure size in mm
dpi = 300                             # Print resolution
ppi = np.sqrt(1920**2+1200**2)/24     # Screen resolution

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=figsize/25.4, constrained_layout=True,dpi=ppi)
# ax1.plot(kin0, label='kinetic cold electrons')
ax1.plot(kin1, label='cold electrons')
ax2.plot(kin0, label='precipitating electrons')
ax3.plot(kin2, label='background ions')
# ax4.plot(kin3,label='kinetic specie 3')
ax1.legend(loc='lower left')
ax2.legend(loc='lower left')
ax3.legend(loc='lower left')
# ax4.legend(loc='lower left')
plt.show()
