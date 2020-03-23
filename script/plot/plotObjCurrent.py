import h5py
import pylab as plt
import numpy as np

hist = h5py.File('../../data/history.xy.h5','r')
curr1 = hist['/current/electrons/dataset']
curr2 = hist['/current/ions/dataset']

hist2 = h5py.File('../../data_nocoll/history.xy.h5','r')
curr3 = hist2['/current/electrons/dataset']
curr4 = hist2['/current/ions/dataset']
print(curr1)

curr1 = curr1[:,1]		# Extract y-axis
curr2 = curr2[:,1];
curr3 = curr3[:,1]		# Extract y-axis
curr4 = curr4[:,1];


#print(len(tot))
#avgEn = np.average(tot)
#maxEn = np.max(tot)
#minEn = np.min(tot)
#absError = max(maxEn-avgEn,avgEn-minEn)
#relError = absError/avgEn;
#print("Relative error: %.2f%%\n"%(relError*100))

#plt.plot(pot,label='potential')
plt.plot(curr1,label='electron current')
plt.plot(curr2,label='ion current')
plt.plot(curr3,label='NOCOLL-electron current')
plt.plot(curr4,label='NOCOLL-ion current')


plt.legend(loc='lower right')
plt.ylabel("I [A]")
plt.xlabel("timesteps")
plt.show()

