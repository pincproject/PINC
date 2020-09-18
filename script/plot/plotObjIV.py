import h5py
import pylab as plt
import numpy as np


def calc_EMA(array,S_val,N_val):

    # S_val: Smothing factor
    # N_val: number of point to use in filter

    new_array=np.zeros(len(array))
    k=S_val/(N_val+1)
    new_array[0]=array[0]
    for t in range(1,len(array)):
        new_array[t]=array[t]*k+new_array[t-1]*(1-k)

    return new_array

data= ''

hist = h5py.File('../../data'+data+'/history.xy.h5','r')
curr1 = hist['/current/electrons/dataset']
curr2 = hist['/current/ions/dataset']
pot=hist['/potential/dataset']

#hist2 = h5py.File('../../data_nocoll/history.xy.h5','r')
#curr3 = hist2['/current/electrons/dataset']
#curr4 = hist2['/current/ions/dataset']
#print(curr1)
pot=pot[:,1]
curr1 = curr1[:,1]		# Extract y-axis
curr2 = curr2[:,1];
#curr3 = curr3[:,1]		# Extract y-axis
#curr4 = curr4[:,1];

# Uncoment these to usie Exponential Moving Average
#curr1=calc_EMA(curr1,2,100)
#curr2=calc_EMA(curr2,2,100)
curr=curr1+curr2

#print(len(tot))
#avgEn = np.average(tot)
#maxEn = np.max(tot)
#minEn = np.min(tot)
#absError = max(maxEn-avgEn,avgEn-minEn)
#relError = absError/avgEn;
#print("Relative error: %.2f%%\n"%(relError*100))

plt.plot(pot,label='potential')
#plt.plot(curr,label='combined current')
#plt.plot(pot,curr,label='I-V')

#plt.plot(curr2,label='ion current')
#plt.plot(curr3,label='NOCOLL-electron current')
#plt.plot(curr4,label='NOCOLL-ion current')


plt.legend(loc='lower right')
plt.ylabel("I [A]")
plt.xlabel("V [Volt]")
plt.savefig("current_plots/pot"+data+"both.png")
plt.show()

