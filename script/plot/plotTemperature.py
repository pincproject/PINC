import h5py
import pylab as plt
import numpy as np

hist = h5py.File('../../data/temperature.xy.h5','r')
#pot = hist['/energy/potential/total']
kinX = hist['/energy/TemperatureX/specie 1']
kinY = hist['/energy/TemperatureY/specie 1']
kinZ = hist['/energy/TemperatureZ/specie 1']
kinTot1 = hist['/energy/TemperatureTot/specie 1']
#print(hist['/energy/kinetic/specie 0'])

kinX = kinX[:,1];		# Extract y-axis
kinY = kinY[:,1];		# Extract y-axis
kinZ = kinZ[:,1];		# Extract y-axis
kinTot1 = kinTot1[:,1];		# Extract y-axis
#pot = pot[:,1];	# Extract y-axis and invert
#tot = pot+kinX;		# Collect total energy

#print(len(tot))
#avgEn = np.average(tot)
#maxEn = np.max(tot)
#minEn = np.min(tot)
#absError = max(maxEn-avgEn,avgEn-minEn)
#relError = absError/avgEn;
#print("Relative error: %.2f%%\n"%(relError*100))

#plt.plot(pot,label='potential')
print(kinX[-1])
plt.plot(kinX,label='Temp x')
plt.plot(kinY,label='Temp y')
plt.plot(kinZ,label='Temp z')
#plt.plot(kinTot,label='Temperature tot')
#plt.plot(tot,label='total')
plt.legend(loc='center right')
plt.title('Ions')
plt.show()


# Electrons

#pot = hist['/energy/potential/total']
kinX = hist['/energy/TemperatureX/specie 0']
kinY = hist['/energy/TemperatureY/specie 0']
kinZ = hist['/energy/TemperatureZ/specie 0']
kinTot0 = hist['/energy/TemperatureTot/specie 0']
#print(hist['/energy/kinetic/specie 0'])

kinX = kinX[:,1];		# Extract y-axis
kinY = kinY[:,1];		# Extract y-axis
kinZ = kinZ[:,1];		# Extract y-axis
kinTot0 = kinTot0[:,1];		# Extract y-axis
#pot = pot[:,1];	# Extract y-axis and invert
#tot = pot+kinX;		# Collect total energy

#print(len(tot))
#avgEn = np.average(tot)
#maxEn = np.max(tot)
#minEn = np.min(tot)
#absError = max(maxEn-avgEn,avgEn-minEn)
#relError = absError/avgEn;
#print("Relative error: %.2f%%\n"%(relError*100))

#plt.plot(pot,label='potential')
plt.plot(kinX,label='Temp x')
plt.plot(kinY,label='Temp y')
plt.plot(kinZ,label='Temp z')
#plt.plot(kinTot,label='Temperature tot')
#plt.plot(tot,label='total')
plt.legend(loc='center right')
plt.title('Electrons')
plt.show()

plt.plot(kinTot0,label='Temperature Electrons')
plt.plot(kinTot1,label='Temperature Ions')
#plt.plot(kinTot,label='Temperature tot')
#plt.plot(tot,label='total')
plt.legend(loc='center right')
plt.title('Electrons vs Ions')
plt.show()
