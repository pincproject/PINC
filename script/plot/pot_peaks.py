# coding: utf-8
from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
pot = np.loadtxt('xline.txt',delimiter='\n')
peakIndex, peaks = find_peaks(pot, height=0)
plt.plot(pot)
plt.xlim(0,128)
plt.plot(peakIndex[0], pot[peakIndex[0]], "o")
#plt.plot(np.zeros_like(pot), "--", color="gray")
potMin = -1 * pot
peakIndexMin, peakMin = find_peaks(potMin,height=0)
plt.plot(peakIndexMin[1], pot[peakIndexMin[1]],"o")


plt.annotate(f"{pot[peakIndexMin[1]]:.2f} V", xy=(peakIndexMin[1],pot[peakIndexMin[1]]),  xycoords='data',
            xytext=(peakIndexMin[1]-16,pot[peakIndexMin[1]]-1), textcoords='data')
plt.annotate(f"{pot[peakIndex[0]]:.2f} V", xy=(peakIndex[0],pot[peakIndex[0]]),  xycoords='data',
            xytext=(peakIndex[0]-16,pot[peakIndex[0]]-1), textcoords='data')
#plt.annotate(f"{pot[peakIndexMin[1]]:.2f} V", (peakIndexMin[1],pot[peakIndexMin[1]]))
#plt.annotate(f"{pot[peakIndex[0]]:.2f} V", (peakIndex[0],pot[peakIndex[0]]))
plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.5)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.show()