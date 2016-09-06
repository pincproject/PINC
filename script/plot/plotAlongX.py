import h5py
import numpy as np
import pylab as plt
from utility import *
import sys as sys

if len(sys.argv) > 1:
    dim = int(sys.argv[1])
    nPlots = int(sys.argv[2])
else:
    dim = 0

for i in range(nPlots):
    n = 0
    m = 0
    f, ax = plt.subplots(3,2)
    for name in ("phi", "sol", "E", "rho", "res", "error"):
        path = '../framework/test_'+name+'_'+ str(i) +'.grid.h5'
        grid = transformData(dim,h5py.File(path,'r'),0)
        plot1DSubgrid(name, grid, ax[n%3,m%2])
        ax[n%3,m%2].set_title(name)
        m+=1
        n+=1-(m%2)


plt.show()
