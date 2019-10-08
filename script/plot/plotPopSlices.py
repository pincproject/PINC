# Compares 2 runs in different, used to debug a langmuir wave disrepancy in x-y
# direction

import h5py
import numpy as np
import pylab as plt
import subprocess
from scipy.stats import gaussian_kde


def simplePlot(rho):
	x = np.arange(rho.shape[0])
	y = np.arange(rho.shape[1])

	X,Y = np.meshgrid(x,y,indexing='ij')

	fig, ax = plt.subplots(1)
	im = ax.contourf(X,Y,rho, 50)

	fig.subplots_adjust(bottom = 0.25)
	cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
	fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")

#cmd = "../../mpinc.sh" + " " + "../../input/langmuir2D.ini"

#X perturb
#perturb = " " + "population:perturbAmplitude=0.001,0,0,0"
#perturb += " " + "population:perturbMode=1,0,0,0"
#subprocess.call(cmd + perturb,shell=True)

#fileRhoX = h5py.File('../../data/rho.grid.h5','r')
filePopX = h5py.File('../../data/pop.pop.h5', 'r')
#subprocess.call("rm *.h5",shell=True)

#Y perturb
#perturb = " " + "population:perturbAmplitude=0,0.001,0,0"
#perturb += " " + "population:perturbMode=0,1,0,0"
#subprocess.call(cmd + perturb,shell=True)

#fileRhoY = h5py.File('test_rho.grid.h5','r')
#filePopY = h5py.File('test_pop.pop.h5', 'r')
#subprocess.call("rm *.h5",shell=True)


#rhoX = fileRhoX['/n=4.0']
#rhoX = np.transpose(rhoX,(3,2,1,0))
#rhoX = np.squeeze(rhoX)

#rhoY = fileRhoY['/n=0.0']
#rhoY = np.transpose(rhoY,(0,2,1))
#rhoY = np.squeeze(rhoY)

popX = filePopX['/pos/specie 0/n=4.0']
pop = popX

#popTemp = []

#for i in range(pop.shape[0]):
#	if pop[i,0] < 32 and  pop[i,0] >-1:
#		popTemp.append(pop[i,:])

#pop = np.array(popTemp)

popTemp = []

for i in range(pop.shape[0]):
	if pop[i,2] < 15 and  pop[i,2] >11: # slice 2 to 6
		popTemp.append(pop[i,:])

pop = np.array(popTemp)

x = pop[:,0]
y = pop[:,1]

fig = plt.figure()

plt.hist2d(x, y, bins=200)
plt.colorbar()


#simplePlot(rhoX[:,:,4]) # slice 4
#simplePlot(rhoY)


plt.show()


