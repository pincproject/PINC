"""
** Transforming PINC output data from HDF5 to VTK format **

Requires: PyEVTK (Python Script to Export VTK)
		  EVTK (Export VTK) package allows exporting data to binary VTK files for visualization
		  and data analysis

@file                HDF52VTK.py
@author              Sayan Adhikari <sayan.adhikari@fys.uio.no>

Instruction for installation:
=============================
git clone https://github.com/paulo-herrera/PyEVTK.git
python3 setup.py install
=============================

Read more about transformation:
https://vtk.org/Wiki/VTK/Writing_VTK_files_using_python
https://github.com/paulo-herrera/PyEVTK
"""
import h5py
import numpy as np
import matplotlib.pyplot as plt
from evtk.hl import pointsToVTK
import os
from time import sleep
from tqdm import tqdm

data_dir = "../../data"#"rhoNeutral" #"P"
file_name = "PINC_particles"
outDir = "animatePointData"

# timesteps:
Nt_start = 100 #1
Nt_stop = 2000 #pos.shape[1]


#========== Data Directory Setup =============
if os.path.exists(outDir):
	os.system('rm '+outDir+'/*')
else:
	os.system('mkdir '+outDir)

#============= Grid denormalization factor ==========
grid = h5py.File(data_dir+'/phi.grid.h5','r')
dimen  = grid.attrs["Axis denormalization factor"][0]
#====================================================
file = h5py.File(data_dir+'/pop.pop.h5','r')
pos0 = file['/pos/specie 0']
vel0 = file['/vel/specie 0']

Np = pos0['n=%.1f'%Nt_start].shape[0]	# Number of particles
Nd = pos0['n=%.1f'%Nt_start].shape[1]	# Number of dimensions

print('number of particles = %i' % Np)
print('number of dimensions = %i' % Nd)

x 	= []
y 	= []
z 	= []
vx 	= []
vy 	= []
vz 	= []
energy 	= []
# for n in range(Nt_start,Nt_stop):
# 	time = str(float(n))
# 	pop = popFile['/pos/specie 0/n=' +time]
# 	x = pop1[:,0]
# 	y = pop1[:,1]
# 	z = pop1[:,2]
# time = int(100)
# time	= str(float(n))
j = 0
for time in tqdm(range(Nt_start,Nt_stop,100)):
	sleep(0.01)
	pos 	= file['/pos/specie 0/n=%d.0'%time]
	vel 	= file['/vel/specie 0/n=%d.5'%time]
	x 	= pos[:,0]*dimen
	y 	= pos[:,1]*dimen
	z 	= pos[:,2]*dimen
	vx 	= vel[:,0]
	vy 	= vel[:,1]
	vz 	= vel[:,2]
	energy = 0.5*(np.square(vx)+np.square(vy)+np.square(vz))
	# pointsToVTK("./"+file_name, x, y, z)
	pointsToVTK("./"+outDir+'/'+file_name+'%d'%j, x, y, z, data = {"energy" : energy})
	j = j+1
