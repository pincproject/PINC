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
from evtk.hl import gridToVTK

file_name = "phi"#"rhoNeutral" #"P"

ppc = 128 # particle per cell (for rho plots)

# timesteps:
start = 500# # Must exist in dataset

nx, ny, nz = 32, 32, 32
lx, ly, lz = 0.05, 0.05, 0.05
dx, dy, dz = lx/nx, ly/ny, lz/nz

h5 = h5py.File('data1/'+file_name+'.grid.h5','r') # Edit the directory as per your need

dimen = h5.attrs["Axis denormalization factor"][0]
denorm = h5.attrs["Quantity denormalization factor"][0]

timesteps=[]
for item in h5:
	timesteps.append(item.split("=")[-1])
timesteps = np.array(timesteps,dtype =float)
timesteps = np.sort(timesteps)
end = timesteps[-1]

for i in range(len(timesteps)):
	if timesteps[i] == start:
		timesteps = timesteps[i:]
		break


DATA = []
for i in timesteps:
	dataset = h5["/n=%.1f"%i]
	data = np.squeeze(dataset)
	data = np.transpose(data)*denorm
	if("rho" in file_name ):
		if("rho_e" in file_name ):
			data = -1*data
		if("rho" == file_name ):
			data = (data/denorm)/(ppc) # Number density
		else: data = (data/denorm)/(ppc)

	DATA.append(data)

DATA = np.array(DATA)
avgData = DATA[0]
for i in range(1,len(DATA)):
	avgData += DATA[i]
avgData /= len(DATA)



 # Coordinates
x = np.linspace(0, lx, nx, dtype='float64')
y = np.linspace(0, ly, ny, dtype='float64')
z = np.linspace(0, lz, nz, dtype='float64')

datavtk = avgData.reshape(nx,ny,nz)

gridToVTK("./"+file_name, x, y, z, pointData = {file_name : datavtk})
