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

data_dir = "../../data"#"rhoNeutral" #"P"
file_name = "PINC"

ppc = 30 # particle per cell (for rho plots)

# timesteps:
start = 100# # Must exist in dataset


phih5  = h5py.File(data_dir+'/phi.grid.h5','r') # Edit the directory as per your need
rhoh5  = h5py.File(data_dir+'/rho.grid.h5','r')
rhoeh5 = h5py.File(data_dir+'/rho_e.grid.h5','r')
dimen  = phih5.attrs["Axis denormalization factor"][0]
denormphi = phih5.attrs["Quantity denormalization factor"][0]
denormrho = rhoh5.attrs["Quantity denormalization factor"][0]
denormrhoe = rhoeh5.attrs["Quantity denormalization factor"][0]

timesteps=[]
for item in phih5:
	timesteps.append(item.split("=")[-1])
timesteps = np.array(timesteps,dtype =float)
timesteps = np.sort(timesteps)
end = timesteps[-1]

for i in range(len(timesteps)):
	if timesteps[i] == start:
		timesteps = timesteps[i:]
		break


dataPHI  = []
dataRHO  = []
dataRHOE = []
for i in timesteps:
	datasetPHI  = phih5["/n=%.1f"%i]
	datasetRHO  = rhoh5["/n=%.1f"%i]
	datasetRHOE = rhoeh5["/n=%.1f"%i]
	dataphi = np.squeeze(datasetPHI)
	dataphi = np.transpose(dataphi)*denormphi
	datarho = np.squeeze(datasetRHO)
	datarho = np.transpose(datarho)*denormrho
	datarhoe = np.squeeze(datasetRHOE)
	datarhoe = np.transpose(datarhoe)*denormrhoe
	# if("rho" in file_name ):
	# 	if("rho_e" in file_name ):
	# 		data = -1*data
	# 	if("rho" == file_name ):
	# 		data = (data/denorm)/(ppc) # Number density
	# 	else:
	dataPHI.append(dataphi)
	datarho = (datarho/denormrho)/(ppc)
	dataRHO.append(dataphi)
	datarhoe = -1*datarhoe
	datarhoe = (datarhoe/denormrhoe)/(ppc)
	dataRHOE.append(dataphi)

dataPHI = np.array(dataPHI)
dataRHO = np.array(dataRHO)
dataRHOE = np.array(dataRHOE)
avgPHI = dataPHI[0]
avgRHO = dataRHO[0]
avgRHOE = dataRHOE[0]
for i in range(1,len(dataPHI)):
	avgPHI += dataPHI[i]
	avgRHO += dataRHO[i]
	avgRHOE += dataRHOE[i]
avgPHI /= len(dataPHI)
avgRHO /= len(avgRHO)
avgRHOE /= len(avgRHOE)

#Coordinates
nx, ny, nz = (len(datasetPHI[0,0,:])), (len(datasetPHI[0,:,0])), (len(datasetPHI[:,0,0]))
lx, ly, lz = dimen*(nx-1), dimen*(ny-1), dimen*(nz-1)
x = np.linspace(0, lx, nx, dtype='float64')
y = np.linspace(0, ly, ny, dtype='float64')
z = np.linspace(0, lz, nz, dtype='float64')

phivtk  = avgPHI.reshape(nx,ny,nz)
rhovtk  = avgRHO.reshape(nx,ny,nz)
rhoevtk = avgRHOE.reshape(nx,ny,nz)

gridToVTK("./"+file_name, x, y, z, pointData = {"phi" : phivtk, "rho" : rhovtk, "rho_e" : rhoevtk})
