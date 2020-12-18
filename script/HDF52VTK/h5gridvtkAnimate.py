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
import os
from time import sleep
from tqdm import tqdm


data_dir = "../../data"#"rhoNeutral" #"P"
file_name = "PINC"
outDir = "animateGridData"

ppc = 30 # particle per cell (for rho plots)

# timesteps:
start = 100# # Must exist in dataset

# nx, ny, nz = 64, 16, 16
# # lx, ly, lz = 50, 0.05, 0.05
# dx, dy, dz = 0.781250,0.781250,0.781250 #lx/nx, ly/ny, lz/nz
# lx, ly, lz = (nx-1)*dx, (ny-1)*dx, (nz-1)*dx

#========== Data Directory Setup =============
if os.path.exists(outDir):
	os.system('rm '+outDir+'/*.vtr')
else:
	os.system('mkdir '+outDir)

phih5  = h5py.File(data_dir+'/phi.grid.h5','r') # Edit the directory as per your need
rhoh5  = h5py.File(data_dir+'/rho.grid.h5','r')
rhoeh5 = h5py.File(data_dir+'/rho_e.grid.h5','r')
poph5 = h5py.File(data_dir+'/pop.pop.h5','r')
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
	dataPHI.append(dataphi)
	datarho = (datarho/denormrho)/(ppc)
	dataRHO.append(dataphi)
	datarhoe = -1*datarhoe
	datarhoe = (datarhoe/denormrhoe)/(ppc)
	dataRHOE.append(dataphi)

dataPHI = np.array(dataPHI)
dataRHO = np.array(dataRHO)
dataRHOE = np.array(dataRHOE)


#Coordinates
nx, ny, nz = (len(datasetPHI[0,0,:])), (len(datasetPHI[0,:,0])), (len(datasetPHI[:,0,0]))
lx, ly, lz = dimen*(nx-1), dimen*(ny-1), dimen*(nz-1)
x = np.linspace(0, lx, nx, dtype='float64')
# XFlag = x[::-1]
# X = np.require(XFlag, dtype=np.float32, requirements=['A', 'O', 'W', 'C'])
y = np.linspace(0, ly, ny, dtype='float64')
z = np.linspace(0, lz, nz, dtype='float64')

for j in tqdm(range(len(dataPHI))):
	sleep(0.01)
	phivtk  = dataPHI[j].reshape(nx,ny,nz)
	rhovtk  = dataRHO[j].reshape(nx,ny,nz)
	rhoevtk = dataRHOE[j].reshape(nx,ny,nz)
	gridToVTK("./"+outDir+'/'+file_name+'%d'%j, x, y, z, pointData = {"phi" : phivtk, "rho" : rhovtk, "rho_e" : rhoevtk})
