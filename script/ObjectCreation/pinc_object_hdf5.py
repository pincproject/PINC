#!/usr/bin/python

"""
** PINC Object Creation Using HDF5 **

@file                PINC_Object_HDF5.py
@author              Sayan Adhikari <sayan.adhikari@fys.uio.no>
@date                31.08.2020

"""

import h5py
import numpy as np
# import pylab as plt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

##################### INPUT ###########################
# PINC DOMAIN
domainSizeX = 32  #Domain Size in X (Number of Cells)
domainSizeY = 16  #Domain Size in Y (Number of Cells)
domainSizeZ = 16  #Domain Size in Z (Number of Cells)

# PINC OBJECT
OBJ = "sphere"    #Enter object type (sphere or box)
# OBJECT Visualization
visualization = 0

# OBJECT LOCATION (Center for Sphere, Origin for Box)
ObjLocX  = 16 #60 #14 Dirichlet #0: Periodic
ObjLocY  = 8  #1 Dirichlet #0: Periodic
ObjLocZ  = 8  #1 Dirichlet #0: Periodic



if OBJ=="sphere":
    ObjR = 6      #Object radius (Number of Cells)
elif OBJ=="box":
    ObjSizeX = 4 #2 #1  #Object size in X (Number of Cells)  #1: Dirichlet  #1: Periodic
    ObjSizeY = 4 #15 #61  #Object size in Y (Number of Cells) #1: Dirichlet #31: Periodic
    ObjSizeZ = 4 #15 #61  #Object size in Z (Number of Cells) #1: Dirichlet #31: Periodic
else:
    raise Exception("Error in input DATA TYPE")


##################################################################
# Simulation domain matrix
domain = np.zeros([domainSizeZ, domainSizeY, domainSizeX],dtype=int)

if OBJ=="sphere":
    # Spherical Object Creation
    rad, phi, theta = np.mgrid[0:ObjR:50j, 0:2*np.pi:50j, 0:np.pi:50j] #Modify the size to increase the objkect resolution
    X = ObjLocX + rad*np.sin(theta)*np.cos(phi)
    Y = ObjLocX + rad*np.sin(theta)*np.sin(phi)
    Z = ObjLocX + rad*np.cos(theta)

        # Round and cast to int
    X = np.round(X).astype(int)
    Y = np.round(Y).astype(int)
    Z = np.round(Z).astype(int)

    # clamp just in case
    X = np.clip(X, 0, domainSizeX-1)
    Y = np.clip(Y, 0, domainSizeY-1)
    Z = np.clip(Z, 0, domainSizeZ-1)

    # print("X.dtype, Y.dtype, Z.dtype:", X.dtype, Y.dtype, Z.dtype)
    # print("sample index types:", type(X[0,0,0]), type(Y[0,0,0]), type(Z[0,0,0]))

    # fill
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            for k in range(X.shape[2]):
                domain[Z[i,j,k], Y[i,j,k], X[i,j,k]] = 1
elif OBJ=="box":
    domain[ObjLocZ:ObjLocZ+ObjSizeZ,ObjLocY:ObjLocY+ObjSizeY,ObjLocX:ObjLocX+ObjSizeX] = 1
else:
    raise Exception("Error in input DATA TYPE")

###############################################################
# Visualization of the object
if visualization==1:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if OBJ=="sphere":
        ax.scatter(X, Y, Z, c='r', marker='o')
        ax.set_xlim(0,domainSizeX)
        ax.set_ylim(0,domainSizeY)
        ax.set_zlim(0,domainSizeZ)
    elif OBJ=="box":
        x, y, z = np.mgrid[ObjLocX:ObjLocX+ObjSizeX:20j, ObjLocY:ObjLocY+ObjSizeY:20j, ObjLocZ:ObjLocZ+ObjSizeZ:20j]
        x = x.astype(int)
        y = y.astype(int)
        z = z.astype(int)
        ax.scatter(x,y,z, c='r', marker='o')
        ax.set_xlim(0,domainSizeX)
        ax.set_ylim(0,domainSizeY)
        ax.set_zlim(0,domainSizeZ)
    plt.show()
else:
    None
##############################################################
# Open HDF5 file for writing
hf = h5py.File('object.grid.h5', 'w')
# Create a dataset named Object and write the domain data
hf.create_dataset('Object', data=domain)
# Close the HDF5 file
hf.close()
os.system('mv object.grid.h5 data/')
print('Object data copied to data/ \nContents of data directory:')
os.system('ls -lh data/')
