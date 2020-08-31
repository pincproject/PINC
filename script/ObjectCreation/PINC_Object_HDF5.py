#!/usr/bin/python
#Object Creation Using HDF5

import h5py
import numpy as np
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D

##################### INPUT ###########################
# PINC DOMAIN
domainSizeX = 16  #Domain Size in X (Number of Cells)
domainSizeY = 16  #Domain Size in Y (Number of Cells)
domainSizeZ = 16  #Domain Size in Z (Number of Cells)

# PINC OBJECT
OBJ = "sphere"    #Enter object type (sphere or box)
# OBJECT Visualization
visualization = 1

# OBJECT LOCATION (Center for Sphere, Origin for Box)
ObjLocX  = 8
ObjLocY  = 8
ObjLocZ  = 8



if OBJ=="sphere":
    ObjR = 6      #Object radius (Number of Cells)
elif OBJ=="box":
    ObjSizeX = 4  #Object size in X (Number of Cells)
    ObjSizeY = 4  #Object size in Y (Number of Cells)
    ObjSizeZ = 4  #Object size in Z (Number of Cells)
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

    # Interpolating points to the nearest integer
    X = X.astype(int)
    Y = Y.astype(int)
    Z = Z.astype(int)
    # Adding object data to simulation matrix
    for i in range(len(X)):
        for j in range(len(Y)):
            for k in range(len(Z)):
                domain[X[i,j,k],Y[i,j,k],Z[i,j,k]] = 1
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
        x, y, z = np.mgrid[ObjLocX:ObjLocX+ObjSizeX:10j, ObjLocY:ObjLocY+ObjSizeY:10j, ObjLocZ:ObjLocZ+ObjSizeZ:50j]
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
hf = h5py.File('data.grid.h5', 'w')
# Create a dataset named Object and write the domain data
hf.create_dataset('Object', data=domain)
# Close the HDF5 file
hf.close()
