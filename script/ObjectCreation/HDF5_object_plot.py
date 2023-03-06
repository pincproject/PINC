#!/usr/bin/python

"""
** PINC Object Visualization **

@file                PINC_object_plot.py
@author              Rinku Mishra
@date                22.01.2023

"""

import h5py
import numpy as np
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

if len(sys.argv)<2:
    raise Exception("Provide the path to object file")
else:
    obj_file = sys.argv[1]

hf = h5py.File(obj_file, 'r')

xdim, ydim, zdim = hf['Object'].shape

data = hf['Object'][:,:,:]
index_obj = np.where(data==1)

x, y, z = index_obj[0], index_obj[1], index_obj[2]
uni_x = list(dict.fromkeys(x))
uni_y = list(dict.fromkeys(y))
uni_z = list(dict.fromkeys(z))
# print(uni_x,uni_y,uni_z)
print("Object size: \n",len(uni_x),len(uni_y),len(uni_z)," in terms of Debye length.")
print("NOTE: In PINC the object is read in reverse order.")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x,y,z, c='r', marker='o')
ax.set_xlim(0,xdim)
ax.set_ylim(0,ydim)
ax.set_zlim(0,zdim)
ax.set_xlabel("x")
ax.set_xlabel("y")
ax.set_xlabel("z")
plt.show()

