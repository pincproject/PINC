# coding: utf-8
import h5py
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fileobj = h5py.File('../../data/object.grid.h5', 'r')
obj = fileobj['Object']
obj = np.transpose(obj,(3,2,1,0))
obj = obj.squeeze()
obj = np.argwhere(obj == 1)
print(obj.shape)
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(obj[:,0],obj[:,1],obj[:,2])
ax.set_xlabel('x axis')
ax.set_ylabel('y axis')
ax.set_zlabel('z axis')
plt.show()