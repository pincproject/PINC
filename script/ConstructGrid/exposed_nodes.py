# coding: utf-8
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
points = np.loadtxt('points2.txt', delimiter=',')
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(points[:,0],points[:,1],points[:,2])
ax.set_xlabel('x axis')
ax.set_ylabel('y axis')
ax.set_zlabel('z axis')
plt.show()