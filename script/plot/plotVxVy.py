#!/usr/bin/python
import h5py
from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt
import matplotlib as mpl

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter








# Loading file
file = h5py.File('../../data/pop.pop.h5','r')
test = file['/pos/specie 1']
Nt = len(test) # timesteps
print(test)

n0 = 20000.5


#print('%.1f, %.1f, %.1f, %.1f'% (n0, n1, n2, n3))
#pos = file['/pos/specie 1/n=%.1f' % n0]
vel0 = file['/vel/specie 1/n=%.1f' % n0] # ions



print('putting particles in bins')

Np = vel0.shape[0]	# Number of particles (only adding 1/10000)
Nd = vel0.shape[1]	# Number of dimensions

print(Np)
resolution = 30

data = zeros([resolution,2])


speed0 = zeros([Np,2])

max_vel = 0.00 #max(vel0[0,:])
min_vel = -0.00 #min(vel0[0,:])



# Compute particle speed
for i in range(Np):
	speed0[i][0] = vel0[i,0] # v_x
	if speed0[i][0] > max_vel: max_vel = speed0[i][0]
	if speed0[i][0] < min_vel: min_vel = speed0[i][0]

for i in range(Np):
	speed0[i][1] = vel0[i,1] # v_y
	if speed0[i][1] > max_vel: max_vel = speed0[i][1]
	if speed0[i][1] < min_vel: min_vel = speed0[i][1]


xedges = []
yedges = []

for i in range(0,resolution):
	xedges.append(min_vel+((max_vel-min_vel)/resolution)*(i))
	yedges.append(min_vel+((max_vel-min_vel)/resolution)*(i))

#print(xedges)


x = speed0[:,0]
y = speed0[:,1]


H, xedges, yedges = histogram2d(x, y, bins=(xedges, yedges), range=[[min_vel, max_vel], [min_vel, max_vel]])

#print(xedges)


x = linspace(min_vel,max_vel,len(H[0,:]))
y = linspace(min_vel,max_vel,len(H[:,0]))



#print(H)

X,Y = meshgrid(x,y,indexing='ij')


#print(H[0][:])


#alternative 1
fig, ax = plt.subplots(1)
im = ax.contourf(x,y,H, resolution)
fig.subplots_adjust(bottom = 0.25)
cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")



"""
#alternative 2
fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(X, Y, H, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

fig.colorbar(surf, shrink=0.5, aspect=5)
"""


"""

fig = plt.figure(figsize=(7, 3))
ax = fig.add_subplot(131)
ax.set_title('imshow: equidistant')
im = plt.imshow(H, interpolation='nearest', origin='low',
	extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
"""
"""
ax = fig.add_subplot(132)
ax.set_title('pcolormesh: exact bin edges')
X, Y = meshgrid(xedges, yedges)
ax.pcolormesh(X, Y, H)
ax.set_aspect('equal')

ax = fig.add_subplot(133)
ax.set_title('NonUniformImage: interpolated')
im = mpl.image.NonUniformImage(ax, interpolation='bilinear')
xcenters = xedges[:-1] + 0.5 * (xedges[1:] - xedges[:-1])
ycenters = yedges[:-1] + 0.5 * (yedges[1:] - yedges[:-1])
im.set_data(xcenters, ycenters, H)
ax.images.append(im)
ax.set_xlim(xedges[0], xedges[-1])
ax.set_ylim(yedges[0], yedges[-1])
ax.set_aspect('equal')
"""
plt.show()


"""
plt.figure()
for d in range(Nd):
	plt.subplot(Nd,1,d+1)
	plt.plot(v2,gaussian)
	plt.hist(vel[:,d], bins=100, normed=True)
	plt.title("Velocity Distribution, component %i"%(d))
	plt.ylabel("Probability Density")
plt.xlabel("Normalized Velocity")
plt.savefig("vel.png")
plt.close()

plt.figure()
for d in range(Nd):
	plt.subplot(Nd,1,d+1)
	plt.hist(pos[:,d], bins=100, normed=True)
	plt.title("Position Distribution, component %i"%(d))
	plt.ylabel("Probability Density")
plt.xlabel("Normalized Position")
plt.savefig("pos.png")
plt.close()
"""

