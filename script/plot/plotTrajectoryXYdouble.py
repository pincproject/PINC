import h5py
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt


file = h5py.File('../../data/pop.pop.h5','r')

pos1 = file['/pos/specie 0']
pos2 = file['/pos/specie 1']

#for item in f.attrs.keys():
#    print item + ":", f.attrs[item]



particleNum = 0
Nt = 100 #pos.shape[1]

Np = pos1['n=1.0'].shape[0]	# Number of particles
Nd = pos1['n=1.0'].shape[1]	# Number of dimensions

print 'number of particles = %i' % Np
print 'number of dimensions = %i' % Nd


x1 = []
y1 = []
z1 = []
x2 = []
y2 = []
z2 = []
for i in range(1,Nt):
	
	positionOfParticle1 = pos1['n=%.1f' % i][particleNum]
	x1.append(positionOfParticle1[:][0])
	y1.append(positionOfParticle1[:][1])
	z1.append(positionOfParticle1[:][2])
	positionOfParticle2 = pos2['n=%.1f' % i][particleNum]
	x2.append(positionOfParticle2[:][0])
	y2.append(positionOfParticle2[:][1])
	z2.append(positionOfParticle2[:][2])



#for i in range(0,Nt-1):
#	print "x1 = %f, y1 = %f z1 = %f" %(x1[i],y1[i],z1[i])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x1, y1, z1, label='Electron')
ax.plot(x2, y2, z2, label='Ion')
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

plt.show()

"""
plt.plot(positionOfParticle[:][2],positionOfParticle[:][1])
plt.title("Position Distribution of particle %i"% particleNum)
plt.ylabel("y/w_pe")
plt.xlabel("x/w_pe")
plt.savefig("pos.png")
plt.close()
"""




