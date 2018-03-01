#!/usr/bin/python
import h5py
from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt






# Loading file
file = h5py.File('../../data/pop.pop.h5','r')
test = file['/pos/specie 1']
Nt = len(test) # timesteps
print(test)

n0 = 6000.5


#print('%.1f, %.1f, %.1f, %.1f'% (n0, n1, n2, n3))
#pos = file['/pos/specie 1/n=%.1f' % n0]
vel0 = file['/vel/specie 1/n=%.1f' % n0]



print('computing speed')

Np = vel0.shape[0]/10	# Number of particles
Nd = vel0.shape[1]	# Number of dimensions

speed0 = zeros(Np)


# Compute particle speed
for i in range(Np):
	speed0[i] = norm(vel0[i,:])


# Analytical distributions
"""
v = linspace(0,6,100)
v2 = linspace(-6,6,100)
vth = 0.02
dv = 1
if Nd==2: dv = 2*pi*v
if Nd==3: dv = 4*pi*(v**2)
gaussian   = ((1/(2*pi*(vth**2)))**(  0.5 ))*exp(-0.5*(v2/vth)**2)
maxwellian = ((1/(2*pi*(vth**2)))**(Nd/2.0))*exp(-0.5*(v /vth)**2)*dv
"""

print('Plotting')

# Plots
plt.figure()
#plt.plot(v,maxwellian)

plt.subplot(4,1,1)
plt.hist(speed0, bins=1000, normed=True)
#plt.ylim([0,3000])
#plt.xlim([0,0.0012])


#plt.title("Speed Distribution")
plt.ylabel("Probability Density")
plt.xlabel("Normalized Speed")
#plt.show()
plt.savefig("speedIons.png")
plt.close()


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

