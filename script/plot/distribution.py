#!/usr/bin/python
import h5py
from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt
from time import time

start = time()
prev = start

def timer(text):
	global prev
	now = time()
	diff = now-prev
	tot = now-start
	diffunit = "s "
	totunit = "s "
	if tot<1:
		tot *= 1000
		totunit = "ms"
	if tot<1:
		tot *= 1000
		totunit = "us"
	if diff<1:
		diff *= 1000
		diffunit = "ms"
	if diff<1:
		diff *= 1000
		diffunit = "us"
	print "[tot=%6.2f%s, diff=%6.2f%s] %s"%(tot,totunit,diff,diffunit,text)
	prev = time()
	return





# Loading file
file = h5py.File('../../data/pop.pop.h5','r')
pos = file['/pos/specie 1/n=1000.0']
vel = file['/vel/specie 1/n=1000.5']


timer("loading H5 files")

Np = vel.shape[0]	# Number of particles
Nd = vel.shape[1]	# Number of dimensions

speed = zeros(Np)

# Compute particle speed
for i in range(Np):
	speed[i] = norm(vel[i,:])

timer("computig particle speed")

# Analytical distributions
v = linspace(0,6,100)
v2 = linspace(-6,6,100)
vth = 1
dv = 1
if Nd==2: dv = 2*pi*v
if Nd==3: dv = 4*pi*(v**2)
gaussian   = ((1/(2*pi*(vth**2)))**(  0.5 ))*exp(-0.5*(v2/vth)**2)
maxwellian = ((1/(2*pi*(vth**2)))**(Nd/2.0))*exp(-0.5*(v /vth)**2)*dv

timer("computing analytical distributions")

# Plots
plt.figure()
#plt.plot(v,maxwellian)
plt.hist(speed, bins=100, normed=True)
plt.title("Speed Distribution")
plt.xlabel("Normalized Speed")
plt.ylabel("Probability Density")
plt.savefig("speed.png")
plt.close()

plt.figure()
for d in range(Nd):
	plt.subplot(Nd,1,d+1)
	#plt.plot(v2,gaussian)
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

timer("plotting")


