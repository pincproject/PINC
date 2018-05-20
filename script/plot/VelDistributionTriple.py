#!/usr/bin/python
import h5py
from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt



#decide

min_Vel = 0.
max_Vel = 0.2


# Loading file
file = h5py.File('../../data/pop.pop.h5','r')
test = file['/pos/specie 1']
Nt = len(test) # timesteps

n0 = 1.5#1.5
#n1 = 1000.5#int(Nt/5) + 0.5
#n2 = 1500.5#int(3*Nt/5) + 0.5
#n3 = 2000.5#Nt - 0.5
n1 = 1000.5
n2 = 2000.5
n3 = 4000.5


print('%.1f, %.1f, %.1f, %.1f'% (n0, n1, n2, n3))
#pos = file['/pos/specie 1/n=%.1f' % n0]
vel0 = file['/vel/specie 1/n=%.1f' % n0]
vel1 = file['/vel/specie 1/n=%.1f' % n1]
vel2 = file['/vel/specie 1/n=%.1f' % n2]
vel3 = file['/vel/specie 1/n=%.1f' % n3]


print('computing speed')

Np = vel0.shape[0]	# Number of particles
Nd = vel0.shape[1]	# Number of dimensions


speed0 = zeros(Np)
speed1 = zeros(Np)
speed2 = zeros(Np)
speed3 = zeros(Np)

# Compute particle speed
for i in range(Np-2):
	speed0[i] = norm(vel0[i,:])
	speed1[i] = norm(vel1[i,:])
	speed2[i] = norm(vel2[i,:])
	speed3[i] = norm(vel3[i,:])

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
plt.hist(speed0,range =(min_Vel,max_Vel), bins=100, normed=True)
#plt.ylim([0,3000])
#plt.xlim([0,0.0012])

plt.subplot(4,1,2)
plt.hist(speed1,range =(min_Vel,max_Vel), bins=100, normed=True)
#plt.ylim([0,3000])
#plt.xlim([0,0.0012])

plt.subplot(4,1,3)
plt.hist(speed2,range =(min_Vel,max_Vel), bins=100, normed=True)
#plt.ylim([0,3000])
#plt.xlim([0,0.0012])
plt.ylabel("Probability Density (part of distr = 1/10^5)")

plt.subplot(4,1,4)
plt.hist(speed3,range =(min_Vel,max_Vel), bins=100, normed=True)
#plt.ylim([0,3000])
#plt.xlim([0,0.0012])

#plt.title("Speed Distribution")

plt.xlabel("Normalized Speed T = %f,%f,%f,%f (*dt)"%(n0,n1,n2,n3))
#plt.show()
plt.savefig("speedIons.png")
plt.close()

############################


#Electrons

min_Vel = 0.
max_Vel = 1.0

test = file['/pos/specie 0']
Nt = len(test) # timesteps

#n0 = 1.5#1.5
#n1 = 1000.5#int(Nt/5) + 0.5
#n2 = 1500.5#int(3*Nt/5) + 0.5
#n3 = 2000.5#Nt - 0.5
#n1 = 1000.5
#n2 = 2000.5
#n3 = 4000.5


print('%.1f, %.1f, %.1f, %.1f'% (n0, n1, n2, n3))
#pos = file['/pos/specie 1/n=%.1f' % n0]
vel0 = file['/vel/specie 0/n=%.1f' % n0]
vel1 = file['/vel/specie 0/n=%.1f' % n1]
vel2 = file['/vel/specie 0/n=%.1f' % n2]
vel3 = file['/vel/specie 0/n=%.1f' % n3]


print('computing speed')

Np = vel0.shape[0]	# Number of particles
Nd = vel0.shape[1]	# Number of dimensions


speed0 = zeros(Np)
speed1 = zeros(Np)
speed2 = zeros(Np)
speed3 = zeros(Np)

# Compute particle speed
for i in range(Np-2):
	speed0[i] = norm(vel0[i,:])
	speed1[i] = norm(vel1[i,:])
	speed2[i] = norm(vel2[i,:])
	speed3[i] = norm(vel3[i,:])

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
plt.hist(speed0,range =(min_Vel,max_Vel), bins=100, normed=True)
#plt.ylim([0,3000])
#plt.xlim([0,0.0012])

plt.subplot(4,1,2)
plt.hist(speed1,range =(min_Vel,max_Vel), bins=100, normed=True)
#plt.ylim([0,3000])
#plt.xlim([0,0.0012])



plt.subplot(4,1,3)
plt.hist(speed2,range =(min_Vel,max_Vel), bins=100, normed=True)
#plt.ylim([0,3000])
#plt.xlim([0,0.0012])

plt.ylabel("Probability Density (part of distr = 1/10^5)")

plt.subplot(4,1,4)
plt.hist(speed3,range =(min_Vel,max_Vel), bins=100, normed=True)
#plt.ylim([0,3000])
#plt.xlim([0,0.0012])

#plt.title("Speed Distribution")

plt.xlabel("Normalized Speed T = %f,%f,%f,%f (*dt)"%(n0,n1,n2,n3))
#plt.show()
plt.savefig("speedElectrons.png")
plt.close()
