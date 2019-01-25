#!/usr/bin/python
import h5py
import numpy as np
from numpy.linalg import *
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 16})


#decide

path = "../"

min_Vel = 0.
max_Vel = 0.0


# Loading file
file = h5py.File('../pop.pop.h5','r')
test = file['/pos/specie 1']
Nt = len(test) # timesteps

velDenorm = file.attrs.__getitem__("Velocity denormalization factor")

step = 500
n0 = 500.5#1.5
n1 = step*int(Nt/3) + 0.5
n2 = step*int(2*Nt/3) + 0.5
n3 = step*(Nt-1) + 0.5
#n0 = n1


print('%.1f, %.1f, %.1f, %.1f'% (n0, n1, n2, n3))
#pos = file['/pos/specie 1/n=%.1f' % n0]
vel0 = file['/vel/specie 1/n=%.1f' % n0]*velDenorm
vel1 = file['/vel/specie 1/n=%.1f' % n1]*velDenorm
vel2 = file['/vel/specie 1/n=%.1f' % n2]*velDenorm
vel3 = file['/vel/specie 1/n=%.1f' % n3]*velDenorm


print('computing speed specie 1')

Np = vel0.shape[0]	# Number of particles
if vel1.shape[0]<Np: Np=vel1.shape[0]
if vel2.shape[0]<Np: Np=vel2.shape[0]
if vel3.shape[0]<Np: Np=vel3.shape[0]
Nd = vel0.shape[1]	# Number of dimensions	
	
print("number of particles = %i"%Np)

speed0 = np.zeros(Np)
speed0x = np.zeros(Np)
speed0y = np.zeros(Np)
speed0z = np.zeros(Np)
speed1 = np.zeros(Np)
speed1x = np.zeros(Np)
speed1y = np.zeros(Np)
speed1z = np.zeros(Np)
speed2 = np.zeros(Np)
speed2x = np.zeros(Np)
speed2y = np.zeros(Np)
speed2z = np.zeros(Np)
speed3 = np.zeros(Np)
speed3x = np.zeros(Np)
speed3y = np.zeros(Np)
speed3z = np.zeros(Np)

# Compute particle speed

temp0 = vel0[:,:]
temp1 = vel1[:,:]
temp2 = vel2[:,:]
temp3 = vel3[:,:]


for i in range(Np-2):
	speed0[i] = np.linalg.norm(temp0[i,:])
	speed0x[i] = temp0[i,0]
	speed0y[i] = temp0[i,1]
	speed0z[i] = temp0[i,2]

	speed1[i] = np.linalg.norm(temp1[i,:])
	speed1x[i] = temp1[i,0]
	speed1y[i] = temp1[i,1]
	speed1z[i] = temp1[i,2]

	speed2[i] = np.linalg.norm(temp2[i,:])
	speed2x[i] = temp2[i,0]
	speed2y[i] = temp2[i,1]
	speed2z[i] = temp2[i,2]

	speed3[i] = np.linalg.norm(temp3[i,:])
	speed3x[i] = temp3[i,0]
	speed3y[i] = temp3[i,1]
	speed3z[i] = temp3[i,2]

	if speed0[i] > max_Vel: max_Vel = speed0[i]
	if speed1[i] > max_Vel: max_Vel = speed1[i]
	if speed2[i] > max_Vel: max_Vel = speed2[i]
	if speed3[i] > max_Vel: max_Vel = speed3[i] 

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
plt.rc('text', usetex = True)
plt.rc('font', family='serif')


plt.subplot(4,1,1)
plt.subplots_adjust(hspace=0.4)
for label in plt.subplot(4,1,1).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.title(r"Ion  Distribution")
plt.hist(speed0,range =(min_Vel,max_Vel), bins=100, normed=False)
	
plt.subplot(4,1,2)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.hist(speed1,range =(min_Vel,max_Vel), bins=100, normed=False)
for label in plt.subplot(4,1,2).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
	
plt.subplot(4,1,3)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.hist(speed2,range =(min_Vel,max_Vel), bins=100, normed=False)
for label in plt.subplot(4,1,3).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)


plt.ylabel(".\hspace{1cm} Time = %.1f \hspace{1cm} %.1f \hspace{1cm} %.1f  \hspace{2cm} %.1f  "%((n3),(n2),(n1),(n0)))
plt.subplot(4,1,4)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.hist(speed3,range =(min_Vel,max_Vel), bins=100, normed=False)
for label in plt.subplot(4,1,4).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)

plt.xlabel(r"Speed [m/s]")
#plt.show()
plt.savefig(path+"speedIons.png")
plt.close()
plt.clf()

### x,y,z

plt.figure()
plt.subplot(3,1,1)
plt.subplots_adjust(hspace=0.4)
plt.title(r"Ion distributin x,y,z-dimension",)
plt.hist(speed0x,range =(-max_Vel,max_Vel), bins=100, normed=False)

for label in plt.subplot(3,1,1).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.subplot(3,1,2)
plt.hist(speed0y,range =(-max_Vel,max_Vel), bins=100, normed=False)
plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x")

for label in plt.subplot(3,1,2).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.subplot(3,1,3)
plt.hist(speed0z,range =(-max_Vel,max_Vel), bins=100, normed=False)
for label in plt.subplot(3,1,3).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.xlabel(r"Velocity [m/s], Time = %.2f $\displaystyle [dt]$"%(n0))
	
#plt.show()
plt.savefig(path+"speedIons0xyz.png")
plt.close()
plt.clf()

###
plt.figure()
plt.subplot(3,1,1)
plt.subplots_adjust(hspace=0.4)
plt.title(r"Ion distributin x,y,z-dimension",)
plt.hist(speed1x,range =(-max_Vel,max_Vel), bins=100, normed=False)

for label in plt.subplot(3,1,1).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.subplot(3,1,2)
plt.hist(speed1y,range =(-max_Vel,max_Vel), bins=100, normed=False)
plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x")

for label in plt.subplot(3,1,2).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.subplot(3,1,3)
plt.hist(speed1z,range =(-max_Vel,max_Vel), bins=100, normed=False)
for label in plt.subplot(3,1,3).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.xlabel(r"Velocity [m/s], Time = %.2f $\displaystyle [dt]$"%(n1))
	
#plt.show()
plt.savefig(path+"speedIons1xyz.png")
plt.close()
plt.clf()

###

plt.figure()

plt.subplot(3,1,1)
plt.subplots_adjust(hspace=0.4)
plt.title(r"Ion distributin x,y,z-dimension",)
plt.hist(speed2x,range =(-max_Vel,max_Vel), bins=100, normed=False)

for label in plt.subplot(3,1,1).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.subplot(3,1,2)
plt.hist(speed2y,range =(-max_Vel,max_Vel), bins=100, normed=False)
plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x")

for label in plt.subplot(3,1,2).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.subplot(3,1,3)
plt.hist(speed2z,range =(-max_Vel,max_Vel), bins=100, normed=False)
for label in plt.subplot(3,1,3).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.xlabel(r"Velocity [m/s], Time = %.2f $\displaystyle [dt]$"%(n2))
#plt.show()
plt.savefig(path+"speedIons2xyz.png")
plt.close()
plt.clf()

###
plt.figure()
plt.subplot(3,1,1)
plt.subplots_adjust(hspace=0.4)
plt.title(r"Ion distributin x,y,z-dimension",)
plt.hist(speed3x,range =(-max_Vel,max_Vel), bins=100, normed=False)

for label in plt.subplot(3,1,1).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.subplot(3,1,2)
plt.hist(speed3y,range =(-max_Vel,max_Vel), bins=100, normed=False)
plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x")

for label in plt.subplot(3,1,2).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.subplot(3,1,3)
plt.hist(speed3z,range =(-max_Vel,max_Vel), bins=100, normed=False)
for label in plt.subplot(3,1,3).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.xlabel(r"Velocity [m/s], Time = %.2f $\displaystyle [dt]$"%(n3))
#plt.show()
plt.savefig(path+"speedIons3xyz.png")
plt.close()
plt.clf()

	
############################
	
	
#Electrons
	
min_Vel = 0.
max_Vel = 0.0

test = file['/pos/specie 0']

print('%.1f, %.1f, %.1f, %.1f'% (n0, n1, n2, n3))
#pos = file['/pos/specie 1/n=%.1f' % n0]
vel0 = file['/vel/specie 0/n=%.1f' % n0]*velDenorm
vel1 = file['/vel/specie 0/n=%.1f' % n1]*velDenorm
vel2 = file['/vel/specie 0/n=%.1f' % n2]*velDenorm
vel3 = file['/vel/specie 0/n=%.1f' % n3]*velDenorm
	
	
print('computing speed specie 0')
	
Np = vel0.shape[0]	# Number of particles
if vel1.shape[0]<Np: Np=vel1.shape[0]
if vel2.shape[0]<Np: Np=vel2.shape[0]
if vel3.shape[0]<Np: Np=vel3.shape[0]
Nd = vel0.shape[1]	# Number of dimensions
	
print("number of particles = %i"%Np)
	
speed0 = np.zeros(Np)
speed0x = np.zeros(Np)
speed0y = np.zeros(Np)
speed0z = np.zeros(Np)
speed1 = np.zeros(Np)
speed1x = np.zeros(Np)
speed1y = np.zeros(Np)
speed1z = np.zeros(Np)
speed2 = np.zeros(Np)
speed2x = np.zeros(Np)
speed2y = np.zeros(Np)
speed2z = np.zeros(Np)
speed3 = np.zeros(Np)
speed3x = np.zeros(Np)
speed3y = np.zeros(Np)
speed3z = np.zeros(Np)
	
# Compute particle speed
temp0 = vel0[:,:]
temp1 = vel1[:,:]
temp2 = vel2[:,:]
temp3 = vel3[:,:]


for i in range(Np-2):
	speed0[i] = np.linalg.norm(temp0[i,:])
	speed0x[i] = temp0[i,0]
	speed0y[i] = temp0[i,1]
	speed0z[i] = temp0[i,2]

	speed1[i] = np.linalg.norm(temp1[i,:])
	speed1x[i] = temp1[i,0]
	speed1y[i] = temp1[i,1]
	speed1z[i] = temp1[i,2]

	speed2[i] = np.linalg.norm(temp2[i,:])
	speed2x[i] = temp2[i,0]
	speed2y[i] = temp2[i,1]
	speed2z[i] = temp2[i,2]

	speed3[i] = np.linalg.norm(temp3[i,:])
	speed3x[i] = temp3[i,0]
	speed3y[i] = temp3[i,1]
	speed3z[i] = temp3[i,2]

	if speed0[i] > max_Vel: max_Vel = speed0[i]
	if speed1[i] > max_Vel: max_Vel = speed1[i]
	if speed2[i] > max_Vel: max_Vel = speed2[i]
	if speed3[i] > max_Vel: max_Vel = speed3[i] 
	
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
plt.subplots_adjust(hspace=0.4)

plt.title(r"Electron Distribution")
plt.hist(speed0,range =(min_Vel,max_Vel), bins=100, normed=False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
for label in plt.subplot(4,1,1).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)

plt.subplot(4,1,2)
plt.hist(speed1,range =(min_Vel,max_Vel), bins=100, normed=False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
for label in plt.subplot(4,1,2).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)

plt.subplot(4,1,3)
plt.hist(speed2,range =(min_Vel,max_Vel), bins=100, normed=False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.ylabel(".\hspace{1cm} Time = %.1f \hspace{1cm} %.1f \hspace{1cm} %.1f  \hspace{2cm} %.1f  "%((n3),(n2),(n1),(n0)))
for label in plt.subplot(4,1,3).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)

plt.subplot(4,1,4)
plt.hist(speed3,range =(min_Vel,max_Vel), bins=100, normed=False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis

for label in plt.subplot(4,1,4).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)

plt.xlabel(r"Speed [m/s]")
#plt.show()
plt.savefig(path+"speedElectrons.png")
plt.close()

### x,y,z	

plt.figure()
plt.subplot(3,1,1)
plt.subplots_adjust(hspace=0.4)
plt.title(r"Electron distributin x,y,z-dimension",)
plt.hist(speed0x,range =(-max_Vel,max_Vel), bins=100, normed=False)

for label in plt.subplot(3,1,1).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.subplot(3,1,2)
plt.hist(speed0y,range =(-max_Vel,max_Vel), bins=100, normed=False)
plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x")

for label in plt.subplot(3,1,2).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.subplot(3,1,3)
plt.hist(speed0z,range =(-max_Vel,max_Vel), bins=100, normed=False)
for label in plt.subplot(3,1,3).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.xlabel(r"Velocity [m/s], Time = %.2f $\displaystyle [dt]$"%(n0))
#plt.show()
plt.savefig(path+"speedElectrons0xyz.png")
plt.close()
plt.clf()

###

plt.figure()
plt.subplot(3,1,1)
plt.subplots_adjust(hspace=0.4)
plt.title(r"Electron distributin x,y,z-dimension",)
plt.hist(speed1x,range =(-max_Vel,max_Vel), bins=100, normed=False)

for label in plt.subplot(3,1,1).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.subplot(3,1,2)
plt.hist(speed1y,range =(-max_Vel,max_Vel), bins=100, normed=False)
plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x")

for label in plt.subplot(3,1,2).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.subplot(3,1,3)
plt.hist(speed1z,range =(-max_Vel,max_Vel), bins=100, normed=False)
for label in plt.subplot(3,1,3).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.xlabel(r"Velocity [m/s], Time = %.2f $\displaystyle [dt]$"%(n1))
#plt.show()
plt.savefig(path+"speedElectrons1xyz.png")
plt.close()
plt.clf()

###

plt.figure()
plt.subplot(3,1,1)
plt.subplots_adjust(hspace=0.4)
plt.title(r"Electron distributin x,y,z-dimension",)
plt.hist(speed2x,range =(-max_Vel,max_Vel), bins=100, normed=False)

for label in plt.subplot(3,1,1).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.subplot(3,1,2)
plt.hist(speed2y,range =(-max_Vel,max_Vel), bins=100, normed=False)
plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x")

for label in plt.subplot(3,1,2).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.subplot(3,1,3)
plt.hist(speed2z,range =(-max_Vel,max_Vel), bins=100, normed=False)
for label in plt.subplot(3,1,3).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.xlabel(r"Velocity [m/s], Time = %.2f $\displaystyle [dt]$"%(n2))
#plt.show()
plt.savefig(path+"speedElectrons2xyz.png")
plt.close()
plt.clf()

###

plt.figure()
plt.subplot(3,1,1)
plt.subplots_adjust(hspace=0.4)
plt.title(r"Electron distributin x,y,z-dimension",)
plt.hist(speed3x,range =(-max_Vel,max_Vel), bins=100, normed=False)

for label in plt.subplot(3,1,1).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.subplot(3,1,2)
plt.hist(speed3y,range =(-max_Vel,max_Vel), bins=100, normed=False)
plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x")

for label in plt.subplot(3,1,2).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.subplot(3,1,3)
plt.hist(speed3z,range =(-max_Vel,max_Vel), bins=100, normed=False)
for label in plt.subplot(3,1,3).yaxis.get_ticklabels()[:-1]:
    label.set_visible(False)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) # scientific y-axis
plt.xlabel(r"Velocity [m/s], Time = %.2f $\displaystyle [dt]$"%(n3))	

#plt.show()
plt.savefig(path+"speedElectrons3xyz.png")
plt.close()
plt.clf()















