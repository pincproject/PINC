#Some rudimentary plotting of distributions

#!/usr/bin/python
import h5py
import numpy as np
import pylab as plt
# from mayavi import mlab
#from numba import jit
import sys, os, inspect

#constants
k_b = 1.38064852*10**(-23)
eps_0 = 8.854*10**(-12)

#@jit
def findAverageESquare(input_File,starttime,endtime):
	"""find what timesteps E exists, and return
	array with avg E values"""
	
	denorm = input_File.attrs.__getitem__("Quantity denormalization factor")
	avg_E_square = []
	max_E_square = []
	time = []
	#data = np.transpose(dataset,(3,2,1,0))
	#data = np.squeeze(data)
	for t in range(starttime,endtime):
		val = 0
		max_E = 0
		try: #dont exit on errors!! =  (-.0)
			
			E = np.asarray(input_File['/n=%.1f'%t])
			print(t)
			size = len(E[:,0,0])*len(E[0,:,0])*len(E[0,0,:])
			##print("%f,%f,%f"%(i,j,k))
			for i in range(len(E[:,0,0])):
				for j in range(len(E[0,:,0])):
					for k in range(len(E[0,0,:])):
						temp = np.sqrt(E[i,j,k,0]*E[i,j,k,0]+E[i,j,k,1]*E[i,j,k,1]+E[i,j,k,2]*E[i,j,k,2])
						val += denorm*temp/(size)
						
						if max_E < temp:
							max_E = temp #max E at t
							#print(val)
							
			#each timestep store val
			avg_E_square.append(val)
			max_E_square.append(max_E)
			time.append(t)
		except:
			0
		
	return avg_E_square,max_E_square,time


def plot_electic_field_in_time(starttime,endtime,step,path,dt,Omega_i,Omega_e):
	print("working the Electic field")
	print("ploting %i timesteps, this may take some time"%((endtime-starttime)/step))
	h5file = h5py.File(path +'E.grid.h5','r')

	#compile
	#a,b,c = findAverageESquare(h5file,0,0)

	#run
	avg_E_square,max_E_square,time = findAverageESquare(h5file,starttime,endtime)
	#print(avg_E_square)
	plt.clf()
	time = np.asarray(time)*dt*Omega_e
	plt.plot(time,avg_E_square,label = "averaged E(t)")
	#plt.plot(time,max_E_square,label = "max E^2")
	plt.legend(loc='lower right')
	plt.rc('text', usetex = True) #latex
	plt.rc('font', family='serif')
	plt.xlabel("Time $\displaystyle[\omega_{pi}]$")
	plt.ylabel("V/m")
	plt.savefig((path +"ElectricFieldAvg.png")) #save to file
	plt.plot(time,max_E_square,label = "max E(t)")
	plt.savefig((path +"ElectricFieldMax.png"))
	plt.gcf().clear()
	#plt.show()

def animate(title,path,subdir,h5,startindex,stopindex,step,dt,Omega_i,Omega_e,dx):
	""" makes plot of data perpendicular to B_0 (assumes in z direction)
	in folder "subdir" (relative to "path") at "step" intervals"""
	plt.clf()
	count = startindex	
	for i in range(startindex,stopindex,step):#start and stop timestep
		dataset = h5["/n=%.1f"%i]
		data = np.transpose(dataset,(3,2,1,0))
		data = np.squeeze(data)
		#print(data.shape)
	
		x = np.arange(data.shape[0])
		y = np.arange(data.shape[1])

		X,Y = np.meshgrid(x,y,indexing='ij')
	
		fig, ax = plt.subplots(1)
		im = ax.contourf(X*dx,Y*dx,data[:,:,4], 100)
		
		fig.subplots_adjust(bottom = 0.25)
		plt.rc('text', usetex = True)
		plt.rc('font', family='serif')
		plt.xlabel(r"$\displaystyle\vec{E_0}\times\vec{B_0}$ (x-direction) [m]}",fontsize = 16)
		plt.ylabel(r"$\displaystyle\vec{E_0} $ (y-direction) [m]",fontsize = 16)	
		plt.title(title+r" perpendicular to $\displaystyle\vec{B_0}, t=$ %.2f $\displaystyle[\omega_{pi}]$"%((i*dt*Omega_e)) ,fontsize = 16);

		cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
		fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")
		cbar_rho.axes.tick_params(labelsize = '16') #colorbar
		ax.tick_params(labelsize = '16')	#figure
		plt.savefig((path +subdir +"dt%i.png")%count) #save to file
		count +=1
		#plt.show()
				
def animate_by_name(name,path,starttime,endtime,step,dt,Omega_i,Omega_e,dx):
	"""given "name" of file uses function animate() 
	to open file "name" and write to subdir "name" """
	print("working the "+name+" field")
	h5 = h5py.File(path +name+'.grid.h5','r')
	if not os.path.exists(path + name+"/"):
		os.mkdir(path + name+"/")
	animate(name,path,name+"/",h5,starttime,endtime,step,dt,Omega_i,Omega_e,dx)
		
def plot_energy(path):
	fig = plt.figure()
	plt.rc('text', usetex = True)
	
	#plt.clf()
	hist = h5py.File('../../data/history.xy.h5','r')
	pot = hist['/energy/potential/total']
	kin = hist['/energy/kinetic/total']


	#kin = kin[:,1];		# Extract y-axis
	#pot = pot[:,1];	# Extract y-axis and invert
	#tot = pot+kin;		# Collect total energy

	#specie0
	pot0 = hist['/energy/potential/specie 0']
	kin0 = hist['/energy/kinetic/specie 0']
	kin0 = kin0[:,1];		# Extract y-axis (x is time)

	#specie1
	pot1 = hist['/energy/potential/specie 1']
	kin1 = hist['/energy/kinetic/specie 1']
	kin1 = kin1[:,1];		# Extract y-axis (x is time)


	#plt.plot(pot,label='potential')
	plt.plot(kin0,label='specie0')
	plt.plot(kin1,label='specie1')
	plt.xlabel("t/dt")
	plt.ylabel("Total  Kinetic Energy")
	plt.legend(loc='lower right')
	plt.savefig((path +"Energy_kinetic.png")) #save to file
	plt.legend([])	
	plt.clf()
	plt.cla()
	plt.close()
	fig = plt.figure()
	plt.plot(pot[0],pot[1],label='specie0')
	#plt.plot(pot1,label='specie1')
	plt.xlabel("t/dt")
	plt.ylabel("potential Energy")
	plt.legend(loc='lower right')
	plt.savefig((path +"Energy_potential.png")) #save to file
	plt.gcf().clear()
	#plt.show()

def plot_temperature(path,dt,Omega_i,Omega_e):
	plt.clf()
	plt.close("all")
	fig = plt.figure()
	plt.rc('text', usetex = True)
	plt.rc('font', family='serif')
	hist = h5py.File('../../data/temperature.xy.h5','r')
	#pot = hist['/energy/potential/total']
	kinX = hist['/energy/TemperatureX/specie 1']
	kinY = hist['/energy/TemperatureY/specie 1']
	kinZ = hist['/energy/TemperatureZ/specie 1']
	kinTot1 = hist['/energy/TemperatureTot/specie 1']
	#print(hist['/energy/kinetic/specie 0'])

	kinX = kinX[:,1];		# Extract y-axis
	kinY = kinY[:,1];		# Extract y-axis
	kinZ = kinZ[:,1];		# Extract y-axis
	kinTot1 = kinTot1[:,1];		# Extract y-axis
	time = np.linspace(0,len(kinX),len(kinX))
	time = time*dt*Omega_e
	#print(time) 
	#pot = pot[:,1];	# Extract y-axis and invert
	#tot = pot+kinX;		# Collect total energy
	
	#print(len(tot))
	#avgEn = np.average(tot)
	#maxEn = np.max(tot)
	#minEn = np.min(tot)
	#absError = max(maxEn-avgEn,avgEn-minEn)
	#relError = absError/avgEn;
	#print("Relative error: %.2f%%\n"%(relError*100))
	
	#plt.plot(pot,label='potential')
	#print(kinX[-1])
	plt.plot(time,kinX,label=r'$\displaystyle T_x$')
	plt.plot(time,kinY,label=r'$\displaystyle T_y$')
	plt.plot(time,kinZ,label=r'$\displaystyle T_z$')
	plt.xlabel(r"Time \ $\displaystyle [\omega_{pi}]$")
	#plt.plot(kinTot,label='Temperature tot')
	#plt.plot(tot,label='total')
	plt.legend(loc='center right')
	plt.title('Ions')
	plt.savefig((path +"Temperature_Ions.png"))
	#plt.show()
	plt.close()
	plt.clf()
	fig = plt.figure()


	# Electrons

	#pot = hist['/energy/potential/total']
	kinX = hist['/energy/TemperatureX/specie 0']
	kinY = hist['/energy/TemperatureY/specie 0']
	kinZ = hist['/energy/TemperatureZ/specie 0']
	kinTot0 = hist['/energy/TemperatureTot/specie 0']
	#print(hist['/energy/kinetic/specie 0'])

	kinX = kinX[:,1];		# Extract y-axis
	kinY = kinY[:,1];		# Extract y-axis
	kinZ = kinZ[:,1];		# Extract y-axis
	kinTot0 = kinTot0[:,1];		# Extract y-axis
	#pot = pot[:,1];	# Extract y-axis and invert
	#tot = pot+kinX;		# Collect total energy

	#print(len(tot))
	#avgEn = np.average(tot)
	#maxEn = np.max(tot)
	#minEn = np.min(tot)
	#absError = max(maxEn-avgEn,avgEn-minEn)
	#relError = absError/avgEn;
	#print("Relative error: %.2f%%\n"%(relError*100))

	#plt.plot(pot,label='potential')
	plt.plot(time,kinX,label=r'$\displaystyle T_x$')
	plt.plot(time,kinY,label=r'$\displaystyle T_y$')
	plt.plot(time,kinZ,label=r'$\displaystyle T_z$')
	#plt.plot(kinTot,label='Temperature tot')
	#plt.plot(tot,label='total')
	plt.xlabel(r"Time \ $\displaystyle [\omega_{pi}]$")
	plt.ylabel(r"Temperature \ [K]")
	plt.legend(loc='center right')
	plt.title(r'Electrons')
	plt.savefig((path +"Temperature_Electrons.png"))
	#plt.show()
	plt.close()
	plt.clf()
	fig = plt.figure()

	plt.plot(time,kinTot0,label=r'Temperature \ Electrons')
	plt.plot(time,kinTot1,label=r'Temperature \ Ions')
	#plt.plot(kinTot,label='Temperature tot')
	#plt.plot(tot,label='total')
	plt.xlabel(r"Time \ $\displaystyle[\omega_{pi}]$")
	plt.ylabel(r"Temperature \ [K]")
	plt.legend(loc='center right')
	plt.title(r'Electrons \ vs \ Ions')
	plt.savefig((path +"Temperature_Tot.png"))
	#plt.show()
	plt.close()
	plt.clf()
	#fig = plt.figure()
	
def plot_velocity_distribution(step,path,dt,Omega_pi,Omega_e):
	#decide

	min_Vel = 0.
	max_Vel = 0.0 #found later


	# Loading file
	file = h5py.File(path+'pop.pop.h5','r')
	test = file['/pos/specie 1']
	Nt = len(test) # timesteps

	n0 = 1.5#1.5
	n1 = step*int(Nt/3) + 0.5
	n2 = step*int(2*Nt/3) + 0.5
	n3 = step*(Nt-1) + 0.5
	


	print('%.1f, %.1f, %.1f, %.1f'% (n0, n1, n2, n3))
	#pos = file['/pos/specie 1/n=%.1f' % n0]
	vel0 = file['/vel/specie 1/n=%.1f' % n0]
	vel1 = file['/vel/specie 1/n=%.1f' % n1]
	vel2 = file['/vel/specie 1/n=%.1f' % n2]
	vel3 = file['/vel/specie 1/n=%.1f' % n3]


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
	for i in range(Np-2):
		speed0[i] = np.linalg.norm(vel0[i,:])
		speed0x[i] = vel0[i,0]
		speed0y[i] = vel0[i,1]
		speed0z[i] = vel0[i,2]

		speed1[i] = np.linalg.norm(vel1[i,:])
		speed1x[i] = vel1[i,0]
		speed1y[i] = vel1[i,1]
		speed1z[i] = vel1[i,2]

		speed2[i] = np.linalg.norm(vel2[i,:])
		speed2x[i] = vel2[i,0]
		speed2y[i] = vel2[i,1]
		speed2z[i] = vel2[i,2]

		speed3[i] = np.linalg.norm(vel3[i,:])
		speed3x[i] = vel3[i,0]
		speed3y[i] = vel3[i,1]
		speed3z[i] = vel3[i,2]

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
	plt.title(r"Velocity distribution Ions (part of distr)",fontsize = 16)
	plt.hist(speed0,range =(min_Vel,max_Vel), bins=100, normed=False)
	
	plt.subplot(4,1,2)
	plt.hist(speed1,range =(min_Vel,max_Vel), bins=100, normed=False)
	
	plt.subplot(4,1,3)
	plt.hist(speed2,range =(min_Vel,max_Vel), bins=100, normed=False)
	
	plt.subplot(4,1,4)
	plt.hist(speed3,range =(min_Vel,max_Vel), bins=100, normed=False)
	
	plt.xlabel(r"Normalized Speed T = %.2f,%.2f,%.2f,%.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi),(n1*dt*Omega_pi),(n2*dt*Omega_pi),(n3*dt*Omega_pi)),fontsize = 16)
	#plt.show()
	plt.savefig(path+"speedIons.png")
	plt.close()
	plt.clf()

	### x,y,z

	plt.figure()
	plt.subplot(3,1,1)
	plt.title(r"Velocity distribution Ions (part of distr)",fontsize = 16)
	plt.hist(speed0x,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.subplot(3,1,2)
	plt.hist(speed0y,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x",fontsize = 16)
	plt.subplot(3,1,3)
	plt.hist(speed0z,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.xlabel(r"Normalized Speed x,y,z-dimension T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)),fontsize = 16)
	
	#plt.show()
	plt.savefig(path+"speedIons0xyz.png")
	plt.close()
	plt.clf()

	###
	plt.figure()
	plt.subplot(3,1,1)
	plt.title(r"Velocity distribution Ions (part of distr)",fontsize = 16)
	plt.hist(speed1x,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.subplot(3,1,2)
	plt.hist(speed1y,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x",fontsize = 16)
	plt.subplot(3,1,3)
	plt.hist(speed1z,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.xlabel(r"Normalized Speed x,y,z-dimension T = %.2f $\displaystyle [\omega_{pi}]$"%((n1*dt*Omega_pi)),fontsize = 16)
	
	#plt.show()
	plt.savefig(path+"speedIons1xyz.png")
	plt.close()
	plt.clf()

	###

	plt.figure()
	plt.subplot(3,1,1)
	plt.title(r"Velocity distribution Ions (part of distr)",fontsize = 16)
	plt.hist(speed2x,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.subplot(3,1,2)
	plt.hist(speed2y,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x",fontsize = 16)
	plt.subplot(3,1,3)
	plt.hist(speed2z,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.xlabel(r"Normalized Speed x,y,z-dimension T = %.2f $\displaystyle [\omega_{pi}]$"%((n2*dt*Omega_pi)),fontsize = 16)
	#plt.show()
	plt.savefig(path+"speedIons2xyz.png")
	plt.close()
	plt.clf()

	###
	plt.figure()
	plt.subplot(3,1,1)
	plt.title(r"Velocity distribution Ions (part of distr)",fontsize = 16)
	plt.hist(speed3x,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.subplot(3,1,2)
	plt.hist(speed3y,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x",fontsize = 16)
	plt.subplot(3,1,3)
	plt.hist(speed3z,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.xlabel(r"Normalized Speed x,y,z-dimension T = %.2f $\displaystyle [\omega_{pi}]$"%((n2*dt*Omega_pi)),fontsize = 16)
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
	vel0 = file['/vel/specie 0/n=%.1f' % n0]
	vel1 = file['/vel/specie 0/n=%.1f' % n1]
	vel2 = file['/vel/specie 0/n=%.1f' % n2]
	vel3 = file['/vel/specie 0/n=%.1f' % n3]
	
	
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
	for i in range(Np-2):
		speed0[i] = np.linalg.norm(vel0[i,:])
		speed0x[i] = vel0[i,0]
		speed0y[i] = vel0[i,1]
		speed0z[i] = vel0[i,2]
		
		speed1[i] = np.linalg.norm(vel1[i,:])
		speed1x[i] = vel1[i,0]
		speed1y[i] = vel1[i,1]
		speed1z[i] = vel1[i,2]

		speed2[i] = np.linalg.norm(vel2[i,:])
		speed2x[i] = vel2[i,0]
		speed2y[i] = vel2[i,1]
		speed2z[i] = vel2[i,2]

		speed3[i] = np.linalg.norm(vel3[i,:])
		speed3x[i] = vel3[i,0]
		speed3y[i] = vel3[i,1]
		speed3z[i] = vel3[i,2]
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
	plt.title(r"Velocity distribution Electrons (part of distr)",fontsize = 16)
	plt.hist(speed0,range =(min_Vel,max_Vel), bins=100, normed=False)
	plt.subplot(4,1,2)
	plt.hist(speed1,range =(min_Vel,max_Vel), bins=100, normed=False)
	plt.subplot(4,1,3)
	plt.hist(speed2,range =(min_Vel,max_Vel), bins=100, normed=False)
	plt.ylabel("Probability Density unormalized")
	plt.subplot(4,1,4)
	plt.hist(speed3,range =(min_Vel,max_Vel), bins=100, normed=False)
	plt.xlabel(r"Normalized Speed T = %.2f,%.2f,%.2f,%.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi),(n1*dt*Omega_pi),(n2*dt*Omega_pi),(n3*dt*Omega_pi)),fontsize = 16)
	#plt.show()
	plt.savefig(path+"speedElectrons.png")
	plt.close()

	### x,y,z
	plt.figure()
	plt.subplot(3,1,1)
	plt.title(r"Velocity distribution Electrons x,y,z-dimension (part of distr)",fontsize = 16)
	plt.hist(speed0x,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.subplot(3,1,2)
	plt.hist(speed0y,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x",fontsize = 16)
	plt.subplot(3,1,3)
	plt.hist(speed0z,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.xlabel(r"Normalized Speed T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)),fontsize = 16)
	#plt.show()
	plt.savefig(path+"speedElectrons0xyz.png")
	plt.close()
	plt.clf()

	###
	plt.figure()
	plt.subplot(3,1,1)
	plt.title(r"Velocity distribution Electrons x,y,z-dimension (part of distr)",fontsize = 16)
	plt.hist(speed1x,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.subplot(3,1,2)
	plt.hist(speed1y,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x",fontsize = 16)
	plt.subplot(3,1,3)
	plt.hist(speed1z,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.xlabel(r"Normalized Speed T = %.2f $\displaystyle [\omega_{pi}]$"%((n1*dt*Omega_pi)),fontsize = 16)
	#plt.show()
	plt.savefig(path+"speedElectrons1xyz.png")
	plt.close()
	plt.clf()

	###
	plt.figure()
	plt.subplot(3,1,1)
	plt.title(r"Velocity distribution Electrons x,y,z-dimension (part of distr)",fontsize = 16)
	plt.hist(speed2x,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.subplot(3,1,2)
	plt.hist(speed2y,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x",fontsize = 16)
	plt.subplot(3,1,3)
	plt.hist(speed2z,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.xlabel(r"Normalized Speed T = %.2f $\displaystyle [\omega_{pi}]$"%((n2*dt*Omega_pi)),fontsize = 16)
	#plt.show()
	plt.savefig(path+"speedElectrons2xyz.png")
	plt.close()
	plt.clf()

	###
	plt.figure()
	plt.subplot(3,1,1)
	plt.title(r"Velocity distribution Electrons x,y,z-dimension (part of distr)",fontsize = 16)
	plt.hist(speed3x,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.subplot(3,1,2)
	plt.hist(speed3y,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x",fontsize = 16)
	plt.subplot(3,1,3)
	plt.hist(speed3z,range =(-max_Vel,max_Vel), bins=100, normed=False)
	plt.xlabel(r"Normalized Speed  T = %.2f $\displaystyle [\omega_{pi}]$"%((n2*dt*Omega_pi)),fontsize = 16)	
	#plt.show()
	plt.savefig(path+"speedElectrons3xyz.png")
	plt.close()
	plt.clf()


def plot_vx_vy(step,path,dt,Omega_pi,Omega_e,resolution = 30):
	""" Function plots vx vs vy for both species."""

	# Loading file
	print("ploting Vx-Vy")
	file = h5py.File(path + 'pop.pop.h5','r')
	test = file['/pos/specie 1']
	Nt = len(test) # timesteps

	n0 = 1.5#1.5
	n1 = step*int(Nt/3) + 0.5
	n2 = step*int(2*Nt/3) + 0.5
	n3 = step*(Nt-1) + 0.5

	#print('%.1f, %.1f, %.1f, %.1f'% (n0, n1, n2, n3))
	#pos = file['/pos/specie 1/n=%.1f' % n0]
	vel0 = file['/vel/specie 0/n=%.1f' % n0] # Electrons
	vel1 = file['/vel/specie 0/n=%.1f' % n1] # 
	vel2 = file['/vel/specie 0/n=%.1f' % n2] # 
	vel3 = file['/vel/specie 0/n=%.1f' % n3] # 
	print('putting particles in bins specie 0')

#	Np = vel0.shape[0]	# Number of particles (only adding 1/10000)
	Np = vel0.shape[0]	# Number of particles
	if vel1.shape[0]<Np: Np=vel1.shape[0]
	if vel2.shape[0]<Np: Np=vel2.shape[0]
	if vel3.shape[0]<Np: Np=vel3.shape[0]
	Nd = vel0.shape[1]	# Number of dimensions
	#data = np.zeros([resolution,2])
	speed0 = np.zeros([Np,2])
	speed1 = np.zeros([Np,2])
	speed2 = np.zeros([Np,2])
	speed3 = np.zeros([Np,2])
	max_vel = 0.0#0128739455343 #max(vel0[0,:])
	min_vel = -0.0#0128739455343 #min(vel0[0,:])

	# Compute particle speed
	for i in range(Np):
		speed0[i][0] = vel0[i,0] # v_x
		if max_vel < vel0[i,0]: max_vel = vel0[i,0]
		if min_vel > vel0[i,0]: min_vel = vel0[i,0]
		
		speed1[i][0] = vel1[i,0] # v_x
		if max_vel < vel1[i,0]: max_vel = vel1[i,0]
		if min_vel > vel1[i,0]: min_vel = vel1[i,0]
		
		speed2[i][0] = vel2[i,0] # v_x
		if max_vel < vel2[i,0]: max_vel = vel2[i,0]
		if min_vel > vel2[i,0]: min_vel = vel2[i,0]

		speed3[i][0] = vel3[i,0] # v_x
		if max_vel < vel3[i,0]: max_vel = vel3[i,0]
		if min_vel > vel3[i,0]: min_vel = vel3[i,0]
		

	for i in range(Np):
		speed0[i][1] = vel0[i,1] # v_y
		if max_vel < vel0[i,0]: max_vel = vel0[i,0]
		if min_vel > vel0[i,0]: min_vel = vel0[i,0]

		speed1[i][1] = vel1[i,1] # v_y
		if max_vel < vel1[i,0]: max_vel = vel1[i,0]
		if min_vel > vel1[i,0]: min_vel = vel1[i,0]

		speed2[i][1] = vel2[i,1] # v_y
		if max_vel < vel2[i,0]: max_vel = vel2[i,0]
		if min_vel > vel2[i,0]: min_vel = vel2[i,0]

		speed3[i][1] = vel3[i,1] # v_y
		if max_vel < vel3[i,0]: max_vel = vel3[i,0]
		if min_vel > vel3[i,0]: min_vel = vel3[i,0]

	### set up figure and plot n0
	xedges = []
	yedges = []
	
	for i in range(0,resolution):
		xedges.append(min_vel+((max_vel-min_vel)/resolution)*(i))
		yedges.append(min_vel+((max_vel-min_vel)/resolution)*(i))

	x = speed0[:,0]
	y = speed0[:,1]

	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_vel, max_vel], [min_vel, max_vel]])
	x = np.linspace(min_vel,max_vel,len(H[0,:]))
	y = np.linspace(min_vel,max_vel,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.rc('text', usetex = True) #latex
	plt.rc('font', family='serif')
	plt.xlabel("Vx")
	plt.ylabel("Vy")
	plt.title(r"Electrons  T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)))
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
	#plt.show()




	plt.savefig(path+"Vx_Vy_Electrons0.png")
	plt.clf()	
	
	### set up figure and plot n1
	plt.figure()	

	x = speed1[:,0]
	y = speed1[:,1]
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_vel, max_vel], [min_vel, max_vel]])
	x = np.linspace(min_vel,max_vel,len(H[0,:]))
	y = np.linspace(min_vel,max_vel,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("Vx")
	plt.ylabel("Vy")
	plt.title(r"Electrons  T = %.2f $\displaystyle [\omega_{pi}]$"%((n1*dt*Omega_pi)))
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
	#plt.show()

	plt.savefig(path+"Vx_Vy_Electrons1.png")
	plt.clf()

	### set up figure and plot n2
	plt.figure()	

	x = speed2[:,0]
	y = speed2[:,1]
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_vel, max_vel], [min_vel, max_vel]])
	x = np.linspace(min_vel,max_vel,len(H[0,:]))
	y = np.linspace(min_vel,max_vel,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("Vx")
	plt.ylabel("Vy")
	plt.title(r"Electrons  T = %.2f $\displaystyle [\omega_{pi}]$"%((n2*dt*Omega_pi)))
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
	#plt.show()

	plt.savefig(path+"Vx_Vy_Electrons2.png")
	plt.clf()


### set up figure and plot n3
	plt.figure()	

	x = speed3[:,0]
	y = speed3[:,1]
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_vel, max_vel], [min_vel, max_vel]])
	x = np.linspace(min_vel,max_vel,len(H[0,:]))
	y = np.linspace(min_vel,max_vel,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("Vx")
	plt.ylabel("Vy")
	plt.title(r"Electrons  T = %.2f $\displaystyle [\omega_{pi}]$"%((n3*dt*Omega_pi)))
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
	#plt.show()

	plt.savefig(path+"Vx_Vy_Electrons3.png")
	plt.clf()

	### electrons
	test = file['/pos/specie 1']
	Nt = len(test) # timesteps

	vel0 = file['/vel/specie 1/n=%.1f' % n0] # ions
	vel1 = file['/vel/specie 1/n=%.1f' % n0] #
	vel2 = file['/vel/specie 1/n=%.1f' % n0] #
	vel3 = file['/vel/specie 1/n=%.1f' % n0] #
	print('putting particles in bins specie 1')
	Np = vel0.shape[0]	# Number of particles (only adding 1/10000)
	Nd = vel0.shape[1]	# Number of dimensions


	#data = np.zeros([resolution,2])
	speed0 = np.zeros([Np,2])
	speed1 = np.zeros([Np,2])
	speed2 = np.zeros([Np,2])
	speed3 = np.zeros([Np,2])
	max_vel = 0.0#0128739455343 #max(vel0[0,:])
	min_vel = -0.0#0128739455343 #min(vel0[0,:])

	# Compute particle speed
	for i in range(Np):
		speed0[i][0] = vel0[i,0] # v_x
		if max_vel < vel0[i,0]: max_vel = vel0[i,0]
		if min_vel > vel0[i,0]: min_vel = vel0[i,0]
		
		speed1[i][0] = vel1[i,0] # v_x
		if max_vel < vel1[i,0]: max_vel = vel1[i,0]
		if min_vel > vel1[i,0]: min_vel = vel1[i,0]
		
		speed2[i][0] = vel2[i,0] # v_x
		if max_vel < vel2[i,0]: max_vel = vel2[i,0]
		if min_vel > vel2[i,0]: min_vel = vel2[i,0]

		speed3[i][0] = vel3[i,0] # v_x
		if max_vel < vel3[i,0]: max_vel = vel3[i,0]
		if min_vel > vel3[i,0]: min_vel = vel3[i,0]
		

	for i in range(Np):
		speed0[i][1] = vel0[i,1] # v_y
		if max_vel < vel0[i,0]: max_vel = vel0[i,0]
		if min_vel > vel0[i,0]: min_vel = vel0[i,0]

		speed1[i][1] = vel1[i,1] # v_y
		if max_vel < vel1[i,0]: max_vel = vel1[i,0]
		if min_vel > vel1[i,0]: min_vel = vel1[i,0]

		speed2[i][1] = vel2[i,1] # v_y
		if max_vel < vel2[i,0]: max_vel = vel2[i,0]
		if min_vel > vel2[i,0]: min_vel = vel2[i,0]

		speed3[i][1] = vel3[i,1] # v_y
		if max_vel < vel3[i,0]: max_vel = vel3[i,0]
		if min_vel > vel3[i,0]: min_vel = vel3[i,0]

	xedges = []
	yedges = []
	

	### set up figure and plot n0

	for i in range(0,resolution):
		xedges.append(min_vel+((max_vel-min_vel)/resolution)*(i))
		yedges.append(min_vel+((max_vel-min_vel)/resolution)*(i))

	x = speed0[:,0]
	y = speed0[:,1]
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_vel, max_vel], [min_vel, max_vel]])
	x = np.linspace(min_vel,max_vel,len(H[0,:]))
	y = np.linspace(min_vel,max_vel,len(H[:,0]))
	plt.figure()
	X,Y = np.meshgrid(x,y,indexing='ij')

	#alternative 1
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("Vx")
	plt.ylabel("Vy")
	plt.title(r"Ions  T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)))
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
	#plt.show()

	plt.savefig(path+"Vx_Vy_Ions0.png")
	plt.clf()
	### set up figure and plot n1
	plt.figure()	

	x = speed1[:,0]
	y = speed1[:,1]
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_vel, max_vel], [min_vel, max_vel]])
	x = np.linspace(min_vel,max_vel,len(H[0,:]))
	y = np.linspace(min_vel,max_vel,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("Vx")
	plt.ylabel("Vy")
	plt.title(r"Ions  T = %.2f $\displaystyle [\omega_{pi}]$"%((n1*dt*Omega_pi)))
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
	#plt.show()

	plt.savefig(path+"Vx_Vy_Ions1.png")
	plt.clf()

	### set up figure and plot n2
	plt.figure()	

	x = speed2[:,0]
	y = speed2[:,1]
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_vel, max_vel], [min_vel, max_vel]])
	x = np.linspace(min_vel,max_vel,len(H[0,:]))
	y = np.linspace(min_vel,max_vel,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("Vx")
	plt.ylabel("Vy")
	plt.title(r"Ions  T = %.2f $\displaystyle [\omega_{pi}]$"%((n2*dt*Omega_pi)))
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
	#plt.show()

	plt.savefig(path+"Vx_Vy_Ions2.png")
	plt.clf()

	### set up figure and plot n3
	plt.figure()	

	x = speed3[:,0]
	y = speed3[:,1]
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_vel, max_vel], [min_vel, max_vel]])
	x = np.linspace(min_vel,max_vel,len(H[0,:]))
	y = np.linspace(min_vel,max_vel,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("Vx")
	plt.ylabel("Vy")
	plt.title(r"Ions  T = %.2f $\displaystyle [\omega_{pi}]$"%((n3*dt*Omega_pi)))
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
	#plt.show()

	plt.savefig(path+"Vx_Vy_Ions3.png")
	plt.clf()

def plot_speed_density(dx,dt,nSteps,path,Omega_pi,step,resolution=30):
	
	#alternative 2 
	"""
	from mpl_toolkits.mplot3d import Axes3D
	from matplotlib import cm
	from matplotlib.ticker import LinearLocator, FormatStrFormatter
	"""

	# Loading file
	print("ploting speed density in space")
	file = h5py.File(path+'pop.pop.h5','r')
	test = file['/pos/specie 1']
	Nt = len(test) # timesteps

	n0 = 1.5#1.5
	n1 = step*int(Nt/3) + 0.5
	n2 = step*int(2*Nt/3) + 0.5
	n3 = step*(Nt-1) + 0.5

	#print('%.1f, %.1f, %.1f, %.1f'% (n0, n1, n2, n3))
	#pos = file['/pos/specie 1/n=%.1f' % n0]

	#specie 0
	vel0 = file['/vel/specie 0/n=%.1f' % n0] # Electrons
	vel1 = file['/vel/specie 0/n=%.1f' % n1] # 
	vel2 = file['/vel/specie 0/n=%.1f' % n2] # 
	vel3 = file['/vel/specie 0/n=%.1f' % n3] # 

	pos0 = file['/pos/specie 0/n=%.1f' %(n0-0.5)] # Electrons
	pos1 = file['/pos/specie 0/n=%.1f' %(n1-0.5)] # 
	pos2 = file['/pos/specie 0/n=%.1f' %(n2-0.5)] # 
	pos3 = file['/pos/specie 0/n=%.1f' %(n3-0.5)] #
	print('putting particles in bins specie 0')

#	Np = vel0.shape[0]	# Number of particles (only adding 1/10000)
	Np = vel0.shape[0]	# Number of particles
	if vel1.shape[0]<Np: Np=vel1.shape[0]
	if vel2.shape[0]<Np: Np=vel2.shape[0]
	if vel3.shape[0]<Np: Np=vel3.shape[0]
	Nd = vel0.shape[1]	# Number of dimensions
	data0 = np.zeros([Np,2])
	data1 = np.zeros([Np,2])
	data2 = np.zeros([Np,2])
	data3 = np.zeros([Np,2])
	max_vel = 0.0#0128739455343 #max(vel0[0,:])
	min_vel = -0.0#0128739455343 #min(vel0[0,:])
	max_pos = 0.0
	min_pos = -0.0

	# Compute particle speed
	for i in range(Np):
		data0[i][0] = pos0[i,0] # x
		if max_pos < pos0[i,0]: max_pos = pos0[i,0]
		if min_pos > pos0[i,0]: min_pos = pos0[i,0]
		
		data1[i][0] = pos1[i,0] # x
		if max_pos < pos1[i,0]: max_pos = pos1[i,0]
		if min_pos > pos1[i,0]: min_pos = pos1[i,0]
		
		data2[i][0] = pos2[i,0] # x
		if max_pos < pos2[i,0]: max_pos = pos2[i,0]
		if min_pos > pos2[i,0]: min_pos = pos2[i,0]

		data3[i][0] = pos3[i,0] # x
		if max_pos < pos3[i,0]: max_pos = pos3[i,0]
		if min_pos > pos3[i,0]: min_pos = pos3[i,0]
		

	for i in range(Np):
		data0[i][1] = vel0[i,0] # v_x
		if max_vel < vel0[i,0]: max_vel = vel0[i,0]
		if min_vel > vel0[i,0]: min_vel = vel0[i,0]

		data1[i][1] = vel1[i,0] # v_x
		if max_vel < vel1[i,0]: max_vel = vel1[i,0]
		if min_vel > vel1[i,0]: min_vel = vel1[i,0]

		data2[i][1] = vel2[i,0] # v_x
		if max_vel < vel2[i,0]: max_vel = vel2[i,0]
		if min_vel > vel2[i,0]: min_vel = vel2[i,0]

		data3[i][1] = vel3[i,0] # v_x
		if max_vel < vel3[i,0]: max_vel = vel3[i,0]
		if min_vel > vel3[i,0]: min_vel = vel3[i,0]

	### set up figure and plot n0
	xedges = []
	yedges = []
	
	for i in range(0,resolution):
		xedges.append(min_pos+((max_pos-min_pos)/resolution)*(i))
		yedges.append(min_vel+((max_vel-min_vel)/resolution)*(i))

	x = data0[:,0]
	y = data0[:,1]
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_pos, max_pos], [min_vel, max_vel]])
	x = np.linspace(min_pos,max_pos,len(H[0,:]))
	y = np.linspace(min_vel,max_vel,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1	
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("x")
	plt.ylabel("Vx")
	plt.title(r"Electrons  T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)))
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
	#plt.show()
	plt.rc('text', usetex = True) #latex
	plt.rc('font', family='serif')

	plt.savefig(path+"x_Vx_Electrons0.png")
	plt.clf()

	### set up figure and plot n1
	plt.figure()	
	x = data1[:,0]
	y = data1[:,1]
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_pos, max_pos], [min_vel, max_vel]])
	x = np.linspace(min_pos,max_pos,len(H[0,:]))
	y = np.linspace(min_vel,max_vel,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1	
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("x")
	plt.ylabel("Vx")
	plt.title(r"Electrons  T = %.2f $\displaystyle [\omega_{pi}]$"%((n1*dt*Omega_pi)))
	fig.subplots_adjust(bottom = 0.25)
	cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
	fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")
	#plt.show()

	plt.savefig(path+"x_Vx_Electrons1.png")
	plt.clf()

	### set up figure and plot n2
	plt.figure()	
	x = data2[:,0]
	y = data2[:,1]
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_pos, max_pos], [min_vel, max_vel]])
	x = np.linspace(min_pos,max_pos,len(H[0,:]))
	y = np.linspace(min_vel,max_vel,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1	
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("x")
	plt.ylabel("Vx")
	plt.title(r"Electrons  T = %.2f $\displaystyle [\omega_{pi}]$"%((n2*dt*Omega_pi)))
	fig.subplots_adjust(bottom = 0.25)
	cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
	fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")
	#plt.show()

	plt.savefig(path+"x_Vx_Electrons2.png")
	plt.clf()

	### set up figure and plot n1
	plt.figure()	
	x = data3[:,0]
	y = data3[:,1]
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_pos, max_pos], [min_vel, max_vel]])
	x = np.linspace(min_pos,max_pos,len(H[0,:]))
	y = np.linspace(min_vel,max_vel,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1	
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("x")
	plt.ylabel("Vx")
	plt.title(r"Electrons  T = %.2f $\displaystyle [\omega_{pi}]$"%((n3*dt*Omega_pi)))
	fig.subplots_adjust(bottom = 0.25)
	cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
	fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")
	#plt.show()

	plt.savefig(path+"x_Vx_Electrons3.png")
	plt.clf()

	#specie 1
	vel0 = file['/vel/specie 1/n=%.1f' % n0] # Electrons
	vel1 = file['/vel/specie 1/n=%.1f' % n1] # 
	vel2 = file['/vel/specie 1/n=%.1f' % n2] # 
	vel3 = file['/vel/specie 1/n=%.1f' % n3] # 

	pos0 = file['/pos/specie 1/n=%.1f' %(n0-0.5)] # Electrons
	pos1 = file['/pos/specie 1/n=%.1f' %(n1-0.5)] # 
	pos2 = file['/pos/specie 1/n=%.1f' %(n2-0.5)] # 
	pos3 = file['/pos/specie 1/n=%.1f' %(n3-0.5)] #
	print('putting particles in bins specie 0')

#	Np = vel0.shape[0]	# Number of particles (only adding 1/10000)
	Np = vel0.shape[0]	# Number of particles
	if vel1.shape[0]<Np: Np=vel1.shape[0]
	if vel2.shape[0]<Np: Np=vel2.shape[0]
	if vel3.shape[0]<Np: Np=vel3.shape[0]
	Nd = vel0.shape[1]	# Number of dimensions
	data0 = np.zeros([Np,2])
	data1 = np.zeros([Np,2])
	data2 = np.zeros([Np,2])
	data3 = np.zeros([Np,2])
	max_vel = 0.0#0128739455343 #max(vel0[0,:])
	min_vel = -0.0#0128739455343 #min(vel0[0,:])
	max_pos = 0.0
	min_pos = -0.0

	# Compute particle speed
	for i in range(Np):
		data0[i][0] = pos0[i,0] # x
		if max_pos < pos0[i,0]: max_pos = pos0[i,0]
		if min_pos > pos0[i,0]: min_pos = pos0[i,0]
		
		data1[i][0] = pos1[i,0] # x
		if max_pos < pos1[i,0]: max_pos = pos1[i,0]
		if min_pos > pos1[i,0]: min_pos = pos1[i,0]
		
		data2[i][0] = pos2[i,0] # x
		if max_pos < pos2[i,0]: max_pos = pos2[i,0]
		if min_pos > pos2[i,0]: min_pos = pos2[i,0]

		data3[i][0] = pos3[i,0] # x
		if max_pos < pos3[i,0]: max_pos = pos3[i,0]
		if min_pos > pos3[i,0]: min_pos = pos3[i,0]
		

	for i in range(Np):
		data0[i][1] = vel0[i,0] # v_x
		if max_vel < vel0[i,0]: max_vel = vel0[i,0]
		if min_vel > vel0[i,0]: min_vel = vel0[i,0]

		data1[i][1] = vel1[i,0] # v_x
		if max_vel < vel1[i,0]: max_vel = vel1[i,0]
		if min_vel > vel1[i,0]: min_vel = vel1[i,0]

		data2[i][1] = vel2[i,0] # v_x
		if max_vel < vel2[i,0]: max_vel = vel2[i,0]
		if min_vel > vel2[i,0]: min_vel = vel2[i,0]

		data3[i][1] = vel3[i,0] # v_x
		if max_vel < vel3[i,0]: max_vel = vel3[i,0]
		if min_vel > vel3[i,0]: min_vel = vel3[i,0]

	### set up figure and plot n0
	xedges = []
	yedges = []
	
	for i in range(0,resolution):
		xedges.append(min_pos+((max_pos-min_pos)/resolution)*(i))
		yedges.append(min_vel+((max_vel-min_vel)/resolution)*(i))

	x = data0[:,0]
	y = data0[:,1]
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_pos, max_pos], [min_vel, max_vel]])
	x = np.linspace(min_pos,max_pos,len(H[0,:]))
	y = np.linspace(min_vel,max_vel,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1	
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("x")
	plt.ylabel("Vx")
	plt.title(r"Ions  T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)))
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
	#plt.show()

	plt.savefig(path+"x_Vx_Ions0.png")
	plt.clf()

	### set up figure and plot n1
	plt.figure()	
	x = data1[:,0]
	y = data1[:,1]
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_pos, max_pos], [min_vel, max_vel]])
	x = np.linspace(min_pos,max_pos,len(H[0,:]))
	y = np.linspace(min_vel,max_vel,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1	
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("x")
	plt.ylabel("Vx")
	plt.title(r"Ions  T = %.2f $\displaystyle [\omega_{pi}]$"%((n1*dt*Omega_pi)))
	fig.subplots_adjust(bottom = 0.25)
	cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
	fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")
	#plt.show()

	plt.savefig(path+"x_Vx_Ions1.png")
	plt.clf()

	### set up figure and plot n2
	plt.figure()	
	x = data2[:,0]
	y = data2[:,1]
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_pos, max_pos], [min_vel, max_vel]])
	x = np.linspace(min_pos,max_pos,len(H[0,:]))
	y = np.linspace(min_vel,max_vel,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1	
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("x")
	plt.ylabel("Vx")
	plt.title(r"Ions  T = %.2f $\displaystyle [\omega_{pi}]$"%((n2*dt*Omega_pi)))
	fig.subplots_adjust(bottom = 0.25)
	cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
	fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")
	#plt.show()

	plt.savefig(path+"x_Vx_Ions2.png")
	plt.clf()

	### set up figure and plot n1
	plt.figure()	
	x = data3[:,0]
	y = data3[:,1]
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_pos, max_pos], [min_vel, max_vel]])
	x = np.linspace(min_pos,max_pos,len(H[0,:]))
	y = np.linspace(min_vel,max_vel,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1	
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("x")
	plt.ylabel("Vx")
	plt.title(r"Ions  T = %.2f $\displaystyle [\omega_{pi}]$"%((n3*dt*Omega_pi)))
	fig.subplots_adjust(bottom = 0.25)
	cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
	fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")
	#plt.show()

	plt.savefig(path+"x_Vx_Ions3.png")
	plt.clf()

def plot_kx_ky(dx,dt,path,endtime,Omega_pi,resolution =64):

	probe = h5py.File(path+'probe.xyz.h5','r')
	x = probe['/X']
	y = probe['/Y']
	
	Fs = 1./dx # not Fs but spatial sample rate
	t = np.arange(len(x[:,0]))*dt
	space = np.arange(len(x[0,:]))*dx
	x_0 = x[0,:]		# probe in x-space, time = 0
	y_0 = y[0,:]		# probe in y-space, time = 0
	
	
	n = len(space)
	#y = np.asarray([0,1,0,3,2,5,0,7,0,9])
	#print(y)
	size = n*dx
	kx = space/size # two sides frequency range
	ky = space/size # two sides frequency range
	kx = kx[range(n//2)] # one side frequency range
	ky = ky[range(n//2)]
	
	X_0 = np.fft.fft(x[0,:])/n
	Y_0 = np.fft.fft(y[0,:])/n
	for i in range(1,10000):
		X_0 += np.fft.fft(x[i,:])/n
		Y_0 += np.fft.fft(y[i,:])/n

	X_0 = X_0/10000
	Y_0 = Y_0/10000
	X_0 = X_0[range(n//2)]
	Y_0 = Y_0[range(n//2)]
	
	##plt.plot(space,x_0)
	#plt.show()
	
	#plt.plot(kx,abs(X_0))
	#plt.show()

	### set up figure and plot n0
	
	X_0 = np.real(X_0)
	Y_0 = np.real(Y_0)	
	min_ky = min(Y_0)
	min_kx = min(X_0)
	max_kx = max(X_0)
	max_ky = max(Y_0)	
	
	xedges = []
	yedges = []
	
	for i in range(0,resolution):
		xedges.append(min_kx+((max_kx-min_kx)/resolution)*(i))
		yedges.append(min_ky+((max_ky-min_ky)/resolution)*(i))

	x = X_0
	y = Y_0
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_kx, max_kx], [min_kx, max_ky]])
	x = np.linspace(min_kx,max_kx,len(H[0,:]))
	y = np.linspace(min_ky,max_ky,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1	
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("kx")
	plt.ylabel("ky")
	plt.title("first 1000 timesteps")
	fig.subplots_adjust(bottom = 0.25)
	cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
	fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")

	plt.savefig(path+"kx_ky_0.png")
	plt.clf()


	


	x = probe['/X']
	y = probe['/Y']

	X_0 = np.fft.fft(x[endtime-10000,:])/n
	Y_0 = np.fft.fft(y[endtime-10000,:])/n
	for i in range(int(endtime-9999),int(endtime)):
		X_0 += np.fft.fft(x[i,:])/n
		Y_0 += np.fft.fft(y[i,:])/n

	X_0 = X_0/10000
	Y_0 = Y_0/10000
	X_0 = X_0[range(n//2)]
	Y_0 = Y_0[range(n//2)]

	##plt.plot(space,x_0)
	#plt.show()
	
	#plt.plot(kx,abs(X_0))
	##plt.show()

	### set up figure and plot n0
	
	X_0 = np.real(X_0)
	Y_0 = np.real(Y_0)	
	min_ky = min(Y_0)
	min_kx = min(X_0)
	max_kx = max(X_0)
	max_ky = max(Y_0)	
	
	xedges = []
	yedges = []
	
	for i in range(0,resolution):
		xedges.append(min_kx+((max_kx-min_kx)/resolution)*(i))
		yedges.append(min_ky+((max_ky-min_ky)/resolution)*(i))

	x = X_0
	y = Y_0
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges), range=[[min_kx, max_kx], [min_kx, max_ky]])
	x = np.linspace(min_kx,max_kx,len(H[0,:]))
	y = np.linspace(min_ky,max_ky,len(H[:,0]))
	X,Y = np.meshgrid(x,y,indexing='ij')
	#alternative 1	
	fig, ax = plt.subplots(1)
	im = ax.contourf(x,y,H, resolution)
	plt.xlabel("kx")
	plt.ylabel("ky")
	plt.title("last 1000 timesteps")
	fig.subplots_adjust(bottom = 0.25)
	cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
	fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")

	plt.savefig(path+"kx_ky_last.png")
	plt.clf()








	"""
	import matplotlib.pyplot as plt
	#import plotly.plotly as py
	import numpy as np
	# Learn about API authentication here: https://plot.ly/python/getting-started
	# Find your api_key here: https://plot.ly/settings/api
	
	Fs = 150.0;  # sampling rate
	Ts = 1.0/Fs; # sampling interval
	t = np.arange(0,1,Ts) # time vector
	
	ff = 5;   # frequency of the signal
	y = np.sin(2*np.pi*ff*t)
	
	n = len(y) # length of the signal
	k = np.arange(n)
	T = n/Fs
	frq = k/T # two sides frequency range
	frq = frq[range(n//2)] # one side frequency range
	
	Y = np.fft.fft(y)/n # fft computing and normalization
	Y = Y[range(n//2)]

	fig, ax = plt.subplots(2, 1)
	ax[0].plot(t,y)
	ax[0].set_xlabel('Time')
	ax[0].set_ylabel('Amplitude')
	ax[1].plot(frq,abs(Y),'r') # plotting the spectrum
	ax[1].set_xlabel('Freq (Hz)')
	ax[1].set_ylabel('|Y(freq)|')

	plt.show()

	"""



def post_process_all(starttime,endtime,step,path,dt,dx):
	"""TLWR, kernel "sort of" """

	# DEFAULT VALUES
	nSteps = 64 # number of spatial steps
	q = 1.602e-19 # charge
	B = 7.5e-6 # mag. field
	M_i = 5e-26 #mass Ions
	M_e = 9.109e-31#4e-29 
	n_0 = 1.*10**(9)
	Omega_i = (q*B)/M_i
	Omega_e = (q*B)/M_e
	Omega_pe = np.sqrt((n_0*q*q)/(M_e*eps_0))
	Omega_pi = np.sqrt((n_0*q*q)/(M_i*eps_0))
	Omega_e = Omega_pi #using ion plasma frq
	
	#plot_electic_field_in_time(starttime,endtime,step,path,dt,Omega_i,Omega_e)

	# grid plots
	
	#animate_by_name("rho",path,starttime,endtime,step,dt,Omega_i,Omega_e,dx)
	#animate_by_name("phi",path,starttime,endtime,step,dt,Omega_i,Omega_e,dx)
	
	#plot_energy(path)
	plot_temperature(path,dt,Omega_i,Omega_e)
	
	plot_velocity_distribution(1000,path,dt,Omega_pi,Omega_e)
	
	plot_vx_vy(1000,path,dt,Omega_pi,Omega_e,resolution = 50)

	plot_speed_density(dx,dt,nSteps,path,Omega_pi,1000,resolution = 50) #step = 1000

	plot_kx_ky(dx,dt,path,Omega_pi,nSteps)


post_process_all(0,20000,100,path="../../data/",dt = 4e-7,dx=0.08)


























