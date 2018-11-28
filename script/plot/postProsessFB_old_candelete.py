#Some rudimentary plotting of distributions

#!/usr/bin/python
import h5py
import numpy as np
import pylab as plt
from matplotlib import colors, ticker, cm
# from mayavi import mlab
#from numba import jit
import sys, os, inspect
from numpy import inf

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
		im = ax.contourf(X*dx,Y*dx,data[:,4,:], 100)

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


def plot_velocity_distribution(step,path,dt,Omega_pi,Omega_e,mem_step=100000):
	#Handles large datasets, use memstep*3*float as a guide on memory usage.

	min_Vel = 0.
	max_Vel = 0.0 #found later


	# Loading file
	file = h5py.File(path+'pop.pop.h5','r')
	test = file['/pos/specie 1']
	Nt = len(test) # timesteps

	n0 = step +0.5#1.5


	#pos = file['/pos/specie 1/n=%.1f' % n0]
	vel0 = file['/vel/specie 1/n=%.1f' % n0]


	print('computing speed specie 1')

	Np = int(vel0.shape[0]/part_distr)	# Number of particles
	Nd = vel0.shape[1]	# Number of dimensions

	print("number of particles = %i"%Np)
	print("ploting %i at a time"%(mem_step))


	mem_step = int(mem_step)
	Np = int(Np)
	real_index = 0 # index of real distribution
	#counter = 0 # index of local distribution
	iterations = 1 # to get percentage done


	numbins = 1000
	real_hist = np.zeros(numbins-1, dtype='int32')

	## first sub distr
	speed0 = np.zeros(mem_step)
	temp_vel = np.zeros((mem_step,3),dtype='float64')
	print("computing speeds %f %% done"%(((float(100*0*mem_step))/Np)))
	temp_vel = vel0[real_index:mem_step,:]

	for counter in range(0,int(mem_step)): #could be for loop

		# Compute particle speed
		temp = np.linalg.norm(temp_vel[counter,:])

		speed0[counter] = temp#np.linalg.norm(temp_vel[counter,:])
		#print(speed0[counter])
		if real_index < mem_step:
			if speed0[counter] > max_Vel: max_Vel = speed0[counter] + speed0[counter]/10.

		#counter += 1
		real_index += 1

	bins = np.linspace(min_Vel, max_Vel, numbins)

	## sucsessive sub distributions
	while real_index < int(Np-mem_step):
		print("computing speeds %f %% done"%(((float(100*iterations*mem_step))/Np)))
		#counter = 0
		#speed0 = np.zeros(mem_step)
		#speed0x = np.zeros(mem_step)
		#speed0y = np.zeros(mem_step)
		#speed0z = np.zeros(mem_step)

		#temp_vel = np.zeros((mem_step,3),dtype='float64')
		temp_vel = vel0[real_index:(iterations+1)*mem_step,:]
		for counter in range(0,int(mem_step)): #could be for loop

			# Compute particle speed
			temp = np.linalg.norm(temp_vel[counter,:])
			speed0[counter] = temp#np.linalg.norm(temp_vel[counter,:])
			#counter += 1
			real_index += 1
		iterations +=1

		htemp, jnk = np.histogram(speed0, bins)
    	#real_hist+=htemp
		real_hist += htemp
		#print(real_hist)

	## LAST asymetric part
	temp_vel = vel0[real_index:Np,:]

	for counter in range(0,Np-real_index): #could be for loop

		# Compute particle speed
		temp = np.linalg.norm(temp_vel[counter,:])
		speed0[counter] = temp#np.linalg.norm(temp_vel[counter,:])
		#counter += 1
		real_index += 1
	iterations +=1

	htemp, jnk = np.histogram(speed0, bins)
	#real_hist+=htemp
	real_hist += htemp
	#print(real_hist)


	# Plots
	plt.figure()
	plt.subplot(1,1,1)
	plt.rc('text', usetex = True)
	plt.rc('font', family='serif')
	#plt.hist(real_hist,range =jnk, bins=numbins, normed=False)
	plt.plot(np.linspace(min_Vel,max_Vel,numbins-1),real_hist)
	plt.title(r"Velocity distribution Ions (part of distr)",fontsize = 16)
	plt.xlabel(r"Normalized Speed T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)),fontsize = 16)
	#plt.show()
	plt.savefig(path+"speedIons.png")
	plt.close()
	plt.clf()

	### x,y,z
	min_Vel = 0.
	max_Vel = 0.0 #found later
	real_index = 0 # index of real distribution
	#counter = 0 # index of local distribution
	iterations = 1 # to get percentage done

	real_histX = np.zeros(numbins-1, dtype='int32')
	real_histY = np.zeros(numbins-1, dtype='int32')
	real_histZ = np.zeros(numbins-1, dtype='int32')


	## first sub distr
	temp_velX = np.zeros((mem_step),dtype='float64')
	temp_velY = np.zeros((mem_step),dtype='float64')
	temp_velZ = np.zeros((mem_step),dtype='float64')
	temp_velX = vel0[real_index:mem_step,0]
	temp_velY = vel0[real_index:mem_step,1]
	temp_velZ = vel0[real_index:mem_step,2]
	speed0x = np.zeros(mem_step)
	speed0y = np.zeros(mem_step)
	speed0z = np.zeros(mem_step)
	for counter in range(0,int(mem_step)):
		#Compute particle speed
		speed0x[counter] = temp_velX[counter] #vel0[real_index,0]
		speed0y[counter] = temp_velY[counter] #vel0[real_index,1]
		speed0z[counter] = temp_velZ[counter] #vel0[real_index,2]
		if real_index < mem_step:
			if speed0x[counter] > max_Vel: max_Vel = speed0x[counter] + speed0x[counter]/10.
			if speed0y[counter] > max_Vel: max_Vel = speed0y[counter] + speed0y[counter]/10.
			if speed0z[counter] > max_Vel: max_Vel = speed0z[counter] + speed0z[counter]/10.

			if speed0x[counter] < min_Vel: min_Vel = speed0x[counter] - speed0x[counter]/10.
			if speed0y[counter] < min_Vel: min_Vel = speed0y[counter] - speed0y[counter]/10.
			if speed0z[counter] < min_Vel: min_Vel = speed0z[counter] - speed0z[counter]/10.
		counter += 1
		real_index += 1

	bins = np.linspace(min_Vel, max_Vel, numbins)
	## sucsessive sub distributions
	while real_index < Np-mem_step:
		print("computing velocities %i %% done"%(((float(100*iterations*mem_step))/Np)))
		counter = 0
		#speed0x = np.zeros(mem_step)
		#speed0y = np.zeros(mem_step)
		#speed0z = np.zeros(mem_step)

		temp_velX = vel0[real_index:(iterations+1)*mem_step,0]
		temp_velY = vel0[real_index:(iterations+1)*mem_step,1]
		temp_velZ = vel0[real_index:(iterations+1)*mem_step,2]
		for counter in range(0,int(mem_step)): #could be for loop
			#Compute particle speed
			speed0x[counter] = temp_velX[counter]#vel0[real_index,0]
			speed0y[counter] = temp_velY[counter]#vel0[real_index,1]
			speed0z[counter] = temp_velZ[counter]#vel0[real_index,2]
			real_index += 1
		iterations +=1
		htempX, jnk = np.histogram(speed0x, bins)
		htempY, jnk = np.histogram(speed0y, bins)
		htempZ, jnk = np.histogram(speed0z, bins)
		real_histX += htempX
		real_histY += htempY
		real_histZ += htempZ
		#print(real_hist)


	## LAST asymetric part
	temp_vel = vel0[real_index:Np,:]

	temp_velX = vel0[real_index:Np,0]
	temp_velY = vel0[real_index:Np,1]
	temp_velZ = vel0[real_index:Np,2]

	for counter in range(0,Np-real_index): #could be for loop

		# Compute particle speed
		speed0x[counter] = temp_velX[counter]#vel0[real_index,0]
		speed0y[counter] = temp_velY[counter]#vel0[real_index,1]
		speed0z[counter] = temp_velZ[counter]#vel0[real_index,2]
		#counter += 1
		real_index += 1
	iterations +=1

	htempX, jnk = np.histogram(speed0x, bins)
	htempY, jnk = np.histogram(speed0y, bins)
	htempZ, jnk = np.histogram(speed0z, bins)
	#real_hist+=htemp
	real_histX += htempX
	real_histY += htempY
	real_histZ += htempZ
	#print(real_hist)

	plt.figure()
	plt.subplot(3,1,1)
	plt.title(r"Velocity distribution Ions (part of distr)",fontsize = 16)
	plt.plot(np.linspace(-max_Vel,max_Vel,numbins-1),real_histX)
	plt.subplot(3,1,2)
	plt.plot(np.linspace(-max_Vel,max_Vel,numbins-1), real_histY)
	plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x",fontsize = 16)
	plt.subplot(3,1,3)
	plt.plot(np.linspace(-max_Vel,max_Vel,numbins-1),real_histZ)
	plt.xlabel(r"Normalized Speed x,y,z-dimension T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)),fontsize = 16)

	#plt.show()
	plt.savefig(path+"speedIons0xyz.png")
	plt.close()
	plt.clf()

	############################


	#Electrons

	min_Vel = 0.
	max_Vel = 0.0

	test = file['/pos/specie 0']


	#pos = file['/pos/specie 1/n=%.1f' % n0]
	vel0 = file['/vel/specie 0/n=%.1f' % n0]


	print('computing speed specie 0')
	Np = vel0.shape[0]/part_distr	# Number of particles
	Nd = vel0.shape[1]	# Number of dimensions
	Np = int(Np)
	print("number of particles = %i"%Np)
	print("ploting %i at a time"%(mem_step))



	real_index = 0 # index of real distribution
	#counter = 0 # index of local distribution
	iterations = 1 # to get percentage done


	numbins = 1000
	real_hist = np.zeros(numbins-1, dtype='int32')
	## first sub distr
	speed0 = np.zeros(mem_step)
	temp_vel = np.zeros((mem_step,3),dtype='float64')
	print("computing speeds %f %% done"%(((float(100*0*mem_step))/Np)))
	temp_vel = vel0[real_index:mem_step,:]

	for counter in range(0,int(mem_step)): #could be for loop

		# Compute particle speed
		temp = np.linalg.norm(temp_vel[counter,:])

		speed0[counter] = temp#np.linalg.norm(temp_vel[counter,:])
		#print(speed0[counter])
		if real_index < mem_step:
			if speed0[counter] > max_Vel: max_Vel = speed0[counter] + speed0[counter]/10.

		#counter += 1
		real_index += 1

	bins = np.linspace(min_Vel, max_Vel, numbins)

	## sucsessive sub distributions
	while real_index < int(Np-mem_step):
		print("computing speeds %f %% done"%(((float(100*iterations*mem_step))/Np)))
		#counter = 0
		#speed0 = np.zeros(mem_step)
		#speed0x = np.zeros(mem_step)
		#speed0y = np.zeros(mem_step)
		#speed0z = np.zeros(mem_step)

		#temp_vel = np.zeros((mem_step,3),dtype='float64')
		temp_vel = vel0[real_index:(iterations+1)*mem_step,:]
		for counter in range(0,int(mem_step)): #could be for loop

			# Compute particle speed
			temp = np.linalg.norm(temp_vel[counter,:])
			speed0[counter] = temp#np.linalg.norm(temp_vel[counter,:])
			#counter += 1
			real_index += 1
		iterations +=1

		htemp, jnk = np.histogram(speed0, bins)
    	#real_hist+=htemp
		real_hist += htemp
		#print(real_hist)

	## LAST asymetric part
	temp_vel = vel0[real_index:Np,:]
	for counter in range(0,Np-real_index): #could be for loop

		# Compute particle speed
		temp = np.linalg.norm(temp_vel[counter,:])
		speed0[counter] = temp#np.linalg.norm(temp_vel[counter,:])
		#counter += 1
		real_index += 1
	iterations +=1

	htemp, jnk = np.histogram(speed0, bins)
	#real_hist+=htemp
	real_hist += htemp
	#print(real_hist)

	# Plots
	plt.figure()
	#plt.plot(v,maxwellian)
	plt.subplot(1,1,1)
	plt.title(r"Velocity distribution Electrons (part of distr)",fontsize = 16)
	plt.plot(np.linspace(min_Vel,max_Vel,numbins-1),real_hist)
	plt.xlabel(r"Normalized Speed T = %.2f $\displaystyle [\omega_{pi}]$"%(n0*dt*Omega_pi),fontsize = 16)
	#plt.show()
	plt.savefig(path+"speedElectrons.png")
	plt.close()



	min_Vel = 0.
	max_Vel = 0.0 #found later
	real_index = 0 # index of real distribution
	counter = 0 # index of local distribution
	iterations = 1 # to get percentage done

	real_histX = np.zeros(numbins-1, dtype='int32')
	real_histY = np.zeros(numbins-1, dtype='int32')
	real_histZ = np.zeros(numbins-1, dtype='int32')

	temp_velX = np.zeros((mem_step),dtype='float64')
	temp_velY = np.zeros((mem_step),dtype='float64')
	temp_velZ = np.zeros((mem_step),dtype='float64')
	temp_velX = vel0[real_index:mem_step,0]
	temp_velY = vel0[real_index:mem_step,1]
	temp_velZ = vel0[real_index:mem_step,2]
	speed0x = np.zeros(mem_step)
	speed0y = np.zeros(mem_step)
	speed0z = np.zeros(mem_step)
	for counter in range(0,int(mem_step)):
		#Compute particle speed
		speed0x[counter] = temp_velX[counter] #vel0[real_index,0]
		speed0y[counter] = temp_velY[counter] #vel0[real_index,1]
		speed0z[counter] = temp_velZ[counter] #vel0[real_index,2]
		if real_index < mem_step:
			if speed0x[counter] > max_Vel: max_Vel = speed0x[counter] + speed0x[counter]/10.
			if speed0y[counter] > max_Vel: max_Vel = speed0y[counter] + speed0y[counter]/10.
			if speed0z[counter] > max_Vel: max_Vel = speed0z[counter] + speed0z[counter]/10.

			if speed0x[counter] < min_Vel: min_Vel = speed0x[counter] - speed0x[counter]/10.
			if speed0y[counter] < min_Vel: min_Vel = speed0y[counter] - speed0y[counter]/10.
			if speed0z[counter] < min_Vel: min_Vel = speed0z[counter] - speed0z[counter]/10.
		counter += 1
		real_index += 1

	bins = np.linspace(min_Vel, max_Vel, numbins)
	## sucsessive sub distributions
	while real_index < Np-mem_step:
		print("computing velocities %i %% done"%(((float(100*iterations*mem_step))/Np)))
		counter = 0
		#speed0x = np.zeros(mem_step)
		#speed0y = np.zeros(mem_step)
		#speed0z = np.zeros(mem_step)

		temp_velX = vel0[real_index:(iterations+1)*mem_step,0]
		temp_velY = vel0[real_index:(iterations+1)*mem_step,1]
		temp_velZ = vel0[real_index:(iterations+1)*mem_step,2]
		for counter in range(0,int(mem_step)): #could be for loop
			#Compute particle speed
			speed0x[counter] = temp_velX[counter]#vel0[real_index,0]
			speed0y[counter] = temp_velY[counter]#vel0[real_index,1]
			speed0z[counter] = temp_velZ[counter]#vel0[real_index,2]
			real_index += 1
		iterations +=1
		htempX, jnk = np.histogram(speed0x, bins)
		htempY, jnk = np.histogram(speed0y, bins)
		htempZ, jnk = np.histogram(speed0z, bins)
		real_histX += htempX
		real_histY += htempY
		real_histZ += htempZ
		#print(real_hist)

	## LAST asymetric part
	temp_vel = vel0[real_index:Np,:]

	temp_velX = vel0[real_index:Np,0]
	temp_velY = vel0[real_index:Np,1]
	temp_velZ = vel0[real_index:Np,2]

	for counter in range(0,Np-real_index): #could be for loop

		# Compute particle speed
		speed0x[counter] = temp_velX[counter]#vel0[real_index,0]
		speed0y[counter] = temp_velY[counter]#vel0[real_index,1]
		speed0z[counter] = temp_velZ[counter]#vel0[real_index,2]
		#counter += 1
		real_index += 1
	iterations +=1

	htempX, jnk = np.histogram(speed0x, bins)
	htempY, jnk = np.histogram(speed0y, bins)
	htempZ, jnk = np.histogram(speed0z, bins)
	#real_hist+=htemp
	real_histX += htempX
	real_histY += htempY
	real_histZ += htempZ
	#print(real_hist)


	### x,y,z
	plt.figure()
	plt.subplot(3,1,1)
	plt.title(r"Velocity distribution Electrons x,y,z-dimension (part of distr)",fontsize = 16)
	plt.plot(np.linspace(-max_Vel,max_Vel,numbins-1),real_histX)
	plt.subplot(3,1,2)
	plt.plot(np.linspace(-max_Vel,max_Vel,numbins-1),real_histY)
	plt.ylabel("z \hspace{4cm}  y \hspace{4cm}  x",fontsize = 16)
	plt.subplot(3,1,3)
	plt.plot(np.linspace(-max_Vel,max_Vel,numbins-1),real_histZ)
	plt.xlabel(r"Normalized Speed T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)),fontsize = 16)
	#plt.show()
	plt.savefig(path+"speedElectrons0xyz.png")
	plt.close()
	plt.clf()


def plot_vx_vy(step,path,dt,Omega_pi,Omega_e,resolution = 30,mem_step=100000):
	""" Function plots vx vs vy for both species."""

	# Loading file


	print("ploting Vx-Vy")
	file = h5py.File(path + 'pop.pop.h5','r')
	test = file['/pos/specie 1']
	Nt = len(test) # timesteps

	n0 = step + 0.5#1.5

	vel0 = file['/vel/specie 0/n=%.1f' % n0] # Electrons
	print('putting particles in bins specie 0')

	Np = vel0.shape[0]/part_distr	# Number of particles
	Nd = vel0.shape[1]	# Number of dimensions

	mem_step = int(mem_step)
	Np = int(Np)
	#data = np.zeros([resolution,2])
	max_vel = 0.0#0128739455343 #max(vel0[0,:])
	min_vel = 0.0#0128739455343 #min(vel0[0,:])

	print("number of particles = %i"%Np)
	print("ploting %i at a time"%(mem_step))

	real_index = 0 # index of real distribution
	counter = 0 # index of local distribution
	iterations = 1 # to get percentage done

	numbins = resolution
	real_hist = np.zeros((numbins-1,numbins-1), dtype='int32')

	temp_velX = np.zeros((mem_step),dtype='float64')
	temp_velY = np.zeros((mem_step),dtype='float64')
	temp_velX = vel0[real_index:mem_step,0]
	temp_velY = vel0[real_index:mem_step,1]
	speed0x = np.zeros(mem_step)
	speed0y = np.zeros(mem_step)

	### FIRST part

	while counter < mem_step: #could be for loop
		#Compute particle speed
		speed0x[counter] = temp_velX[counter] #vel0[real_index,0]
		speed0y[counter] = temp_velY[counter] #vel0[real_index,1]

		if real_index < mem_step:
			if max_vel < speed0x[counter]: max_vel = speed0x[counter]+speed0x[counter]/10
			if speed0x[counter] < min_vel: min_vel = speed0x[counter]-speed0x[counter]/10
			if max_vel < speed0y[counter]: max_vel = speed0y[counter]+speed0y[counter]/10
			if min_vel > speed0y[counter]: min_vel = speed0y[counter]-speed0y[counter]/10
		counter += 1
		real_index += 1

	## sucsessive part

	while real_index < Np-mem_step:
		print("computing velocities %i %% done"%(((float(100*iterations*mem_step))/Np)))
		counter = 0
		#speed0x = np.zeros(mem_step)
		#speed0y = np.zeros(mem_step)
		#print("allocated")
		temp_velX = vel0[real_index:mem_step*(iterations+1),0]
		temp_velY = vel0[real_index:mem_step*(iterations+1),1]
		#print(temp_velX )
		while counter < mem_step: #could be for loop
			#Compute particle speed
			speed0x[counter] = temp_velX[counter] #vel0[real_index,0]
			speed0y[counter] = temp_velY[counter] #vel0[real_index,1]
			counter += 1
			real_index += 1
		iterations +=1
		#print(min_vel)
		#print(max_vel)
		#print(speed0y)
		bins = np.linspace(min_vel, max_vel, numbins)
		htemp, jnk1,jnk2 = np.histogram2d(speed0x,speed0y, bins=(bins, bins), range=[[min_vel, max_vel], [min_vel, max_vel]])
		#print(bins)

		real_hist += htemp
		#print(len(bins))
		#print(len(real_hist[:,0]))
		#print(len(real_hist[0,:]))
		#print(htemp)

	### LAST asymetric part
	temp_velX = vel0[real_index:Np,0]
	temp_velY = vel0[real_index:Np,1]
	while counter < mem_step: #could be for loop
		#Compute particle speed
		speed0x[counter] = temp_velX[counter] #vel0[real_index,0]
		speed0y[counter] = temp_velY[counter] #vel0[real_index,1]
		counter += 1
		real_index += 1
	iterations +=1
	bins = np.linspace(min_vel, max_vel, numbins)
	htemp, jnk1,jnk2 = np.histogram2d(speed0x,speed0y, bins=(bins, bins), range=[[min_vel, max_vel], [min_vel, max_vel]])
	real_hist += htemp


	### set up figure and plot n0
	real_hist = np.transpose(real_hist)
	bins = np.linspace(min_vel, max_vel, numbins-1)
	fig, ax = plt.subplots(1)
	im = ax.contourf(bins,bins,real_hist, resolution)
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

	plt.savefig(path+"Vx_Vy_Electrons.png")
	plt.clf()

	### logplot
	bins = np.linspace(min_vel, max_vel, numbins-1)
	fig, ax = plt.subplots(1)
	real_hist = np.log10(real_hist)
	real_hist[real_hist<=-inf] =0
	real_hist[real_hist>=inf] =0
	#print(real_hist)
	im = ax.contourf(bins,bins,real_hist, resolution) #locator=ticker.LogLocator()
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

	plt.savefig(path+"Vx_Vy_Electrons_log.png")
	plt.clf()

	### Ions
	test = file['/pos/specie 1']
	Nt = len(test) # timesteps

	vel0 = file['/vel/specie 1/n=%.1f' % n0] # ions
	print('putting particles in bins specie 1')
	Np = vel0.shape[0]/part_distr	# Number of particles
	Nd = vel0.shape[1]	# Number of dimensions
	#data = np.zeros([resolution,2])
	max_vel = 0.0#0128739455343 #max(vel0[0,:])
	min_vel = 0.0#0128739455343 #min(vel0[0,:])

	print("number of particles = %i"%Np)
	print("ploting %i at a time"%(mem_step))

	real_index = 0 # index of real distribution
	counter = 0 # index of local distribution
	iterations = 1 # to get percentage done

	numbins = resolution
	real_hist = np.zeros((numbins-1,numbins-1), dtype='int32')

	temp_velX = np.zeros((mem_step),dtype='float64')
	temp_velY = np.zeros((mem_step),dtype='float64')
	temp_velX = vel0[real_index:mem_step,0]
	temp_velY = vel0[real_index:mem_step,1]
	speed0x = np.zeros(mem_step)
	speed0y = np.zeros(mem_step)

	### FIRST part

	while counter < mem_step: #could be for loop
		#Compute particle speed
		speed0x[counter] = temp_velX[counter] #vel0[real_index,0]
		speed0y[counter] = temp_velY[counter] #vel0[real_index,1]

		if real_index < mem_step:
			if max_vel < speed0x[counter]: max_vel = speed0x[counter]+speed0x[counter]/10
			if speed0x[counter] < min_vel: min_vel = speed0x[counter]-speed0x[counter]/10
			if max_vel < speed0y[counter]: max_vel = speed0y[counter]+speed0y[counter]/10
			if min_vel > speed0y[counter]: min_vel = speed0y[counter]-speed0y[counter]/10
		counter += 1
		real_index += 1

	## sucsessive part

	while real_index < Np-mem_step:
		print("computing velocities %i %% done"%(((float(100*iterations*mem_step))/Np)))
		counter = 0
		#speed0x = np.zeros(mem_step)
		#speed0y = np.zeros(mem_step)
		#print("allocated")
		temp_velX = vel0[real_index:mem_step*(iterations+1),0]
		temp_velY = vel0[real_index:mem_step*(iterations+1),1]
		#print(temp_velX )
		while counter < mem_step: #could be for loop
			#Compute particle speed
			speed0x[counter] = temp_velX[counter] #vel0[real_index,0]
			speed0y[counter] = temp_velY[counter] #vel0[real_index,1]
			counter += 1
			real_index += 1
		iterations +=1
		#print(min_vel)
		#print(max_vel)
		#print(speed0y)
		bins = np.linspace(min_vel, max_vel, numbins)
		htemp, jnk1,jnk2 = np.histogram2d(speed0x,speed0y, bins=(bins, bins), range=[[min_vel, max_vel], [min_vel, max_vel]])
		#print(bins)

		real_hist += htemp
		#print(len(bins))
		#print(len(real_hist[:,0]))
		#print(len(real_hist[0,:]))
		#print(htemp)

	### LAST asymetric part
	temp_velX = vel0[real_index:Np,0]
	temp_velY = vel0[real_index:Np,1]
	while counter < mem_step: #could be for loop
		#Compute particle speed
		speed0x[counter] = temp_velX[counter] #vel0[real_index,0]
		speed0y[counter] = temp_velY[counter] #vel0[real_index,1]
		counter += 1
		real_index += 1
	iterations +=1
	bins = np.linspace(min_vel, max_vel, numbins)
	htemp, jnk1,jnk2 = np.histogram2d(speed0x,speed0y, bins=(bins, bins), range=[[min_vel, max_vel], [min_vel, max_vel]])
	real_hist += htemp

	### set up figure and plot n0
	real_hist = np.transpose(real_hist)
	bins = np.linspace(min_vel, max_vel, numbins-1)
	plt.figure()
	#alternative 1
	fig, ax = plt.subplots(1)
	im = ax.contourf(bins,bins,real_hist, resolution)
	plt.xlabel("Vx")
	plt.ylabel("Vy")
	plt.title(r"Ions  T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)))
	fig.subplots_adjust(bottom = 0.25)
	cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
	fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")
	#plt.show()
	plt.savefig(path+"Vx_Vy_Ions.png")
	plt.clf()

	### logplot
	plt.figure()
	#alternative 1
	fig, ax = plt.subplots(1)
	real_hist = np.log10(real_hist)
	real_hist[real_hist<=-inf] =0
	real_hist[real_hist>=inf] =0
	#print(real_hist)
	im = ax.contourf(bins,bins,real_hist, resolution)#locator=ticker.LogLocator()
	plt.xlabel("Vx")
	plt.ylabel("Vy")
	plt.title(r"Ions  T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)))
	fig.subplots_adjust(bottom = 0.25)
	cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
	fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")
	#plt.show()
	plt.savefig(path+"Vx_Vy_Ions_log.png")
	plt.clf()


def plot_speed_density(dx,dt,nSteps,path,Omega_pi,step,resolution=30,mem_step=1000):

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

	n0 = step + 0.5#1.5
	#print('%.1f, %.1f, %.1f, %.1f'% (n0, n1, n2, n3))
	#pos = file['/pos/specie 1/n=%.1f' % n0]
	#specie 0
	vel0 = file['/vel/specie 0/n=%.1f' % n0] # Electrons
	pos0 = file['/pos/specie 0/n=%.1f' %(n0-0.5)] # Electrons


	Np = vel0.shape[0]/part_distr	# Number of particles
	Nd = vel0.shape[1]	# Number of dimensions
	max_vel = 0.0#0128739455343 #max(vel0[0,:])
	min_vel = -0.0#0128739455343 #min(vel0[0,:])
	max_pos = nSteps#max(pos0[:,0])#0.0
	min_pos = 0.0#min(pos0[:,0])#-0.0

	mem_step = int(mem_step)
	Np = int(Np)


	print("number of particles = %i"%Np)
	print("ploting %i at a time"%(mem_step))

	real_index = 0 # index of real distribution
	counter = 0 # index of local distribution
	iterations = 1 # to get percentage done

	numbins = resolution
	real_hist = np.zeros((numbins-1,numbins-1), dtype='int32')

	#alloc
	temp_velX = np.zeros((mem_step),dtype='float64')
	temp_posX = np.zeros((mem_step),dtype='float64')
	speedX = np.zeros(mem_step)
	posX = np.zeros(mem_step)
	#read part
	temp_velX = vel0[real_index:mem_step,0]
	temp_posX = pos0[real_index:mem_step,1]


	### FIRST part

	while counter < mem_step: #could be for loop
		#Compute particle speed
		speedX[counter] = temp_velX[counter]
		posX[counter] = temp_posX[counter]
		if real_index < mem_step:
			#if max_pos < posX[counter]: max_pos = posX[counter]+posX[counter]/10.
			#if min_pos > posX[counter]: min_pos = posX[counter]-posX[counter]/10.

			if max_vel < speedX[counter]: max_vel = speedX[counter]+speedX[counter]/10.
			if min_vel > speedX[counter]: min_vel = speedX[counter]-speedX[counter]/10.
		counter += 1
		real_index += 1

	### Successive parts

	while real_index < Np-mem_step:
		print("computing velocities %i %% done"%(((float(100*iterations*mem_step))/Np)))
		counter = 0
		#speedX = np.zeros(mem_step)
		#posX = np.zeros(mem_step)
		temp_velX = vel0[real_index:mem_step*(iterations+1),0]
		temp_posX = pos0[real_index:mem_step*(iterations+1),1]
		while counter < mem_step: #could be for loop
			#Compute particle speed
			speedX[counter] = temp_velX[counter]
			posX[counter] = temp_posX[counter]
			counter += 1
			real_index += 1
		iterations +=1
		#print(min_vel)
		#print(max_vel)
		#print(speed0y)
		binsx = np.linspace(min_pos, max_pos, numbins)
		binsy = np.linspace(min_vel, max_vel, numbins)
		htemp, jnk1,jnk2 = np.histogram2d(posX,speedX, bins=(binsx, binsy), range=[[min_pos, max_pos], [min_vel, max_vel]])
		real_hist += htemp
		#print(real_hist)

	### LAST asymetric part
	temp_velX = vel0[real_index:Np,0]
	temp_posX = pos0[real_index:Np,1]
	while counter < mem_step: #could be for loop
		#Compute particle speed
		speedX[counter] = temp_velX[counter]
		posX[counter] = temp_posX[counter]
		counter += 1
		real_index += 1
	iterations +=1

	binsx = np.linspace(min_pos, max_pos, numbins)
	binsy = np.linspace(min_vel, max_vel, numbins)
	htemp, jnk1,jnk2 = np.histogram2d(posX,speedX, bins=(binsx, binsy), range=[[min_pos, max_pos], [min_vel, max_vel]])
	real_hist += htemp
	#print(real_hist)



	### set up figure and plot n0
	binsx = np.linspace(min_pos, max_pos, numbins-1)
	binsy = np.linspace(min_vel, max_vel, numbins-1)

	# ### set up figure and plot n0
	real_hist = np.transpose(real_hist)
	fig, ax = plt.subplots(1)
	im = ax.contourf(binsx,binsy,real_hist, resolution)
	plt.xlabel("x")
	plt.ylabel("Vx")
	plt.title(r"Electrons  T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)))
	fig.subplots_adjust(bottom = 0.25)
	cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
	fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")

	#plt.show()
	plt.rc('text', usetex = True) #latex
	plt.rc('font', family='serif')

	plt.savefig(path+"x_Vx_Electrons.png")
	plt.clf()

	#### LOG plot
	fig, ax = plt.subplots(1)
	real_hist = np.log10(real_hist)
	real_hist[real_hist<=-inf] =0
	real_hist[real_hist>=inf] =0
	im = ax.contourf(binsx,binsy,real_hist, resolution)
	plt.xlabel("x")
	plt.ylabel("Vx")
	plt.title(r"Electrons  T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)))
	fig.subplots_adjust(bottom = 0.25)
	cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
	fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")

	#plt.show()
	plt.rc('text', usetex = True) #latex
	plt.rc('font', family='serif')

	plt.savefig(path+"x_Vx_Electrons_LOG.png")
	plt.clf()


	#specie 1
	vel0 = file['/vel/specie 1/n=%.1f' % n0] # Ions
	pos0 = file['/pos/specie 1/n=%.1f' %(n0-0.5)] # Ions
	print('putting particles in bins specie 0')

	Np = vel0.shape[0]/part_distr	# Number of particles
	Nd = vel0.shape[1]	# Number of dimensions
	max_vel = 0.0#0128739455343 #max(vel0[0,:])
	min_vel = -0.0#0128739455343 #min(vel0[0,:])
	max_pos = nSteps#0.0
	min_pos = -0.0


	print("number of particles = %i"%Np)
	print("ploting %i at a time"%(mem_step))

	real_index = 0 # index of real distribution
	counter = 0 # index of local distribution
	iterations = 1 # to get percentage done

	numbins = resolution
	real_hist = np.zeros((numbins-1,numbins-1), dtype='int32')



	#alloc
	temp_velX = np.zeros((mem_step),dtype='float64')
	temp_posX = np.zeros((mem_step),dtype='float64')
	speedX = np.zeros(mem_step)
	posX = np.zeros(mem_step)
	#read part
	temp_velX = vel0[real_index:mem_step,0]
	temp_posX = pos0[real_index:mem_step,1]


	### FIRST part

	while counter < mem_step: #could be for loop
		#Compute particle speed
		speedX[counter] = temp_velX[counter]
		posX[counter] = temp_posX[counter]
		if real_index < mem_step:
			#if max_pos < posX[counter]: max_pos = posX[counter]+posX[counter]/10.
			#if min_pos > posX[counter]: min_pos = posX[counter]-posX[counter]/10.

			if max_vel < speedX[counter]: max_vel = speedX[counter]+speedX[counter]/10.
			if min_vel > speedX[counter]: min_vel = speedX[counter]-speedX[counter]/10.
		counter += 1
		real_index += 1

	### Successive parts

	while real_index < Np-mem_step:
		print("computing velocities %i %% done"%(((float(100*iterations*mem_step))/Np)))
		counter = 0
		#speedX = np.zeros(mem_step)
		#posX = np.zeros(mem_step)
		temp_velX = vel0[real_index:mem_step*(iterations+1),0]
		temp_posX = pos0[real_index:mem_step*(iterations+1),1]
		while counter < mem_step: #could be for loop
			#Compute particle speed
			speedX[counter] = temp_velX[counter]
			posX[counter] = temp_posX[counter]
			counter += 1
			real_index += 1
		iterations +=1
		#print(min_vel)
		#print(max_vel)
		#print(speed0y)
		binsx = np.linspace(min_pos, max_pos, numbins)
		binsy = np.linspace(min_vel, max_vel, numbins)
		htemp, jnk1,jnk2 = np.histogram2d(posX,speedX, bins=(binsx, binsy), range=[[min_pos, max_pos], [min_vel, max_vel]])
		real_hist += htemp
		#print(real_hist)

	### LAST asymetric part
	temp_velX = vel0[real_index:Np,0]
	temp_posX = pos0[real_index:Np,1]
	while counter < mem_step: #could be for loop
		#Compute particle speed
		speedX[counter] = temp_velX[counter]
		posX[counter] = temp_posX[counter]
		counter += 1
		real_index += 1
	iterations +=1

	binsx = np.linspace(min_pos, max_pos, numbins)
	binsy = np.linspace(min_vel, max_vel, numbins)
	htemp, jnk1,jnk2 = np.histogram2d(posX,speedX, bins=(binsx, binsy), range=[[min_pos, max_pos], [min_vel, max_vel]])
	real_hist += htemp
	#print(real_hist)



	### set up figure and plot n0
	real_hist = np.transpose(real_hist)
	binsx = np.linspace(min_pos, max_pos, numbins-1)
	binsy = np.linspace(min_vel, max_vel, numbins-1)

	fig, ax = plt.subplots(1)
	im = ax.contourf(binsx,binsy,real_hist, resolution)
	plt.xlabel("x")
	plt.ylabel("Vx")
	plt.title(r"Ions  T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)))
	fig.subplots_adjust(bottom = 0.25)
	cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
	fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")

	#plt.show()

	plt.savefig(path+"x_Vx_Ions.png")
	plt.clf()

	######## LOG plot

	### set up figure and plot n0
	binsx = np.linspace(min_pos, max_pos, numbins-1)
	binsy = np.linspace(min_vel, max_vel, numbins-1)

	fig, ax = plt.subplots(1)
	real_hist = np.log10(real_hist)
	real_hist[real_hist<=-inf] =0
	real_hist[real_hist>=inf] =0
	im = ax.contourf(binsx,binsy,real_hist, resolution)
	plt.xlabel("x")
	plt.ylabel("Vx")
	plt.title(r"Ions  T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)))
	fig.subplots_adjust(bottom = 0.25)
	cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
	fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")

	#plt.show()

	plt.savefig(path+"x_Vx_Ions_LOG.png")
	plt.clf()


	print(real_index)



def plot_kx_ky_grid(dx,dt,path,Omega_pi,nSteps,timeStep = 0,resolution =64,grid = "E"):

	print("using timestep %i"%(timeStep))
	if grid == "E":
		print("using E grid")
		h5file = h5py.File(path +'E.grid.h5','r')
		E = np.asarray(h5file['/n=%.1f'%(timeStep)])
		#print(E.shape)

		#print(E.shape)
		Nx = len(E[:,0,0])
		Ny = len(E[0,:,0])
		Nz = len(E[0,:,0])

		array = np.zeros((Nx,Ny,Nz))
		##print("%f,%f,%f"%(i,j,k))
		for i in range(Nx): # x
			for j in range(Ny): #y
				for k in range(Nz):
					array[i,j,k] = np.sqrt(E[i,j,k,0]*E[i,j,k,0]+E[i,j,k,1]*E[i,j,k,1]+E[i,j,k,2]*E[i,j,k,2])

		array = np.transpose(array, (2,1,0))
		array = np.squeeze(array)

	if grid == "rho":
		h5file = h5py.File(path +'rho.grid.h5','r')
		rho = np.asarray(h5file['/n=%.1f'%(timeStep)])
		array = np.transpose(rho, (3,2,1,0))
		array = np.squeeze(array)
		print("using rho grid")

	if grid == "phi":
		h5file = h5py.File(path +'phi.grid.h5','r')
		phi = np.asarray(h5file['/n=%.1f'%(timeStep)])
		array = np.transpose(phi, (3,2,1,0))
		array = np.squeeze(array)
		print("using rho phi")
	#FFT_array = np.fft.fft2(array[:,:,Nz/2])
	#FFT_array = np.fft.fftshift(FFT_array)
	###plt.plot(np.abs(FFT_array))
	#plt.contourf((np.abs(FFT_array)))
	###plt.hist2d(FFT_array,FFT_array)
	#plt.show()
	#print((abs(np.real(FFT_array))))


	nx, ny     = (nSteps, nSteps)
	xmax, ymax = nx*dx, ny*dx
	x          = np.linspace(-xmax, xmax, nx)
	y          = np.linspace(-ymax, ymax, ny)
	#dx         = x[1] - x[0]
	dy         = dx#y[1] - y[0]
	X, Y       = np.meshgrid(x, y)
	Z          = array[:,:,0]#(X*X+Y*Y)<1**2          # circular hole

	#Z          = np.exp(-(X*X + Y*Y)/1**2)   # Gauss
	#Z           = (np.abs(X)<0.5 ) * (np.abs(Y)<2 )  # rectangle
	#Z          = (((X+1)**2+Y**2)<0.25**2) + (((X-1)**2+Y**2)<0.25**2)   # two circular holes

	ZFT    = np.fft.fft2(Z)        # compute 2D-FFT
	ZFT    = np.fft.fftshift(ZFT)  # Shift the zero-frequency component to the center of the spectrum.
	kx     = (-nx/2 + np.arange(0,nx))*2*np.pi/(2*xmax)
	ky     = (-ny/2 + np.arange(0,ny))*2*np.pi/(2*ymax)
	KX, KY = np.meshgrid(kx, ky)
	"""
	nx, ny     = (500, 500)
	xmax, ymax = 20, 20
	x          = np.linspace(-xmax, xmax, nx)
	y          = np.linspace(-ymax, ymax, ny)
	dx         = x[1] - x[0]
	dy         = y[1] - y[0]
	X, Y       = np.meshgrid(x, y)
	Z          = (X*X+Y*Y)<1**2          # circular hole
	#Z          = np.exp(-(X*X + Y*Y)/1**2)   # Gauss
	#Z           = (np.abs(X)<0.5 ) * (np.abs(Y)<2 )  # rectangle
	#Z          = (((X+1)**2+Y**2)<0.25**2) + (((X-1)**2+Y**2)<0.25**2)   # two circular holes

	ZFT    = np.fft.fft2(Z)        # compute 2D-FFT
	ZFT    = np.fft.fftshift(ZFT)  # Shift the zero-frequency component to the center of the spectrum.
	kx     = (-nx/2 + np.arange(0,nx))*2*np.pi/(2*xmax)
	ky     = (-ny/2 + np.arange(0,ny))*2*np.pi/(2*ymax)
	KX, KY = np.meshgrid(kx, ky)


	"""
	# plot Z and it's Fourier transform ZFT
	fig, (ax0, ax1) = plt.subplots(ncols=2)
	plot1 = ax0.pcolormesh(X,Y,Z)
	plot2 = ax1.pcolormesh(KX,KY,np.abs(np.log10(ZFT)))
	ax0.set_title('real space')
	ax1.set_title('Fourier space')
	ax0.set_xlabel('x')
	ax0.set_ylabel('y')
	ax1.set_xlabel('kx')
	ax1.set_ylabel('ky')
	plt.show()

	"""
	FFTx = np.linspace(0,1./(2*(dx)),N/2)
	X_0 = np.real(np.fft.fft(x[0,:]))
	Y_0 = (np.fft.fft(y[0,:]))
	for i in range(len(t)-1000,len(t)):
		X_0 += np.real(np.fft.fft(x[i,:]))
		Y_0 += np.real(np.fft.fft(y[i,:]))


	def testdata(N):
		data = np.zeros(N)
		k = 1.
		for i in range(N):
			data[i] = np.sin(k*i)
		fftdata = np.fft.fft(data)
		return fftdata, data
	"""







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


part_distr = 1 #part of distribution to plot, for debugging
def post_process_all(starttime,endtime,step,path,dt,dx):
	"""TLWR, kernel "sort of" """

	# DEFAULT VALUES
	nSteps = 128 # number of spatial steps
	q = 1.602e-19 # charge
	B = 0.000015#7.5e-6 # mag. field
	M_i = 5e-26 #mass Ions
	M_e = 4e-29 #9.109e-31#
	n_0 = 1.*10**(9)
	Omega_i = (q*B)/M_i
	Omega_e = (q*B)/M_e
	Omega_pe = np.sqrt((n_0*q*q)/(M_e*eps_0))
	Omega_pi = np.sqrt((n_0*q*q)/(M_i*eps_0))
	Omega_e = Omega_pi #using ion plasma frq
	final_timestep = endtime#30000

	#plot_electic_field_in_time(starttime,endtime,step,path,dt,Omega_i,Omega_e) # high cost

	# grid plots

	#animate_by_name("phi",path,starttime,endtime,step,dt,Omega_i,Omega_e,dx)
	animate_by_name("rho",path,starttime,endtime,step,dt,Omega_i,Omega_e,dx)
	
	#plot_electic_field_in_time(starttime,endtime,step,path,dt,Omega_i,Omega_e) # high cost

	#plot_energy(path)
	#plot_temperature(path,dt,Omega_i,Omega_e)

	import time

	start = time.time()
	#plot_velocity_distribution(final_timestep,path,dt,Omega_pi,Omega_e,mem_step=1e5)
	end = time.time()
	print("plot_velocity_distibution used %f sec"%(end - start))

	start = time.time()
	#plot_vx_vy(final_timestep,path,dt,Omega_pi,Omega_e,resolution = 50,mem_step=1e5)
	end = time.time()
	print("plot_vx_vy used %f sec"%(end - start))

	start = time.time()
	#plot_speed_density(dx,dt,nSteps,path,Omega_pi,final_timestep,resolution = 50,mem_step=1e5) #step = 1000
	end = time.time()
	print("plot_speed_density used %f sec"%(end - start))

	#plot_kx_ky_grid(dx,dt,path,Omega_pi,nSteps,timeStep = final_timestep, grid = "E")


post_process_all(0,50000,100,path="../../data/",dt = 3e-6,dx=0.08)
