#Some rudimentary plotting of distributions

import h5py
import numpy as np
import pylab as plt
# from mayavi import mlab
#from numba import jit
import sys, os, inspect

#constants
k_b = 1.38064852*10**(-23)

#@jit
def findAverageESquare(input_File,starttime,endtime):
	"""find what timesteps E exists, and return
	array with E^2 values"""
	
	
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
						temp = (E[i,j,k,0]*E[i,j,k,0]+E[i,j,k,1]*E[i,j,k,1]+E[i,j,k,2]*E[i,j,k,2])
						val += temp/(size*size)
						
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


def plot_electic_field_in_time(starttime,endtime,step,path):
	print("working the Electic field")
	print("ploting %i timesteps, this may take some time"%((endtime-starttime)/step))
	h5file = h5py.File(path +'E.grid.h5','r')

	#compile
	#a,b,c = findAverageESquare(h5file,0,0)

	#run
	avg_E_square,max_E_square,time = findAverageESquare(h5file,starttime,endtime)
	plt.clf()
	plt.plot(time,avg_E_square,label = "averaged E^2")
	#plt.plot(time,max_E_square,label = "max E^2")
	plt.legend(loc='lower right')
	plt.xlabel("Time")
	plt.ylabel("V^2/m^2")
	plt.savefig((path +"ElecticFieldAvg.png")) #save to file
	plt.plot(time,max_E_square,label = "max E^2")
	plt.savefig((path +"ElecticFieldMax.png"))
	plt.gcf().clear()
	#plt.show()

def animate(title,path,subdir,h5,startindex,stopindex,step):
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
		im = ax.contourf(X,Y,data[:,:,4], 100)
	
		fig.subplots_adjust(bottom = 0.25)
		cbar_rho = fig.add_axes([0.10, 0.05, 0.8, 0.10])
		fig.colorbar(im, cax=cbar_rho, orientation = "horizontal")
		plt.title(title+" perpendicular to B_0, t=%i"%i);
		
		plt.savefig((path +subdir +"dt%i.png")%count) #save to file
		count +=1
		#plt.show()
				
def animate_by_name(name,path,starttime,endtime,step):
	"""given "name" of file uses function animate() 
	to open file "name" and write to subdir "name" """
	print("working the "+name+" field")
	h5 = h5py.File('../../data/'+name+'.grid.h5','r')
	if not os.path.exists(path + name+"/"):
		os.mkdir(path + name+"/")
	animate(name,path,name+"/",h5,starttime,endtime,step)
		
def plot_energy(path):
	fig = plt.figure()
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
	plt.plot(pot,label='specie0')
	#plt.plot(pot1,label='specie1')
	plt.xlabel("t/dt")
	plt.ylabel("potential Energy")
	plt.legend(loc='lower right')
	plt.savefig((path +"Energy_potential.png")) #save to file
	plt.gcf().clear()
	#plt.show()

def plot_temperature(path):
	plt.clf()
	hist = h5py.File('../../data/history.xy.h5','r')
	pot = hist['/energy/potential/total']
	kin = hist['/energy/kinetic/total']


	#kin = kin[:,1];		# Extract y-axis
	#pot = pot[:,1];	# Extract y-axis and invert
	#tot = pot+kin;		# Collect total energy

	#specie0
	kin0 = hist['/energy/kinetic/specie 0']/k_b
	kin0 = kin0[:,1];		# Extract y-axis (x is time)

	#specie1
	kin1 = hist['/energy/kinetic/specie 1']/k_b
	kin1 = kin1[:,1];		# Extract y-axis (x is time)


	#plt.plot(pot,label='potential')
	plt.plot(kin0,label='specie0')
	plt.plot(kin1,label='specie1')
	plt.xlabel("t/dt")
	plt.ylabel("Total  Kinetic Energy")
	plt.legend(loc='lower right')
	plt.savefig((path +"Energy.png")) #save to file
	#plt.show()


def post_process_all(starttime,endtime,step,path):
	"""TLWR, kernel "sort of" """

	#plot_electic_field_in_time(starttime,endtime,step,path)

	# grid plots
	
	#animate_by_name("rho",path,starttime,endtime,step)
	#animate_by_name("phi",path,starttime,endtime,step)
	
	plot_energy(path)
	
	#rho





post_process_all(0,40000,100,path="../../data/")


























