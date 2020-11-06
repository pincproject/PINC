
#!/usr/bin/python
import h5py
import numpy as np
import pylab as plt
import sys, os, inspect


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


def find_omega_k(starttime,endtime,step,path,dt,dx):
	"""TLWR, kernel "sort of" """

	# DEFAULT VALUES
	nSteps = 16 # number of spatial steps
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


	plot_kx_ky(dx,dt,path,Omega_pi,nSteps)


find_omega_k(0,10000,100,path="../../data/",dt = 8e-7,dx=0.2)

