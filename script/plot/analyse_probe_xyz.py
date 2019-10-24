
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
import sys, os, inspect
from matplotlib import ticker, cm
import matplotlib.colors as colors

#constants
k_b = 1.38064852*10**(-23)
eps_0 = 8.854*10**(-12)
q = 1.602e-19 # charge
B = 7.5e-6 # mag. field
M_i = 5e-26 #mass Ions
M_e = 9.109e-31#4e-29 
n_0 = 1.*10**(9)

# DEFAULT (GLOBAL) VALUES
nSteps = 128 # number of spatial steps
Omega_i = (q*B)/M_i
Omega_e = (q*B)/M_e
Omega_pe = np.sqrt((n_0*q*q)/(M_e*eps_0))
Omega_pi = np.sqrt((n_0*q*q)/(M_i*eps_0))


def constuct_test_array(starttime,endtime,dx,dt):
	
	testfunct = lambda A,B, omega, t, k, x: A*np.sin(2*np.pi*omega*t - 2*np.pi*k*x) + B*np.sin(2.*2*np.pi*omega*t - 2.*2*np.pi*k*x) + A*np.sin(3.*2*np.pi*omega*t - 3.*2*np.pi*k*x)
	x = np.linspace(0,nSteps+1,nSteps)*dx
	omega = 1000
	k = 1/3.
	Nt = endtime-starttime
	N = nSteps

	A = 1. 
	B = 2.
	testarray = np.zeros((Nt,N))
	#print(testarray)
	time = []	
	for t in range(0,endtime-starttime):
		testarray[t,:] = testfunct(A,B, omega, t*dt, k, x)
		time.append(t+starttime)
	
	plt.plot(x,testarray[0,:])
	plt.title("k = (2*pi*) %f, %f, %f, lambda= %f, %f, %f"%(k, 2*k, 3*k,1./k, 1./(2*k), 1./(3*k)))
	plt.xlabel("x(m)")
	plt.ylabel("constructed data (rho)")
	plt.show()
	plt.plot(time,testarray[:,0])
	plt.title("omega = %f, %f, %f"%(omega, 2.*omega, 3.*omega))
	plt.xlabel("time (s)")
	plt.ylabel("constructed data (rho)")
	plt.show()
	plt.contourf(testarray)
	plt.xlabel("time (s)")
	plt.ylabel("k(2*pi/m)")
	plt.show()
	
	return (testarray)

def construct_omega_k_testarr(starttime,endtime,dx,dt,array):

	#probe = h5py.File(path+'probe.xyz.h5','r')
	#data = probe['/X']

	
	Nt = endtime-starttime
	N = nSteps
	
	
	FFTx = np.linspace(0,1./(2*(dx)),N/2)
	FFTtime = np.linspace(0,1./(2*(dt)),Nt/2)
	

	testdataX = np.fft.fft(array[0,:])[:N//2]
	testdataT = np.fft.fft(array[:,0])[:Nt//2]
	
	plt.plot(FFTx,abs(np.real(testdataX)))
	plt.xlabel("k = 1/x(m)")
	plt.ylabel("FFT(x) Magnitude")
	plt.show()

	plt.plot(FFTtime,abs(np.real(testdataT)))
	plt.xlabel("omega = 1/t(s)")
	plt.ylabel("FFT(time) Magnitude")
	plt.show()
	

	#FFTx = np.linspace(0,1./(2*(dx)),N)
	#FFTtime = np.linspace(0,1./(2*(dt)),Nt)

	#FFTarray = np.fft.fft2(np.real(array))
	FFTarray = construct_omega_k(array)
	FFTarray = np.real(FFTarray)


	fig, ax = plt.subplots()
	print("FFTarray t = %f, x = %f"%(len(FFTarray[:,0]),len(FFTarray[0,:])))
	print("t = %f, x = %f"%(len(FFTtime),len(FFTx)))
	cs = ax.contourf(FFTx,FFTtime,abs(FFTarray[:Nt//2,:N//2])) #,locator=ticker.LogLocator()
	cbar = fig.colorbar(cs)
	plt.show()	
	"""
	testdata = np.fft.fft(data[(time_start_index),:])[:N//2]
	#print(np.real(testdata))
	w_k_arr = np.zeros((Nt/2,N/2))
	used_k = 0 #bogus
	ftdataX = np.zeros(N/2,dtype=complex)
	
	for filter_k in range(N/2):
		omegadata = np.zeros(Nt/2)
		for i in range(0, Nt/2):
			#print("filter_k = %i, i = %i"%(filter_k,(time_start_index+i)))
			#print(Nt+time_start_index-11)
			k = np.fft.fft(data[i,:])[:N//2] #transform x at timestep i
			fftdataX, used_k = bandpass(data[(time_start_index+i),:],1,filter_k,dx,k) #returnes filtered k data
			#print("used_k = %i"%used_k)
			ifftdataX = np.fft.ifft(fftdataX)
			ifftdataX = ifftdataX#/max(ifftdataX)
			omegadata[i] = np.real(ifftdataX[1]) #pick one point, could do every
		w_k_arr[:,filter_k] = np.real(np.fft.fft(omegadata)[:Nt//2])

		 
		#plt.plot(k,2.0/N * abs(fftdataX[:N//2]))
		#plt.plot(x,ifftdataX)



	plt.contourf(w_k_arr)
	plt.show()
	"""

def construct_omega_k_data(starttime,endtime,dx,dt,path,Cs,Csi):



	probe = h5py.File(path+'probe.xyz.h5','r')
	data = probe['/X']
	
	Nt = endtime-starttime
	N = nSteps
	array = data[starttime:endtime,:]
	
	FFTx = np.linspace(0,(1.)/(2*(dx)),N)
	FFTtime = np.linspace(0,(1.)/(2*(dt)),Nt)
	
	

	testdataX = np.fft.fft(array[0,:])[:N//2]
	testdataT = np.fft.fft(array[:,0])[:Nt//2]
	
	##plt.plot(FFTx,abs(np.real(testdataX)))
	###plt.xlabel("k = 1/x(m)")
	#plt.ylabel("FFT(x) Magnitude")
	#plt.show()

	#plt.plot(FFTtime,abs(np.real(testdataT)))
	#plt.xlabel("omega = 1/t(s)")
	#plt.ylabel("FFT(time) Magnitude")
	#plt.show()
	

	#FFTx = np.linspace(0,1./(2*(dx)),N)
	#FFTtime = np.linspace(0,1./(2*(dt)),Nt)

	#FFTarray = np.fft.fft2(array)
	FFTarray = construct_omega_k(array)
	FFTarray = np.real(FFTarray)

	"""
	fig, ax = plt.subplots()
	print("FFTarray t = %f, x = %f"%(len(FFTarray[:,0]),len(FFTarray[0,:])))
	print("t = %f, x = %f"%(len(FFTtime),len(FFTx)))
	cs = ax.pcolormesh(FFTx,FFTtime,abs(FFTarray[:Nt//2,:N//2]),norm=colors.Normalize(vmin=10000, vmax=100000)) # ,locator=ticker.LogLocator()
	cbar = fig.colorbar(cs)
	
	plt.show()
	"""
	fig, ax = plt.subplots()
	analytic_omega = FFTx[:N/16]*Cs#*(1+Psi)
	analytic_omega_iso = FFTx[:N/16]*Csi#*(1+Psi)
	print("FFTarray t = %f, x = %f"%(len(FFTarray[:,0]),len(FFTarray[0,:])))
	print("t = %f, x = %f"%(len(FFTtime),len(FFTx)))
	cs = ax.contourf(FFTx[:N//2],FFTtime[:Nt//2],(abs(FFTarray[:Nt//2,:N//2])),1000) # ,locator=ticker.LogLocator()
	plt.plot(FFTx[:N/16],analytic_omega,"black")
	plt.plot(FFTx[:N/16],analytic_omega_iso,"white")
	cbar = fig.colorbar(cs)
	#plt.rc('text', usetex = True) #latex
	#plt.rc('font', family='serif')
	#plt.title(r"Electrons  T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)))
	#plt.xlabel(r"k $\displaystyle [2\pi / m]$", fontsize = 16)
	#plt.ylabel(r"$\displaystyle  \Omega [2\pi /t] $", fontsize = 16)
	plt.xlabel("k [2pi/m]", fontsize = 16)
	plt.ylabel("Omega [2pi/t]", fontsize = 16)
	plt.show()

	fig, ax = plt.subplots()
	analytic_omega = FFTx[:N/16]*Cs#*(1+Psi)
	analytic_omega_iso = FFTx[:N/16]*Csi#*(1+Psi)
	print("FFTarray t = %f, x = %f"%(len(FFTarray[:,0]),len(FFTarray[0,:])))
	print("t = %f, x = %f"%(len(FFTtime),len(FFTx)))
	cs = ax.contourf(FFTx[:N//2],FFTtime[:Nt//2],np.log10(abs(FFTarray[:Nt//2,:N//2])),1000) # ,locator=ticker.LogLocator()
	plt.plot(FFTx[:N/16],analytic_omega,"black")
	plt.plot(FFTx[:N/16],analytic_omega_iso,"white")
	cbar = fig.colorbar(cs)
	#plt.rc('text', usetex = True) #latex
	#plt.rc('font', family='serif')
	#plt.title(r"Electrons  T = %.2f $\displaystyle [\omega_{pi}]$"%((n0*dt*Omega_pi)))
	#plt.xlabel(r"k $\displaystyle [2\pi / m]$", fontsize = 16)
	#plt.ylabel(r"$\displaystyle  \Omega [2\pi /t] $", fontsize = 16)
	plt.xlabel("k [2pi/m]", fontsize = 16)
	plt.ylabel("Omega [2pi/t]", fontsize = 16)
	plt.show()		
	"""
	testdata = np.fft.fft(data[(time_start_index),:])[:N//2]
	#print(np.real(testdata))
	w_k_arr = np.zeros((Nt/2,N/2))
	used_k = 0 #bogus
	ftdataX = np.zeros(N/2,dtype=complex)
	
	for filter_k in range(N/2):
		omegadata = np.zeros(Nt/2)
		for i in range(0, Nt/2):
			#print("filter_k = %i, i = %i"%(filter_k,(time_start_index+i)))
			#print(Nt+time_start_index-11)
			k = np.fft.fft(data[i,:])[:N//2] #transform x at timestep i
			fftdataX, used_k = bandpass(data[(time_start_index+i),:],1,filter_k,dx,k) #returnes filtered k data
			#print("used_k = %i"%used_k)
			ifftdataX = np.fft.ifft(fftdataX)
			ifftdataX = ifftdataX#/max(ifftdataX)
			omegadata[i] = np.real(ifftdataX[1]) #pick one point, could do every
		w_k_arr[:,filter_k] = np.real(np.fft.fft(omegadata)[:Nt//2])

		 
		#plt.plot(k,2.0/N * abs(fftdataX[:N//2]))
		#plt.plot(x,ifftdataX)



	plt.contourf(w_k_arr)
	plt.show()
	"""

def construct_omega_k(array):

	def bandpass(data,eps,center,k,BW = 1):
	
		""" Center value is the value of k that is kept, k is the 
		fourier transformed array in one dim
		"""
		used_k = 0
		fftdata = np.fft.fft(data)
		newdata = np.zeros(len(fftdata),dtype=complex)
		index = 0
		diff = 1e64
		largest_data = 0
		#print("centered index = %f"%center)
		for i in range(len(k)):
			if abs(i*eps - center) < diff :
				#print("lekflf = %i"%(abs(i*eps - center)))
				#print(i)
				index = i
				diff = abs(eps*i-center) 
		#print("found index %i, k = %f"%(index,k[index]))
		used_k = k[index]
		newdata[index] = fftdata[index]
		#for i in range(center-BW,center+BW+1):
		#print("index = %i, i = %i"%(index,i))
		try:
			#print("FOUND")
			newdata[index-1] = fftdata[index-1]
			newdata[index+1] = fftdata[index+1]
		except:
			None
		return newdata, used_k


	N = len(array[0,:]) #space
	Nt = len(array[:,0]) #time

	w_k_arr = np.zeros((Nt,N))
	used_k = 0 #bogus
	ftdataX = np.zeros(N,dtype=complex)
	
	for filter_k in range(N): # every k
		print("processing %i of %i FFTs"%(filter_k,N))
		for j in range(N/8): # length of box (x dim)
			
			omegadata = np.zeros(Nt)
			for i in range(0, Nt): # construct filtered array in time
				#print("filter_k = %i, i = %i"%(filter_k,(time_start_index+i)))
				#print(Nt+time_start_index-11)
				k = np.fft.fft(array[i,:]) #transform x at timestep i
				fftdataX, used_k = bandpass(array[i,:],1,filter_k,k,1) #returns filtered k data
				#print("used_k = %i"%used_k)
				ifftdataX = np.fft.ifft(fftdataX)
				omegadata[i] = np.real(ifftdataX[j]) #pick one point, could do every # j is every
			w_k_arr[:,filter_k] += np.real(np.fft.fft(omegadata))
	return w_k_arr






def construct_time_x_data(starttime,endtime,dx,dt,path,Cs,Csi,Omega_i):


	import scipy
	from scipy import ndimage

	probe = h5py.File(path+'probe.xyz.h5','r')
	data = probe['/X']
	h5file = h5py.File(path +'phi.grid.h5','r')
	denorm = h5file.attrs.__getitem__("Quantity denormalization factor")
	
	
	#N = nSteps
	array = np.array(denorm*data[starttime:endtime,:])[:,::-1]
	#array = np.flip(array, 0)

	Nt = len(array[:,0])
	nSteps = len(array[0,:])
	print(Nt)

	nx, nt     = (nSteps, Nt)
	xmax, ymax = nx*dx, nt*dt
	x          = np.linspace(0, xmax, nx)
	y          = np.linspace(0, ymax, nt)
	#dx         = x[1] - x[0]
	dy         = dt#y[1] - y[0]
	X, Y       = np.meshgrid(x, y)
	Z          = array[:,:]#(X*X+Y*Y)<1**2          # circular hole

	#Z          = np.exp(-(X*X + Y*Y)/1**2)   # Gauss
	#Z           = (np.abs(X)<0.5 ) * (np.abs(Y)<2 )  # rectangle
	#Z          = (((X+1)**2+Y**2)<0.25**2) + (((X-1)**2+Y**2)<0.25**2)   # two circular holes

	ZFT    = np.fft.fft2(Z)        # compute 2D-FFT
	#ZFT = np.transpose(ZFT)
	#print(ZFT)
	ZFT    = np.fft.fftshift(ZFT)  # Shift the zero-frequency component to the center of the spectrum.
	kx     = np.linspace(-(1.)/((2*dx)),(1.)/((2*dx))+dx,nx+1)#(np.arange(-nx/2,nx/2))*((1.)/(2*(dx)))
	ky     = np.linspace(-(1.)/((2*dt)),(1.)/((2*dt))+dt,nt+1)#(np.arange(-nt/2,nt/2))*((1.)/(2*(dt)))
	#print(np.arange(0,nx/2))
	KX, KY = np.meshgrid( kx, ky)

	analytic_omega = kx*Cs#*(1+Psi)
	analytic_omega_iso = kx*Csi#*(1+Psi)

	# plot Z and it's Fourier transform ZFT

	fig = plt.figure(1)
	fig.clf()
	#plotsingle = ax.pcolormesh(np.transpose(Z))
	ax1 = fig.add_subplot(1,1,1)
	plot1 = ax1.pcolormesh(Omega_i*y,x,np.rot90(Z))	

	plt.rc('text', usetex = True)
	plt.rc('font', family='serif')
	
	ax1.set_title(r"Keogram of $\displaystyle \phi$")	
	ax1.set_ylabel('x-space [m]')	
	ax1.set_xlabel(r"Time \ $\displaystyle[\Omega_{i}]$")
	cbar = plt.colorbar(plot1, ax = ax1)
	cbar.set_label(r'electric potential [V]', rotation=270)
	#cbar.set_label(r'Log scaled', rotation=270)	
	plt.show()

	fig = plt.figure(2)
	#plot1 = ax0.pcolormesh(X,Y*Omega_i,Z)
	cmap = plt.get_cmap('YlOrRd')
	fig.clf()
	ax1 = fig.add_subplot(1,1,1)
	plot2 = ax1.pcolormesh((KX),(KY),((np.abs(np.log10(ZFT)))),cmap=cmap)
	#plt.gca().invert_xaxis()
	#ax0.set_title('real space')
	ax1.set_title(r'Fourier space')
	ax1.plot(kx,analytic_omega,"blue",linestyle='--',label = "Cs")
	plt.plot(kx,analytic_omega_iso,"blue",label= "C_{si}")
	cbar = plt.colorbar(plot2, ax = ax1)
	cbar.set_label(r'Log scaled magnitude', rotation=270)	
	#ax0.set_xlabel('x')
	#ax0.set_ylabel('t')
	ax1.set_xlabel(r"$\displaystyle k_x \ [2 \pi / m]$")
	ax1.set_ylabel(r"$\displaystyle \omega \ [2 \pi / s]$")
	plt.legend(fontsize=12)
	plt.show()
	#plt.savefig((path +"omega_k_Filtered_zoomed_last.png"))





def find_omega_k(starttime,endtime,path,dt,dx):
	"""TLWR, kernel "sort of" """
	
	## TEST suite
	#testarray = constuct_test_array(starttime,endtime,dx,dt)
	#construct_omega_k_testarr(starttime,endtime,dx,dt,testarray)
	
	q = 1.602e-19 # charge
	M_i = 5e-26 #mass Ions
	M_e = 4e-29 
	gamma_a = 1.
	gamma_iso = 5/3.
	T_i = k_b*300 # K
	T_e = k_b*300 # K
	B0 = 0.00005 # mag. field
	E0 = 0.1
	nu_e = 2800
	nu_i = 1800
	Omega_i = (q*B0)/M_i
	Omega_e = (q*B0)/M_e
	Psi = (nu_e*nu_i)/(Omega_e*Omega_i)
	Csa = np.sqrt((gamma_a*T_i+gamma_a*T_e)/M_i)
	Csi = np.sqrt((gamma_iso*T_i+gamma_iso*T_e)/M_i)
	Vd = E0/B0
	Cs = Vd/(1+Psi)
	print(Cs)

	#construct_omega_k_data(starttime,endtime,dx,dt,path,Cs,Csi)

	construct_time_x_data(starttime,endtime,dx,dt,path,Cs,Csi,Omega_i)


if __name__ == "__main__":
	find_omega_k(12500,25000,path="../../data/",dt = 3e-6,dx=0.04)

