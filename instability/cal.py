import numpy as np
import math
from langmuir import *
import scipy.constants as sc


Q          = 1.6e-19
KB         = 1.38e-23
EPS0       = 8.85418782e-12
eV = 1.6e-19

mAr	       = 1.6724828e-27
mE         = 9.10938356e-31
speedLight = 3e8
tEcK        = 11600
tEhK        = 1160452
tIK        = 1000
pDen       = 1e10
tEheV      = 10
tEceV      = 1.0
tEbeV      = 5e3
tIeV	   = 0.026
vthE       = 3e4
# vthE       = np.sqrt(2*(tEceV*sc.e)/sc.m_e)

Ni         = pDen
Ne     	   = pDen
plasmaFrequency = 2.*sc.pi*9000.*np.sqrt(2.e-06*Ne)
ds = 0.005 * speedLight / plasmaFrequency
dx = ds 
dy = ds
dz = ds
driftE     = 1e7
# tE         = np.sqrt(2*(tEeV*sc.e)/sc.m_e)  # m/s sqrt(2*(Te*Q)/mE))
vthI         = np.sqrt(2*(tIeV*sc.e)/mAr)   # m/s sqrt(2*(Ti*Q)/mAr)
# dL         = np.sqrt((EPS0*tEceV*KB)/(Ne*Q*Q))  #// [m]
# DL_th      = np.sqrt((sc.epsilon_0*tEK*sc.Boltzmann)/(np.square(sc.e)*Ne))
# DL_lan     = Electron(n=pDen, vth=vthE).debye
# wpE        = np.sqrt((Ne*Q*Q)/(EPS0*mE))
# wpE_1 	   = (8.98E3)*np.sqrt(Ne*1E-6)
# wc         = Q*B0/mE
# gyroRad    = (mE*tE)/(Q*B0)
d = 1. / np.sqrt( 1./(dx*dx) + 1./(dy*dy) + 1./(dz*dz))
timeStep = 0.3 * d / driftE



########## new calculation
t = 10000
n0 = pDen
vthEc       = np.sqrt(2*(tEceV*sc.e)/sc.m_e)
vthEh       = np.sqrt(2*(tEheV*sc.e)/sc.m_e)
vthEb       = np.sqrt(2*(tEbeV*sc.e)/sc.m_e)
# vthEb       = np.sqrt(2*(tEceV*sc.e)/sc.m_e)
dL    = np.sqrt((EPS0*tEceV*KB)/(Ne*Q*Q))
DL_lan     = Electron(n=n0, eV=tEceV).debye
# DL_lan     = Electron(n=pDen, vth=vthEc).debye
LD  = np.sqrt((EPS0*tEceV*eV)/(n0*Q*Q))

alp = 0.04
beta = 1
# nec0 = n0/(1+alp+beta)
# neh0 = alp*nec0
# neb0 = beta*nec0
ni0 = n0
nec0 = 0.4e10
# n0_cal = nec0 + neh0 + neb0
ud  = 20*vthEh
wp = np.sqrt((nec0*Q*Q)/(mE*EPS0))
# timeStep = 0.01/wpE
# t = 32*1.66053907e-27
t = wp*3e-10

# print('nec0 = ', nec0)
# print('neh0 = ', neh0)
# print('neb0 = ', neb0)
# print('n0_cal = ', n0_cal)
print('ud = ', ud)
print('thermal velocity ion= ', vthI)
print('thermal velocity electron cold= ', vthEc)
print('thermal velocity electron hot= ', vthEh)
print('thermal velocity electron beam= ', vthEb)
print('timeStep = ', timeStep)
print('wp = ',wp)
print('t= ',t)
# print('ds = ',ds)
# print('d = ',d)
print('Debye-dL: ', dL)
print('Debye: ', DL_lan)
print('Debye- LD: ', LD)
# print(plasmaFrequency, wpE, wpE_1)


# e_gyro_vs_DL = gyroRad/DL_th
# # dx         = 0.007
# Nx         = 64
# Ny         = 128
# Nz         = 512
# sys_lenx = Nx* dx*DL_th
# sys_leny = Ny* dx*DL_th
# sys_lenz = Nz* dx*DL_th
# probe_size_x = 1*dx*DL_th
# probe_size_y = 59*dx*DL_th
# probe_size_z = 1*dx*DL_th

# print('tE = ',tE)
# print('tI = ',tI)
# print('dL = ',dL)
# print('DL_th = ',DL_th)
# print('wpE = ',wpE)
# print('wpE_1= ',wpE_1)
# print('gyroRad = ',gyroRad)
# print('e_gyro_vs_DL = ',e_gyro_vs_DL)
# print('sys_lenx = ',sys_lenx)
# print('sys_leny = ',sys_leny)
# print('sys_lenz = ',sys_lenz)
# print('probe_size_x = ',probe_size_x)
# print('probe_size_y = ',probe_size_y)
# print('probe_size_z = ',probe_size_z)