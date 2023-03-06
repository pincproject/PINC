import numpy as np
import math
from langmuir import *
import scipy.constants as sc


Q          = 1.6e-19
KB         = 1.38e-23
EPS0       = 8.85418782e-12
mAr	       = 1.6724828e-27
mE         = 9.10938356e-31
tEeV       = 1.6
tIeV	   = 0.1
B0         = 5e-5
# tEK        = tEeV*11604.525
# tIK        = tIeV*11604.525
tEK        = 1000
tIK        = 1000
pDen       = 5.8977e9
Ni         = pDen
Ne     	   = pDen
driftE     = 2873
tE         = np.sqrt(2*(tEeV*sc.e)/sc.m_e)  # m/s sqrt(2*(Te*Q)/mE)
tI         = np.sqrt(2*(tIeV*sc.e)/mAr)   # m/s sqrt(2*(Ti*Q)/mAr)
dL         = np.sqrt((EPS0*tEK*KB)/(Ne*Q*Q))  #// [m]
DL_th      = np.sqrt((sc.epsilon_0*tEK*sc.Boltzmann)/(np.square(sc.e)*Ne))
DL_lan     = Electron(n=pDen, eV=tEeV).debye
wpE        = np.sqrt((Ne*Q*Q)/(EPS0*mE))
wpE_1 	   = (8.98E3)*np.sqrt(Ne*1E-6)
wc         = Q*B0/mE
gyroRad    = (mE*tE)/(Q*B0)


e_gyro_vs_DL = gyroRad/DL_th
dx         = 0.007
Nx         = 64
Ny         = 128
Nz         = 512
sys_lenx = Nx* dx*DL_th
sys_leny = Ny* dx*DL_th
sys_lenz = Nz* dx*DL_th
probe_size_x = 1*dx*DL_th
probe_size_y = 59*dx*DL_th
probe_size_z = 1*dx*DL_th

print('tE = ',tE)
print('tI = ',tI)
print('dL = ',dL)
print('DL_th = ',DL_th)
print('wpE = ',wpE)
print('wpE_1= ',wpE_1)
print('gyroRad = ',gyroRad)
print('e_gyro_vs_DL = ',e_gyro_vs_DL)
print('sys_lenx = ',sys_lenx)
print('sys_leny = ',sys_leny)
print('sys_lenz = ',sys_lenz)
print('probe_size_x = ',probe_size_x)
print('probe_size_y = ',probe_size_y)
print('probe_size_z = ',probe_size_z)
