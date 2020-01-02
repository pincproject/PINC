# @file			SolarWind.py
# @author		Trym Erik Nielsen <trymen@math.uio.no>
# @copyright	University of Oslo, Norway
# @brief		Solar wind plasma parameter tool
# @date			24.10.19
#
# This tool computes plasma parameters based on the article
# "Radial dependence of solar wind parameters in the ecliptic"
# https://doi.org/10.1007/BF00153841

import sys
import math as m
import scipy.constants as c
    
rad_dist = float(sys.argv[1]) #0.406821#0.607373#
a1d = 0.7766; a2d = -1.934; a3d = 0.01823; a4d = -2.245
a1v = 2.651; a2v = 0.0; a3v = -0.0239; a4v = -1.836
a1t = 4.848; a2t = -0.668; a3t = -4.69e-42; a4t = -40.91

X = m.log10(rad_dist)
Yd = a1d + a2d * X + a3d * m.exp(a4d * X)
Yv = a1v + a2v * X + a3v * m.exp(a4v * X)
Yt = a1t + a2t * X + a3t * m.exp(a4t * X)
#temp = (Yt**2 * c.electron_mass) / c.Boltzmann #in Kelvin

num_density = m.pow(10,Yd) * 1e6 #density in m^-3
drift = m.pow(10,Yv) #drift velocity in km/s
temp = m.pow(10,Yt) #temp in kelvin

thermal_vel_e = m.sqrt((3 * c.Boltzmann * temp) / c.electron_mass) #in m/s
thermal_vel_p = m.sqrt((3 * c.Boltzmann * temp) / c.proton_mass) #in m/s

gyro_freq = m.sqrt((num_density * c.elementary_charge**2)/(c.electron_mass * c.epsilon_0)) #in Hertz
debye = m.sqrt((c.epsilon_0 * c.Boltzmann * temp) / (c.elementary_charge**2 * 4 * m.pi * num_density))

print "Numerical density of plasma at %f AU is %f m^-3 \n" % (rad_dist, num_density)
print "Electron plasma frequency at %f AU is %f sec" % (rad_dist, 1./gyro_freq)
print "Drift velocity of plasma at %f AU is %f km/s \n" % (rad_dist, drift)
print "Thermal velocity of plasma (electrons) at %f AU is %f km/s" % (rad_dist, thermal_vel_e/1000)
print "Thermal velocity of plasma (protons) at %f AU is %f km/s" % (rad_dist, thermal_vel_p/1000)
print "Debye length of plasma at %f AU is %f meters \n" % (rad_dist, debye)
print "Plasma temperature at %f AU is %f Kelvin \n" %(rad_dist, temp)