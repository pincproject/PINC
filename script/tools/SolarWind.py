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
    
rad_dist = 0.607373 #float(sys.argv[1])
a1d = 0.7766; a2d = -1.934; a3d = 0.01823; a4d = -2.245
a1v = 2.651; a2v = 0.0; a3v = -0.0239; a4v = -1.836;
a1t = 4.848; a2t = -0.668; a3t = -4.69e-42; a4t = -40.91;

X = m.log10(rad_dist)
Yd = a1d + a2d * X + a3d * m.exp(a4d * X)
Yv = a1v + a2v * X + a3v * m.exp(a4v * X)
Yt = a1t + a2t #X + a3t * m.exp(a4t * X)

#temp = (Yt**2 * c.electron_mass) / c.Boltzmann #in Kelvin

num_density = m.pow(10,Yd) #density in cm^-3
drift = m.pow(10,Yv) #drift velocity in km/s
temp = m.pow(10,Yt) #temp in kelvin

thermal_vel_e = m.sqrt((2 * c.Boltzmann * temp) / c.electron_mass) #in m/s
thermal_vel_p = m.sqrt((2 * c.Boltzmann * temp) / c.proton_mass) #in m/s

debye = m.sqrt((c.epsilon_0 * c.Boltzmann * temp) / (c.elementary_charge**2 * num_density * 1e6))
print "Numerical density of plasma at %f AU is %f cm^-3 \n" % (rad_dist, num_density)
print "Drift velocity of plasma at %f is %f km/s \n" % (rad_dist, drift)
print "Thermal velocity of plasma (electrons) at %f is %f m/s" % (rad_dist, thermal_vel_e)
print "Thermal velocity of plasma (protons) at %f is %f m/s" % (rad_dist, thermal_vel_p)
print "Debye length of plasma at %f AU is %f meters \n" % (rad_dist, debye)
print "Plasma temperature at %f AU is %f Kelvin \n" %(rad_dist, temp)