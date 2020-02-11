import math as m
import scipy.constants as c

#light stuff
temp = 5800 #sun blackbody temp
sigma = 50000#36364 #cutoff wavenumber
area = 1 #m^2
tstep = 1e-9
dist_sun = 9.0862e10 #m
sun_area = 6.1e18 #m^2
reflectance = 0.0

sr = area / (dist_sun**2) #steradians
tstep = 1e-9

planck = c.Planck
boltz = c.Boltzmann
light_speed = c.speed_of_light
light_speed_sq = light_speed**2

c1 = planck * light_speed / boltz
x = c1 * 100 * sigma / temp
x2 = x**2

iterations = 2.0 + 20.0/x
iters = iterations if iterations<512 else 512
iters = int(iters)
print(iters)

total = 0
for n in range(1,iters):
    dn = 1.0/n
    total += m.exp(-n*x) * (x2 + 2.0*(x + dn)*dn)*dn

ktohc = boltz*temp/(planck*light_speed)
c2 = 2.0*pow(ktohc,3)*light_speed
c2 *= total

radiance = c2 * sr * sun_area * tstep * (1 - reflectance)
print('{:.5e}'.format(radiance))