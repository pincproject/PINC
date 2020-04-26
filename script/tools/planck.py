import scipy.constants as c
from math import exp

time_step = 1e-9
steradians = 5.013025497e-22 
sun_area = 6.088e18
sigma = 33898.31
temp = 5800
reflectance = 0.95


c1 = c.Planck * c.speed_of_light / c.Boltzmann
x = c1 * 100 * sigma / temp
x2 = x * x

iterations = 2.0 + 20.0/x

iterations = iterations if iterations<512 else 512
iterations = int(iterations)

total = 0
for n in range(1, iterations):
    dn = 1.0 / n
    total += exp(-n*x) * (x2 + 2.0*(x + dn)*dn)*dn

kTohc = c.Boltzmann*temp/(c.Planck * c.speed_of_light)
c2 =  2.0 * pow(kTohc,3)*c.speed_of_light
c2 *= total
c2 *= steradians * sun_area * time_step * (1-reflectance)



current_flux = ((c2/time_step) / 6.242e18) / 1.62

print(f"Photon flux per timestep = {c2:e}")
print(f"Current flux in A/m^2 = {current_flux}")