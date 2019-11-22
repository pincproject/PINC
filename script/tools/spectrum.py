"""
This tool computes the blackbody solar irradiance in vacuum 
in terms of photons/m^2/s 
"""
import scipy.constants as c
import numpy as np

k = c.Boltzmann
c = c.speed_of_light
h = c.Planck

def plancs(wav,temp):
    a = 2.0*h*c**2
    b = h*c/(wav*k*T)
    energy_photon = (h*c)/wav
    intensity = a / ( (wav**5) * (np.exp(b) - 1.0) )
    return intensity

def main():
  wave_length = np.linspace(0,2000,20000) #in nanometer
  temp = 5000 #in kelvin



if __name__ == '__main__':
  main()
