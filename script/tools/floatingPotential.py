'''
author: Trym Erik Nielsen
date: 15.05.20
purpose: Compute the floating potential
of a cylindrical spacecraft for a given
photoelemission current and plasma density
assumptions: Charging contribution from
ambient ions can be ignored (usually can be
for positively charged spacecraft)
'''
import scipy.constants as c
import math as m

def main():
    #Cylinder input
    a = 0.9
    L = 0.9
    phCurrent = 0.0261 #Amp


    #plasma input
    n = 1e7
    T = 3.667e6 #K
    v = 1.5e5 #drift m/s

    I_0 = 2 * m.pi * a * n * (-1 * c.elementary_charge) * v * L 
    print(I_0)
    floatV = ((c.Boltzmann * T)/(-1 * c.elementary_charge)) * (1 - m.pow((0.0261/I_0),2));
    print(floatV)



if __name__ == '__main__':
    main()