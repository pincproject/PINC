import numpy as np 
import scipy.constants as c
from math import *

n = 7e9
plasma_f = sqrt((n * c.elementary_charge**2) / (c.epsilon_0 * c.electron_mass))
print(1 / plasma_f)