import numpy as np

# print(p)
eV = 1.6e-19
Q   = 1.60217662e-19
eps0= 8.85418782e-12
kb  = 1.38e-23
mE  = 9.11E-31
mI  = mE*1836
r = 1e-3  # object radius



nI  = 1e10
nE  = nI

tEeV  = 1.24                  
# tEK   = tEeV*11604.525  
tEK   = 2900       
vthE  = np.sqrt(2*(tEeV*Q)/mE)   
tIeV  = 1.24                
# tIK   = tIeV*11604.525   
tIK = 1800    
vthI  = np.sqrt(2*(tIeV*Q)/mI)   

#To convert workfunction (WF) from electron volts (eV) to wave number (cm⁻¹), the following conversion factor is used: 1 eV = 8065.5 cm⁻¹. Therefore: 
# 4.06 eV * 8065.5 cm⁻¹/eV = 32749.7 cm⁻¹
WF = 4.06 * 8065.5     # = 32749.7 cm⁻¹

surface_area = 4*int(np.pi)*r*r


print('tEK =',tEK )
print('tIK =',tIK )
print('vthE =',vthE )
print('vthI =',vthI )
print('surface_area =',surface_area )


