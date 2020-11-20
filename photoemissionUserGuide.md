Procedure for photoemission experiments 
========================================
Test case input file: decaetal.ini based on the work by [Deca et al., *Physics of Plasmas* **20**, 102902 (2013)](https://doi.org/10.1063/1.4826951).

1. Git pull trymen/photoemission branch (up to date now)

2. Verify oMode is set as "mode" in methods in .ini


3. IF current density and average energy of ph electron is known:
	3a. Set currentDensity under "Object" in .ini file (amp/m^2) for each object
	3b. Set phCurrentOn switch to 1 under "Object" in .ini file
	3c. Set averageEnergyPH under "Object" in .ini file (Joule) for each object
	3d. Make sure the following two lines are commented out in the objoAlloc function:
		oPlanckPhotonIntegral(ini, units, obj);
    		oPlanckEnergyIntegral(ini, units, obj);
   	3e. Make sure the following line in the objoAlloc function is un-commented:
   		oPhotoElectronCurrent(ini, units, obj);


4. IF only distance from the sun is known (not ph electron temp or ph current density)
    4a. Set distance from sun under "Object" in .ini file (m) for each object
    4b. Set phCurrentOn switch to 0 under "Object" in .ini file
    4c. Set workFunction under "Object" in .ini file (Joule) for material of each object
    4d. Make sure the following two lines are un-commented in the objoAlloc function:
        - oPlanckPhotonIntegral(ini, units, obj);
    	- oPlanckEnergyIntegral(ini, units, obj);
    4e. Make sure the following line in the objoAlloc function is commented:
   	- oPhotoElectronCurrent(ini, units, obj);
    4f. In the function pPhotoElectrons in population.c, set (the parameters in the .ini file haven't been coded yet..)
   		the phYield (# of outgoing ph.electrons per incoming photon)
   		the reflectance (average reflectance of bandwidth as a %)

5. Run mpinc.sh or run_sage.sh

Note: Cell filling algorithm is currently selected for filling adjacent sunlit cells. Swapping  injection algorithms (if wanted) requires a bit more work in the "pPhotoelectrons" function
