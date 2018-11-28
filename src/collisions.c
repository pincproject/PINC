/**
* @file		 	collisions.c
* @brief		Collisional module, Montecarlo Collision Method.
* @author		Steffen Mattias Brask <steffen.brask@fys.uio.no>
*
* Main module collecting functions conserning collisions
*
* MCC details......
*
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_randist.h>

#include "core.h"
#include "collisions.h"
#include "multigrid.h"
#include "pusher.h"

static void mccSanity(dictionary *ini, const char* name, int nSpecies){
	// check error, see vahedi + surendra p 181
	// check for v_ion + v_neutral > maxVel (moves beyond a cell in one timestep)
	int countSpecies = iniGetInt(ini,"population:nSpecies");
	if(countSpecies!=nSpecies){
		msg(ERROR,"%s only supports 2 species",name);
	}
}


/*************************************************
*	Normalization MCC specific variables
************************************************/

static void mccNormalize(dictionary *ini,const Units *units){

	//make shure *units is already normalized!


	//collision frequency given in 1/s
	double collFrqElectronElastic = iniGetDouble(ini,"collisions:collFrqElectronElastic");
	double collFrqIonElastic = iniGetDouble(ini,"collisions:collFrqIonElastic");
	double collFrqCEX = iniGetDouble(ini,"collisions:collFrqCEX");

	collFrqElectronElastic /= units->frequency;
	collFrqIonElastic /= units->frequency;
	collFrqCEX /= units->frequency;

	iniSetDouble(ini,"collisions:collFrqElectronElastic",collFrqElectronElastic);
	iniSetDouble(ini,"collisions:collFrqIonElastic",collFrqIonElastic);
	iniSetDouble(ini,"collisions:collFrqCEX",collFrqCEX);

	//numberdensityneutrals given as numberofparticles/m^3
	double nt = iniGetDouble(ini,"collisions:numberDensityNeutrals");

	nt /= units->density; // assumes same for elecron and ion
	nt /= units->weights[1];
	//msg(STATUS,"neutral density = %f",nt);
	//we use computational particles that contain many real particles.
	iniSetDouble(ini,"collisions:numberDensityNeutrals",nt);

	// in m/s
	double NvelThermal = iniGetDouble(ini,"collisions:thermalVelocityNeutrals");
	NvelThermal /= (units->length/units->time);
	iniSetDouble(ini,"collisions:thermalVelocityNeutrals",NvelThermal);

	// cross sections given in m^2
	double StaticSigmaCEX = iniGetDouble(ini,"collisions:sigmaCEX");
	double StaticSigmaIonElastic = iniGetDouble(ini,"collisions:sigmaIonElastic");
	double StaticSigmaElectronElastic = iniGetDouble(ini,"collisions:sigmaElectronElastic");

	StaticSigmaCEX /= (units->length*units->length);
	StaticSigmaIonElastic /= (units->length*units->length);
	StaticSigmaElectronElastic /= (units->length*units->length);
	//cross section for computational particles
	StaticSigmaCEX *= units->weights[1];
	StaticSigmaIonElastic *= units->weights[1];
	StaticSigmaElectronElastic *= units->weights[1];

	iniSetDouble(ini,"collisions:sigmaCEX",StaticSigmaCEX);
	iniSetDouble(ini,"collisions:sigmaIonElastic",StaticSigmaIonElastic);
	iniSetDouble(ini,"collisions:sigmaElectronElastic",StaticSigmaElectronElastic);

	// a parameter is max crossect (m^2)
	// b parameter decides velocity to center about.
	// b given as 1/v^2

	double CEX_a = iniGetDouble(ini,"collisions:CEX_a");
	double CEX_b = iniGetDouble(ini,"collisions:CEX_b");
	double ion_elastic_a = iniGetDouble(ini,"collisions:ion_elastic_a");
	double ion_elastic_b = iniGetDouble(ini,"collisions:ion_elastic_b");
	double electron_a = iniGetDouble(ini,"collisions:electron_a");
	double electron_b = iniGetDouble(ini,"collisions:electron_b");

	CEX_a /= (units->length*units->length);
	ion_elastic_a /= (units->length*units->length);
	electron_a /= (units->length*units->length);

	CEX_a *= units->weights[1];
	ion_elastic_a *= units->weights[1];
	electron_a *= units->weights[1];

	CEX_b /= ((units->length/units->time)*(units->length/units->time));
	ion_elastic_b /= ((units->length/units->time)*(units->length/units->time));
	electron_b /= ((units->length/units->time)*(units->length/units->time));

	iniSetDouble(ini,"collisions:CEX_a",CEX_a);
	iniSetDouble(ini,"collisions:CEX_b",CEX_b);
	iniSetDouble(ini,"collisions:ion_elastic_a",ion_elastic_a);
	iniSetDouble(ini,"collisions:ion_elastic_b",ion_elastic_b);
	iniSetDouble(ini,"collisions:electron_a",electron_a);
	iniSetDouble(ini,"collisions:electron_b",electron_b);


}

/*************************************************
*	Search functions, can probably use other Functions
************************************************/

double mccGetMaxVel(const Population *pop, int species){

	// iterate over pop and keep max velocity
	double *vel = pop->vel;
	double MaxVelocity = 0;
	double NewVelocity = 0;
	int nDims = pop->nDims;

	long int iStart = pop->iStart[species];
	long int iStop  = pop->iStop[species];
	for(int i=iStart; i<iStop; i++){
		NewVelocity = sqrt(vel[i*nDims]*vel[i*nDims]+vel[i*nDims+1]\
			*vel[i*nDims+1]+vel[i*nDims+2]*vel[i*nDims+2]);
		if(NewVelocity>MaxVelocity){
			MaxVelocity=NewVelocity;
		}
		// if(NewVelocity>32){
		// 	vel[i*nDims] = 1;
		// 	vel[i*nDims+1] = 1;
		// 	vel[i*nDims+2] = 1;
		// 	msg(WARNING,"velocities greater than 32 encountered. velocities suppressed");
		// }
	}
	return MaxVelocity;
}
double mccGetMinVel(const Population *pop, int species){

	// iterate over pop and keep min velocity
	double *vel = pop->vel;
	double MinVelocity = 100000000000000;
	double NewVelocity = 0;
	int nDims = pop->nDims;

	long int iStart = pop->iStart[species];
	long int iStop  = pop->iStop[species];
	for(int i=iStart; i<iStop; i++){
		NewVelocity = sqrt(vel[i*nDims]*vel[i*nDims]+vel[i*nDims+1]\
			*vel[i*nDims+1]+vel[i*nDims+2]*vel[i*nDims+2]);
		if(NewVelocity<MinVelocity){
			MinVelocity=NewVelocity;
		}
	}
	return MinVelocity;
}


/*************************************************
*		Inline functions
************************************************/

inline double mccSigmaCEXFunctional(double a,double b,double v){
	// educated guess on funtional form
	double sigma = a*exp(-b*v*v);
	//msg(STATUS,"in mccSigmaCexF a=%f, b=%f, v = %f",a,b,v );
	return sigma;//2*0.00121875*pow(10,8);//0.5*pow(10,-1);//temp;
}
inline double mccSigmaIonElasticFunctional(double a,double b,double v){
	// educated guess on funtional form
	double sigma = a*exp(-b*v*v);
	return sigma;
}
inline double mccSigmaElectronElasticFunctional(double a,double b,double v){
	// educated guess on funtional form
	double sigma = a*exp(-b*v*v);
	return sigma;
}

inline double mccGetMyCollFreqStatic(double sigma_T, double vx,double vy,double vz, double nt){
	// collision freq given a static (double) cross section
	// this means nu \propto v
	double v = sqrt(vx*vx + vy*vy + vz*vz);
	double MyFreq = v*sigma_T*nt; //multiplied by constant sigma
	return MyFreq;
}
inline double mccGetMyCollFreqFunctional(double (*sigma)(double, double,double),
						double vx,double vy,double vz,double nt,
						double a,double b){

	// collision freq given a cross section function
	double v = sqrt(vx*vx + vy*vy + vz*vz);
	//we use v instead of Energy
	double newSigma = sigma(a,b,v);
	double newFreq = v*newSigma*nt; //multiplied by functional sigma(v)
	return newFreq;
}

// inline double mccGetMyCollFreq(const dictionary *ini, double (*sigma)(const dictionary *ini, double),
// 						double m, double vx,double vy,double vz, double dt, double nt){
// 	// collision freq given a cross section function
// 	double v = sqrt(vx*vx + vy*vy + vz*vz);
// 	double eps_i = 0.5*m*v*v;
// 	double sigma_T = sigma(ini,eps_i);
// 	double MyFreq = v*sigma_T*nt;
// 	//msg(STATUS,"Myfreq =  %f", MyFreq);
// 	//double P = 1-exp(-(MyFreq*0.04));
// 	//msg(STATUS,"My P =  %f", P);
// 	return MyFreq;
// }


/*************************************************
*		Allocation mcc variables
************************************************/


MccVars *mccAlloc(const dictionary *ini, const Units *units){
	MccVars *mccVars = malloc(sizeof(*mccVars));
	double pMaxElectron = 0;
	double maxFreqElectron = 0;
	double pMaxIon = 0;
	double maxFreqIon = 0;
	double electronMassRatio = 0;
	// double artificialLoss = 1.;
	mccVars->pMaxElectron = pMaxElectron;
	mccVars->maxFreqElectron = maxFreqElectron;
	mccVars->pMaxIon = pMaxIon;
	mccVars->maxFreqIon = maxFreqIon;
	mccVars->artificialLoss = iniGetDouble(ini,"collisions:artificialLoss");
	// double nt = iniGetDouble(ini,"collisions:numberDensityNeutrals"); //constant for now
	// double NvelThermal = iniGetDouble(ini,"collisions:thermalVelocityNeutrals");
	// double mccSigmaElectronElastic = iniGetDouble(ini,"collisions:sigmaElectronElastic");
	// mccVars->nt = nt;
	// mccVars->NvelThermal = NvelThermal;
	// mccVars->mccSigmaElectronElastic = mccSigmaElectronElastic;
	int nSpecies = iniGetInt(ini, "population:nSpecies");
	double *mass = iniGetDoubleArr(ini, "population:mass", nSpecies);
	electronMassRatio = mass[0]*units->mass/units->weights[0];
	electronMassRatio /= iniGetDouble(ini,"collisions:realElectronMass");
	mccVars->nt = iniGetDouble(ini,"collisions:numberDensityNeutrals"); //constant for now
	mccVars->NvelThermal = iniGetDouble(ini,"collisions:thermalVelocityNeutrals");
	mccVars->mccSigmaElectronElastic = iniGetDouble(ini,"collisions:sigmaElectronElastic");
	mccVars->mccSigmaCEX= iniGetDouble(ini,"collisions:sigmaCEX");
	mccVars->mccSigmaIonElastic = iniGetDouble(ini,"collisions:sigmaIonElastic");
	mccVars->energyConvFactor = (units->energy/6.24150913*pow(10,18)); //J/(J/eV)

	mccVars->collFrqCex = iniGetDouble(ini,"collisions:collFrqCex");
	mccVars->collFrqIonElastic = iniGetDouble(ini,"collisions:collFrqIonElastic");
	mccVars->collFrqElectronElastic = iniGetDouble(ini,"collisions:collFrqElectronElastic");

	mccVars->CEX_a = iniGetDouble(ini,"collisions:CEX_a");
	mccVars->CEX_b = iniGetDouble(ini,"collisions:CEX_b");
	mccVars->ion_elastic_a = iniGetDouble(ini,"collisions:ion_elastic_a");
	mccVars->ion_elastic_b = iniGetDouble(ini,"collisions:ion_elastic_b");
	mccVars->electron_a = iniGetDouble(ini,"collisions:electron_a");
	mccVars->electron_b = iniGetDouble(ini,"collisions:electron_b");
	//msg(STATUS,"weigh = %.29f",electronMassRatio);
	//msg(STATUS,"factor = %f",mccVars->energyConvFactor);

	free(mass);
	return mccVars;
}

void mccFreeVars(MccVars *mccVars){
	free(mccVars);

}

/*************************************************
*		Max collision probability Functions
************************************************/

void mccGetPmaxElectronConstantFrq(const dictionary *ini,
	MccVars *mccVars,Population *pop,MpiInfo *mpiInfo){
	//
	double collFrqElectronElastic = mccVars->collFrqElectronElastic;

	double max_v = mccGetMaxVel(pop,0); //2.71828*thermalVel; // e*thermalVel, needs to be max_velocity function
	//msg(STATUS,"maxVelocity Electron =  %f", max_v);
	mccVars->pMaxElectron = 1-exp(-((collFrqElectronElastic)));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "P_coll Electron = %f \n", mccVars->pMaxElectron);
		fMsg(ini, "collision", "max velocity Electron = %f \n", max_v);

	}
}

void mccGetPmaxIonConstantFrq(const dictionary *ini,MccVars *mccVars,
	Population *pop,MpiInfo *mpiInfo){
	//
	double collFrqIonElastic = mccVars->collFrqIonElastic;
	double collFrqCEX = mccVars->collFrqCex;
	double max_v = mccGetMaxVel(pop,1); //2.71828*thermalVel; // e*thermalVel, needs to be max_velocity function
	//msg(STATUS,"maxVelocity Ion =  %f", max_v);
	mccVars->pMaxIon = 1-exp(-((collFrqIonElastic+collFrqCEX)));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "P_coll Ion = %f \n", mccVars->pMaxIon);
		fMsg(ini, "collision", "max velocity Ion = %f \n", max_v);
	}
}

void mccGetPmaxIonFunctional(const dictionary *ini,
	MccVars *mccVars, Population *pop,
	MpiInfo *mpiInfo){
	// determine max local number desity / max_x(n_t(x_i)) (for target species, neutrals)
	// Get å functional form of sigma_T (total crossect as funct. of energy)
	// determine the speed, v(eps_i) = sqrt(2*eps_i *m_s)
	// determine, max_eps(sigma_T *v)

	double nt = mccVars->nt;
	double CEX_a = mccVars->CEX_a;
	double CEX_b = mccVars->CEX_b;
	double elastic_a = mccVars->ion_elastic_a;
	double elastic_b = mccVars->ion_elastic_b;

	double max_v = mccGetMaxVel(pop,1);

	double *vel = pop->vel;
	double NewFreq = 0;
	double NewVelocity = 0;
	int nDims = pop->nDims;
	long int iStart = pop->iStart[1]; //ions as specie 1
	long int iStop  = pop->iStop[1];
	for(int i=iStart; i<iStop; i++){
		NewVelocity = sqrt(vel[i*nDims]*vel[i*nDims]+vel[i*nDims+1]\
			*vel[i*nDims+1]+vel[i*nDims+2]*vel[i*nDims+2]);
		NewFreq = (mccSigmaCEXFunctional(CEX_a,CEX_b,NewVelocity)
			+mccSigmaIonElasticFunctional(elastic_a,elastic_b,NewVelocity))*NewVelocity*nt;
		if(NewFreq>mccVars->maxFreqIon){
			mccVars->maxFreqIon=NewFreq;
		}
	}
	mccVars->pMaxIon = 1-exp(-(mccVars->maxFreqIon));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision","getPmax Ion =  %f \n", mccVars->pMaxIon);
		fMsg(ini, "collision", "max velocity Ion = %f \n", max_v);
		//fMsg(ini, "collision","dt =  %f \n", dt);
	}
}

void mccGetPmaxElectronFunctional(const dictionary *ini,
	MccVars *mccVars, Population *pop,
	MpiInfo *mpiInfo){
	// determine max local number desity / max_x(n_t(x_i)) (for target species, neutrals)
	// Get a functional form of sigma_T (total crossect as funct. of energy)
	// determine the speed, v(eps_i) = sqrt(2*eps_i *m_s)
	// determine, max_eps(sigma_T *v)

	double nt = mccVars->nt;
	double a = mccVars->electron_a;
	double b = mccVars->electron_b;

	double max_v = mccGetMaxVel(pop,0);
	double *vel = pop->vel;
	double NewFreq = 0;
	double NewVelocity = 0;
	int nDims = pop->nDims;
	long int iStart = pop->iStart[0]; //Electrons as specie 0
	long int iStop  = pop->iStop[0];
	for(int i=iStart; i<iStop; i++){
		NewVelocity = sqrt(vel[i*nDims]*vel[i*nDims]+vel[i*nDims+1]\
			*vel[i*nDims+1]+vel[i*nDims+2]*vel[i*nDims+2]);
		NewFreq = mccSigmaElectronElasticFunctional(a,b,NewVelocity)*NewVelocity*nt;
		if(NewFreq>mccVars->maxFreqElectron){
			mccVars->maxFreqElectron=NewFreq;
		}
	}
	mccVars->pMaxElectron = 1-exp(-(mccVars->maxFreqElectron));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision","getPmax electron =  %f \n", mccVars->pMaxElectron);
		fMsg(ini, "collision", "max velocity electron = %f \n", max_v);
	}
}

void mccGetPmaxIonStatic(const dictionary *ini,MccVars *mccVars,
	Population *pop, MpiInfo *mpiInfo){

	// Faster static version. uses static cross sections

	// determine max local number desity / max_x(n_t(x_i)) (for target species, neutrals)
	// Get å functional form of sigma_T (total crossect as funct. of energy)
	// determine the speed, v(eps_i) = sqrt(2*eps_i *m_s)
	// determine, max_eps(sigma_T *v)

	double nt = mccVars->nt;
	double StaticSigmaCEX = mccVars->mccSigmaCEX;
	double StaticSigmaIonElastic = mccVars->mccSigmaIonElastic;
	double max_v = mccGetMaxVel(pop,1);
	mccVars->maxFreqIon = (StaticSigmaCEX +StaticSigmaIonElastic)*max_v*nt;
	mccVars->pMaxIon = 1-exp(-(mccVars->maxFreqIon));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision","getPmax Ion =  %f \n", mccVars->pMaxIon);
		fMsg(ini, "collision", "max velocity Ion = %f \n", max_v);
		//fMsg(ini, "collision","dt =  %f \n", dt);
	}
}

void mccGetPmaxElectronStatic(const dictionary *ini,
	MccVars *mccVars, Population *pop, MpiInfo *mpiInfo){

	// Faster static version. uses static cross sections

	// determine max local number desity / max_x(n_t(x_i)) (for target species, neutrals)
	// Get a functional form of sigma_T (total crossect as funct. of energy)
	// determine the speed, v(eps_i) = sqrt(2*eps_i *m_s)
	// determine, max_eps(sigma_T *v)

	double nt = mccVars->nt;//iniGetDouble(ini,"collisions:numberDensityNeutrals"); //constant for now
	double StaticSigmaElectronElastic = mccVars->mccSigmaElectronElastic;//iniGetDouble(ini,"collisions:sigmaElectronElastic");
	double max_v = mccGetMaxVel(pop,0);//2.71828*thermalVel; // e*thermalVel, needs to be max_velocity
	double min_v = mccGetMinVel(pop,0);
	mccVars->maxFreqElectron=StaticSigmaElectronElastic*max_v*nt;
	mccVars->pMaxElectron = 1-exp(-(mccVars->maxFreqElectron));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision","getPmax electron =  %f \n", mccVars->pMaxElectron);
		fMsg(ini, "collision", "max velocity electron = %f \n", max_v);
		fMsg(ini, "collision", "min velocity electron = %f \n", min_v);
	}
}

/*************************************************
*		Collision Functions
************************************************/

void scatterElectron(double *vx_point, double *vy_point,double *vz_point,
	const gsl_rng *rng, MccVars *mccVars, Population *pop){

	double artificialLoss = mccVars->artificialLoss;
	double *mass = pop->mass;

	double angleChi = 0; //scattering angle in x, y
	double anglePhi = 0; //scattering angle in j
	double angleTheta = 0;
	double A = 0; //scaling Factor used for precomputing (speedup)
	double B = 0;
	double Ekin, newEkin;
	//R1 = gsl_rng_uniform_pos(rng); // New random number per particle. maybe print?
	double R2 = gsl_rng_uniform_pos(rng);
	double R3 = gsl_rng_uniform_pos(rng);
	double vx = *vx_point;
	double vy = *vy_point;
	double vz = *vz_point;
	double vx_,vy_,vz_;

	// unit velocities
	//double velscatter = 0.0;
	double velchange = 0.0;// change of energy incident particle
	double velsquare = sqrt(vx*vx+vy*vy+vz*vz);
	//msg(STATUS, "velsquare = %f",velsquare);
	if(velsquare<0.0000000000000001){
		velsquare = 0.0000000000000001;
	}

	Ekin = mccVars->energyConvFactor*0.5*(vx*vx + vy*vy + vz*vz)*mass[0]; //eV

	//make unit vector
	vx = vx/velsquare;
	vy = vy/velsquare;
	vz = vz/velsquare;

	// angles
	double argument = (2+Ekin-2*pow((1+Ekin),R2))/(Ekin); // This one has the posibility to be greater than 1!!
	Ekin = Ekin/(mccVars->energyConvFactor); //PINC energy

	if(sqrt(argument*argument)>1.0){
		double newargument = sqrt(argument*argument)/argument;
		//fMsg(ini, "collision", "WARNING: old argument acos = %.32f new arg = %.64f, R = %.32f \n", argument, newargument,R2);
		msg(WARNING,"old argument acos = %.32f new arg = %.64f, R = %.32f", argument, newargument,R2);
		argument = newargument;
	}
	//msg(STATUS, "energy is %f energy in eV is %f",Ekin,(Ekin*mccVars->energyConvFactor));
	angleChi =  acos(argument); // gives nan value if abs(argument) > 1
	//msg(STATUS, "angleChi = %f",angleChi);
	newEkin = Ekin*(1-((2.0*mass[0])/(artificialLoss*mass[1]))*(1-(cos(angleChi))));

	if(((2.0*(newEkin))/(mass[0]))<0.0){
		velchange = sqrt( abs((2.0*(newEkin))/(mass[0])));
		msg(STATUS, "corrected velchange = %.32f",velchange);
		msg(STATUS, "Ekin = %f,newEkin = %f",Ekin,newEkin);
	}else{
		velchange = sqrt( (2.0*(newEkin))/(mass[0] ));
	}

	anglePhi = 2*PI*R3;
	if(vx>0.0){
		if(vx<0.0000000000000001){vx=0.0000000000000001;}
		if(vx>0.99999999999999){vx=0.99999999999999;}
	}
	if(vx<=0.0){
		if(vx>-0.0000000000000001){vx=-0.0000000000000001;}
		if(vx<-0.99999999999999){vx=-0.99999999999999;}
	}
	angleTheta = acos(vx);
	//msg(STATUS, "angleTheta = %f",angleTheta);
	// inc - scat vector relation
	//vector relation
	A = (sin(angleChi)*sin(anglePhi))/sin(angleTheta);
	B = (sin(angleChi)*cos(anglePhi))/sin(angleTheta);
	vx_ = (vx*cos(angleChi)+B*(vy*vy + vz*vz) ); //Vx
	vy_ = (vy*cos(angleChi)+A*vz-B*vx*vy); //Vy
	vz_ = (vz*cos(angleChi)-A*vy-B*vx*vz); //Vz

	//msg(STATUS, "should be unity %f",sqrt(vx_*vx_+ vy_*vy_+ vz_*vz_));
	// change in energy
	vx_ *= velchange;
	vy_ *= velchange;
	vz_ *= velchange;

	*vx_point=vx_;
	*vy_point=vy_;
	*vz_point=vz_;

	//msg(STATUS, "last vx_,vy_vz_ %f,%f,%f",vx_,vy_,vz_);

	// unittest

	double energydiffexpr = Ekin*(((2.0*mass[0])/(artificialLoss*mass[1]))*(1-(cos(angleChi))));
	double energydiff = Ekin-0.5*(vx_*vx_ + vy_*vy_ + vz_*vz_)*mass[0];
	if(abs(energydiffexpr-energydiff)>1e-31){
		msg(WARNING,"too large energy error in collide electrons");
	}
	if(Ekin<newEkin){
		msg(WARNING,"energy increased in electron collission!");
	}

}

void scatterIon(double *vx_point, double *vy_point,double *vz_point,
	const gsl_rng *rng, MccVars *mccVars, Population *pop){


	double NvelThermal = mccVars->NvelThermal;
	double *mass = pop->mass;

	double angleChi = 0; //scattering angle in x, y
	double anglePhi = 0; //scattering angle in j
	double angleTheta = 0;
	double A = 0; //scaling Factor used for precomputing (speedup)
	double B = 0;
	double Ekin, newEkin;
	double velchange;
	//double R1 = gsl_rng_uniform_pos(rng);

	double vx = *vx_point;
	double vy = *vy_point;
	double vz = *vz_point;
	double vx_=0;
	double vy_=0;
	double vz_=0;
	double R1 = gsl_rng_uniform_pos(rng);
	double R2 = gsl_rng_uniform_pos(rng);

	//The scattering happens in netral rest frame
	//store (to transfer back)
	double vxMW, vyMW, vzMW;
	vxMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal); //maxwellian dist?
	vyMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal); //yes, gaussian in 3 dim
	vzMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal); // is maxwellian


	//transfer to neutral stationary frame
	double vxTran = vx-vxMW;
	double vyTran = vy-vyMW;
	double vzTran = vz-vzMW;
	//msg(STATUS,"vx =%f,vy =%f,vz =%f",vx,vy,vz);
	//msg(STATUS,"vxTran =%f,vyTran =%f,vzTran =%f",vxTran,vyTran,vzTran);

	double velsquare = sqrt(vxTran*vxTran+vyTran*vyTran+vzTran*vzTran);
	if(velsquare<0.00000000001){
		velsquare = 0.00000000001;
	}

	//msg(STATUS,"velocity before = %f",sqrt(vx*vx+vy*vy+vz*vz));
	//make unit vector
	vx = vxTran/velsquare;
	vy = vyTran/velsquare;
	vz = vzTran/velsquare;

	// angles
	Ekin = 0.5*(vxTran*vxTran+vyTran*vyTran+vzTran*vzTran)*mass[1]; //PINC neutral frame

	// verboncoeur says use THETHA (CM frame)
	angleChi =  acos(sqrt(1-R1)); //labframe goes from 0 to pi/2
	//angleChi =  0.5*acos(1-2*R1); //labframe goes from 0 to pi/2
	//angleChi =  acos(1-2*R1); //CM-frame goes from 0 to pi
	//msg(STATUS,"angleChi = %f",angleChi);

	newEkin = Ekin*(cos(angleChi)*cos(angleChi));
	//newEkin = Ekin*(0.5+0.5*cos(angleChi) );

	//msg(STATUS,"old enrgy = %f, new energy = %f, difference should be %f, and is %f",Ekin,newEkin,Ekin*(((2.0*mass[0])/(artificialLoss*mass[1]))*(1-(cos(angleChi)))),(Ekin-newEkin));
	//msg(STATUS, "loss factor =  %f",(((2.0*mass[0])/(artificialLoss*mass[1]))*(1-(cos(angleChi)))));
	if(((2.0*(newEkin))/(mass[1]))<0.0){
		velchange = sqrt( abs((2.0*(newEkin))/(mass[1])));
		msg(STATUS, "corrected velchange = %.32f",velchange);
		msg(STATUS, "Ekin = %f,newEkin = %f",Ekin,newEkin);
	}else{
		velchange = sqrt( (2.0*newEkin)/mass[1] );
	}

	//velchange = sqrt(velsquare*velsquare*(cos(angleChi)*cos(angleChi)));

	anglePhi = 2*PI*R2;
	if(vx>0.0){
		if(vx<0.0000000000000001){vx=0.0000000000000001;}
		if(vx>0.99999999999999){vx=0.99999999999999;}
	}
	if(vx<=0.0){
		if(vx>-0.0000000000000001){vx=-0.0000000000000001;}
		if(vx<-0.99999999999999){vx=-0.99999999999999;}
	}
	angleTheta = acos(vx);

	//vector relation
	// vx, vy, vz must form a unit vector!
	A = (sin(angleChi)*sin(anglePhi))/sin(angleTheta);
	B = (sin(angleChi)*cos(anglePhi))/sin(angleTheta);
	vx_ = (vx*cos(angleChi)+B*(vy*vy + vz*vz) ); //Vx
	vy_ = (vy*cos(angleChi)+A*vz-B*vx*vy); //Vy
	vz_ = (vz*cos(angleChi)-A*vy-B*vx*vz); //Vz

	// unittest
	//msg(STATUS, "should be unity old %f",sqrt(vx*vx+vy*vy+vz*vz));
	//msg(STATUS, "should be unity %f",sqrt(vx_*vx_+vy_*vy_+vz_*vz_));

	//New velocities (change in energy)
	vx_ *= velchange;
	vy_ *= velchange;
	vz_ *= velchange;

	//msg(STATUS,"velocity after = %f",sqrt(vx_*vx_+vy_*vy_+vz_*vz_));
	//ekinafter = 0.5*(vel[q]*vel[q]+vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2])*pop->mass[1] ;

	//double vectrAngle = acos((*vx_point*vx_+*vy_point*vy_+*vz_point*vz_)/(sqrt(vx_*vx_ + vy_*vy_ + vz_*vz_)*sqrt(*vx_point* *vx_point + *vy_point* *vy_point + +*vz_point* *vz_point)));
	//msg(STATUS,"used angle = %f, computed angle = %f",angleChi,vectrAngle);

	double ekinafter = 0.5*(vx_*vx_ + vy_*vy_ + vz_*vz_)*mass[1]; //before back ttransform
	//transform back to lab frame
	vx_ += vxMW;
	vy_ += vyMW;
	vz_ += vzMW;

		// give value back to pointer
	*vx_point=vx_;
	*vy_point=vy_;
	*vz_point=vz_;

	// unittest
	double energydiffexpr = Ekin-Ekin*(cos(angleChi)*cos(angleChi)); //ekin aftr
	//double energydiffexpr = Ekin-Ekin*(0.5+0.5*cos(angleChi) ); //ekin aftr
	double energydiff = Ekin-ekinafter;
	if(abs(energydiffexpr-energydiff)>1e-31){
		msg(WARNING,"too large energy error in collide Ion = %e",(abs(energydiffexpr-energydiff)));
	}
	if(Ekin<newEkin){
		msg(WARNING,"energy increased in Ion collission!");
	}
}

void mccCollideElectronStatic(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng,
	MpiInfo *mpiInfo){

	// uses static CROSS-Sections, collfreq is proportional to v
	//msg(STATUS,"colliding Electrons");

	mccGetPmaxElectronStatic(ini,mccVars,pop,mpiInfo);

	double nt = mccVars->nt;//iniGetDouble(ini,"collisions:numberDensityNeutrals"); //constant for now
	//double NvelThermal = mccVars->NvelThermal;//iniGetDouble(ini,"collisions:thermalVelocityNeutrals");
	double mccSigmaElectronElastic = mccVars->mccSigmaElectronElastic;//iniGetDouble(ini,"collisions:sigmaElectronElastic");

	double maxfreqElectron = mccVars->maxFreqElectron;
	double Pmax = mccVars->pMaxElectron;


	int nDims = pop->nDims;
	double *vel = pop->vel;
	//double *mass = pop->mass;

	//double R1,R2,R3;
	double R = gsl_rng_uniform_pos(rng);
	double Rp = gsl_rng_uniform(rng);

	long int q = 0;
	double* vx;
	double* vy;
	double* vz;
	long int last_i = 0;

	long int errorcounter = 0;

	//double Ekin = 0;
	long int iStart = pop->iStart[0];		// *nDims??
	long int iStop = pop->iStop[0];      // make shure theese are actually electrons!
	long int NparticleColl = (Pmax)*(pop->iStop[0]-pop->iStart[0]); // 1 dim ?, Number of particles too collide
	long int mccStepSize = floor((double)((iStop-iStart))\
	/ (double)(NparticleColl));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "colliding %i of %i electrons\n", (NparticleColl), ((iStop-iStart)));
	}
	//msg(STATUS,"colliding %i of %i electrons", NparticleColl, (iStop-iStart));

	long int mccStop = iStart + mccStepSize*NparticleColl; //enshure non out of bounds

	for(long int i=iStart;i<mccStop;i+=mccStepSize){  //check this statement #segfault long int p=pStart;p<pStart+Pmax*(pStop-pStart);p++)

		R = gsl_rng_uniform_pos(rng); // New random number per particle. maybe print?
		Rp = gsl_rng_uniform_pos(rng); // separate rand num. for prob.
		q = ((i + floor(R*mccStepSize))*nDims);

		vx = &vel[q];
		vy = &vel[q+1];
		vz = &vel[q+2];
		// need n,n+1,n+2 for x,y,j ?

		double MyCollFreq = mccGetMyCollFreqStatic(mccSigmaElectronElastic,  \
			vel[q],vel[q+1],vel[q+2],nt); //prob of coll for particle i

		if (Rp<(MyCollFreq/maxfreqElectron)){
		// Pmax gives the max possible colls. but we will "never" reach this
		// number, so we need to test for each particles probability.
			errorcounter += 1;
			double testvel = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);

			scatterElectron(vx,vy,vz,rng,mccVars,pop);
			double testvel2 = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);
			if(testvel<testvel2){
				msg(STATUS, "vel bef  = %f",testvel);
				msg(STATUS, "vel aftr  = %f",testvel2);
			}
			//msg(STATUS, "vel aftr = %f,%f,%f",vel[q],vel[q+1],vel[q+2]);
			//msg(STATUS, "Ekin = %f, NewEkin = %f",Ekin,newEkin);

		}
		last_i = i;
	}
	// Special handling of last box to let every particle have posibillity to collide

	R = gsl_rng_uniform_pos(rng); // New random number per particle. maybe print?
	Rp = gsl_rng_uniform_pos(rng); // separate rand num. for prob.
	q = ( (last_i) + floor(R*(((iStop)-last_i))) )*nDims; // particle q collides
	// msg(STATUS,"          electrons ");
	// msg(STATUS,"istop = %i, mccstop = %i total pop size = %i",(iStop),mccStop,(iStop-iStart) );
	// msg(STATUS,"mccstep = %i ,last_i = %i, last_q = %i,prev last_q = %i, last box = %i", mccStepSize,last_i,q,last_q,(iStop-last_q));
	// msg(STATUS,"           ");
	// vx = vel[q];
	// vy = vel[q+1];
	// vz = vel[q+2];
	// need n,n+1,n+2 for x,y,j ?
	vx = &vel[q];
	vy = &vel[q+1];
	vz = &vel[q+2];

	double MyCollFreq = mccGetMyCollFreqStatic(mccSigmaElectronElastic,  \
		vel[q],vel[q+1],vel[q+2],nt); //prob of coll for particle i

	if (Rp<(MyCollFreq/maxfreqElectron)){
		errorcounter += 1;

		double testvel = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);

		scatterElectron(vx,vy,vz,rng,mccVars,pop);
		double testvel2 = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);
		if(testvel<testvel2){
			msg(STATUS, "vel bef  = %f",testvel);
			msg(STATUS, "vel aftr  = %f",testvel2);
		}
	}
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "counted  %i Electron collisions on one MPI node \n",  errorcounter);
	}
	//msg(STATUS,"counted  %i Electron collisions on one MPI node", errorcounter);
}

void mccCollideIonStatic(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng, MpiInfo *mpiInfo){

	// uses static CROSS-Sections, collfreq is proportional to v

	mccGetPmaxIonStatic(ini,mccVars,pop,mpiInfo);

	double nt = mccVars->nt;//iniGetDouble(ini,"collisions:numberDensityNeutrals"); //constant for now
	double NvelThermal = mccVars->NvelThermal;//iniGetDouble(ini,"collisions:thermalVelocityNeutrals");
	double mccSigmaCEX= mccVars->mccSigmaCEX;//iniGetDouble(ini,"collisions:sigmaCEX");
	double mccSigmaIonElastic = mccVars->mccSigmaIonElastic;//iniGetDouble(ini,"collisions:sigmaIonElastic");


	double maxfreqIon = mccVars->maxFreqIon;
	double Pmax = mccVars->pMaxIon;

	//int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *vel = pop->vel;
	//double *mass = pop->mass;
	//double velTh = NvelThermal;//(timeStep/stepSize)*NvelThermal; //neutrals
	//msg(STATUS,"neutral thermal vel used in Ion is %f", velTh);
	double Rp, Rq;

	// pointers to pass to scatter function
	double* vx;
	double* vy;
	double* vz;

	long int q = 0;

	long int errorcounter = 0; // count collisions
	long int errorcounter1 = 0; // count cex collisions

	//double vxMW, vyMW, vzMW;

	long int last_i = 0;

	long int iStart = pop->iStart[1];			// *nDims??
	long int iStop = pop->iStop[1];      // make shure theese are actually ions!
	long int NparticleColl = (Pmax)*(pop->iStop[1]-pop->iStart[1]); // 1 dim ?, Number of particles too collide
	long int mccStepSize = (floor((double)((iStop-iStart))\
	/ (double)(NparticleColl)));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "colliding %i of %i ions \n", (NparticleColl), ((iStop-iStart)));
	}
	///msg(STATUS,"colliding %i of %i ions", NparticleColl, (iStop-iStart));
	//msg(STATUS,"Particles per box to pick one from = %i", (mccStepSize));

	long int mccStop = iStart + mccStepSize*NparticleColl; //enshure non out of bounds

	for(long int i=iStart;i<mccStop;i+=mccStepSize){  //check this statement #segfault long int p=pStart;p<pStart+Pmax*(pStop-pStart);p++)

		Rp = gsl_rng_uniform_pos(rng); //decides type of coll.
		Rq = gsl_rng_uniform_pos(rng);
		q = (i + floor(Rq*mccStepSize))*nDims; // particle q collides

		//using particle vx,vy,vz, but should be transformed??
		double MyCollFreq1 = mccGetMyCollFreqStatic(mccSigmaIonElastic,vel[q],
			vel[q+1],vel[q+2], nt);
		double MyCollFreq2 = mccGetMyCollFreqStatic(mccSigmaCEX,vel[q],vel[q+1],
			vel[q+2], nt); //prob of coll for particle i
		//msg(STATUS,"is Rp = %f < MyCollFreq/maxfreqElectron = %f", Rp, ((MyCollFreq1+MyCollFreq2)/ *maxfreqIon));
		if (Rp<( (MyCollFreq1+MyCollFreq2)/maxfreqIon)){ // MyCollFreq1+MyCollFreq2 defines total coll freq.
			//msg(STATUS,"collided");
			if(Rp < MyCollFreq1/maxfreqIon){ // if test slow, but enshures randomness...
				// elastic:
				//msg(STATUS,"elastic");
				errorcounter += 1;

				//point to new velocity
				vx = &vel[q];
				vy = &vel[q+1];
				vz = &vel[q+2];

				//double testvel = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);
				scatterIon(vx,vy,vz,rng,mccVars,pop);
				//double testvel2 = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);

			}else{
				errorcounter += 1;
				errorcounter1 += 1;
				// Charge exchange:
				// flip ion and neutral. Neutral is new ion.
				//msg(STATUS,"ch-ex");
				vel[q] = gsl_ran_gaussian_ziggurat(rng,NvelThermal);; //Vx
				vel[q+1] = gsl_ran_gaussian_ziggurat(rng,NvelThermal);; //Vy
				vel[q+2] = gsl_ran_gaussian_ziggurat(rng,NvelThermal);; //Vz

			}
		}
		last_i = i;
	}

	// Special handling of last box to let every particle have posibillity to collide

	Rp = gsl_rng_uniform_pos(rng); //decides type of coll.
	Rq = gsl_rng_uniform_pos(rng);
	q = ( (last_i) + floor(Rq*(((iStop)-last_i))) )*nDims; // particle q collides

	double MyCollFreq1 = mccGetMyCollFreqStatic(mccSigmaIonElastic,vel[q],vel[q+1],vel[q+2], nt); //duplicate untill real cross sections are implemented
	double MyCollFreq2 = mccGetMyCollFreqStatic(mccSigmaCEX,vel[q],vel[q+1],vel[q+2], nt); //prob of coll for particle i

	if (Rp<( (MyCollFreq1+MyCollFreq2)/maxfreqIon)){ // MyCollFreq1+MyCollFreq2 defines total coll freq.
		if(Rp < MyCollFreq1/maxfreqIon){ // if test slow, but enshures randomness...
			// elastic:
			//msg(STATUS,"elastic");
			errorcounter += 1;

			//point to new velocity
			vx = &vel[q];
			vy = &vel[q+1];
			vz = &vel[q+2];

			scatterIon(vx,vy,vz,rng,mccVars,pop);

		}else{
			// Charge exchange:
			// flip ion and neutral. Neutral is new ion.
			errorcounter += 1;
			errorcounter1 += 1;
			//msg(STATUS,"ch-ex");
			vel[q] = gsl_ran_gaussian_ziggurat(rng,NvelThermal);  //Vx
			vel[q+1] = gsl_ran_gaussian_ziggurat(rng,NvelThermal);  //Vy
			vel[q+2] = gsl_ran_gaussian_ziggurat(rng,NvelThermal); //Vz

		}
	}
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "counted  %i ION collisions on one MPI node, %i as CEX \n", errorcounter,errorcounter1);
	}
	//msg(STATUS,"counted  %i ION collisions on one MPI node, %i as CEX", errorcounter,errorcounter1);
}

void mccCollideElectronFunctional(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng,
	MpiInfo *mpiInfo){

	// uses static CROSS-Sections, collfreq is proportional to v
	//msg(STATUS,"colliding Electrons");

	mccGetPmaxElectronFunctional(ini,mccVars,pop,mpiInfo);

	double nt = mccVars->nt; //constant for now
	double electron_a = mccVars->electron_a;
	double electron_b = mccVars->electron_b;
	double maxfreqElectron = mccVars->maxFreqElectron;
	double Pmax = mccVars->pMaxElectron;
	int nDims = pop->nDims;
	double *vel = pop->vel;

	double R = gsl_rng_uniform_pos(rng);
	double Rp = gsl_rng_uniform(rng);

	long int q = 0;
	double* vx;
	double* vy;
	double* vz;

	long int last_i = 0;
	long int errorcounter = 0;

	long int iStart = pop->iStart[0];		// *nDims??
	long int iStop = pop->iStop[0];      // make shure theese are actually electrons!
	long int NparticleColl = (Pmax)*(pop->iStop[0]-pop->iStart[0]); // 1 dim ?, Number of particles too collide
	long int mccStepSize = floor((double)((iStop-iStart))\
	/ (double)(NparticleColl));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "colliding %i of %i electrons\n", (NparticleColl), ((iStop-iStart)));
	}
	//msg(STATUS,"colliding %i of %i electrons", NparticleColl, (iStop-iStart));

	long int mccStop = iStart + mccStepSize*NparticleColl; //enshure non out of bounds

	for(long int i=iStart;i<mccStop;i+=mccStepSize){  //check this statement #segfault long int p=pStart;p<pStart+Pmax*(pStop-pStart);p++)

		R = gsl_rng_uniform_pos(rng); // New random number per particle. maybe print?
		Rp = gsl_rng_uniform_pos(rng); // separate rand num. for prob.
		q = ((i + floor(R*mccStepSize))*nDims);

		vx = &vel[q];
		vy = &vel[q+1];
		vz = &vel[q+2];
		// need n,n+1,n+2 for x,y,j ?


		double MyCollFreq = mccGetMyCollFreqFunctional(mccSigmaElectronElasticFunctional,  \
			vel[q],vel[q+1],vel[q+2],nt,electron_a,electron_b); //prob of coll for particle i

		if (Rp<(MyCollFreq/maxfreqElectron)){
			errorcounter += 1;

			double testvel = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);

			scatterElectron(vx,vy,vz,rng,mccVars,pop);
			double testvel2 = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);
			if(testvel<testvel2){
				msg(STATUS, "vel bef  = %f",testvel);
				msg(STATUS, "vel aftr  = %f",testvel2);
			}
		}
		if(mpiInfo->mpiRank==0){
			fMsg(ini, "collision", "counted  %i Electron collisions on one MPI node \n",  errorcounter);
		}
		last_i = i;
	}
	// Special handling of last box to let every particle have posibillity to collide

	R = gsl_rng_uniform_pos(rng); // New random number per particle. maybe print?
	Rp = gsl_rng_uniform_pos(rng); // separate rand num. for prob.
	q = ( (last_i) + floor(R*(((iStop)-last_i))) )*nDims; // particle q collides

	vx = &vel[q];
	vy = &vel[q+1];
	vz = &vel[q+2];
	// need n,n+1,n+2 for x,y,j ?

	double MyCollFreq = mccGetMyCollFreqFunctional(mccSigmaElectronElasticFunctional,  \
		vel[q],vel[q+1],vel[q+2],nt,electron_a,electron_b); //prob of coll for particle i


		if (Rp<(MyCollFreq/maxfreqElectron)){
			errorcounter += 1;

			double testvel = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);

			scatterElectron(vx,vy,vz,rng,mccVars,pop);
			double testvel2 = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);
			if(testvel<testvel2){
				msg(STATUS, "vel bef  = %f",testvel);
				msg(STATUS, "vel aftr  = %f",testvel2);
			}
		}
		if(mpiInfo->mpiRank==0){
			fMsg(ini, "collision", "counted  %i Electron collisions on one MPI node \n",  errorcounter);
		}
		//msg(STATUS,"counted  %i Electron collisions on one MPI node", errorcounter);
}

void mccCollideIonFunctional(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng,MpiInfo *mpiInfo){

	// uses static CROSS-Sections, collfreq is proportional to v



	double nt = mccVars->nt;
	double NvelThermal = mccVars->NvelThermal;
	double CEX_a = mccVars->CEX_a;
	double CEX_b = mccVars->CEX_b;
	double ion_elastic_a =mccVars->ion_elastic_a;
	double ion_elastic_b =mccVars->ion_elastic_b;

	mccGetPmaxIonFunctional(ini,mccVars,pop,mpiInfo);
	double maxfreqIon = mccVars->maxFreqIon;
	double Pmax = mccVars->pMaxIon;
	//int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *vel = pop->vel;

	double Rp, Rq;
	long int q = 0;
	long int last_i = 0;

	long int errorcounter = 0; // count collisions
	long int errorcounter1 = 0; // count cex collisions

	// pointers to pass to scatter function
	double* vx;
	double* vy;
	double* vz;

	long int iStart = pop->iStart[1];			// *nDims??
	long int iStop = pop->iStop[1];      // make shure theese are actually ions!
	long int NparticleColl = (Pmax)*(pop->iStop[1]-pop->iStart[1]); // 1 dim ?, Number of particles too collide
	long int mccStepSize = floor((double)((iStop-iStart))\
	/ (double)(NparticleColl));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "colliding %i of %i ions \n", (NparticleColl), ((iStop-iStart)));
	}
	///msg(STATUS,"colliding %i of %i ions", NparticleColl, (iStop-iStart));
	//msg(STATUS,"Particles per box to pick one from = %i", (mccStepSize));

	long int mccStop = iStart + mccStepSize*NparticleColl; //enshure non out of bounds
	//check energies
	//double ekin = 0; //test value

	//double ekinafter =0;

	for(long int i=iStart;i<mccStop;i+=mccStepSize){  //check this statement #segfault long int p=pStart;p<pStart+Pmax*(pStop-pStart);p++)

		Rp = gsl_rng_uniform_pos(rng); //decides type of coll.
		//R = gsl_rng_uniform_pos(rng); // New random number per particle. maybe print?
		Rq = gsl_rng_uniform_pos(rng);
		q = (i + floor(Rq*mccStepSize))*nDims; // particle q collides
		//msg(STATUS,"R=%g, p = %i, q = %i", R,p,q );



		double MyCollFreq1 = mccGetMyCollFreqFunctional(
			mccSigmaIonElasticFunctional,vel[q],
				vel[q+1],vel[q+2],nt,ion_elastic_a,ion_elastic_b); //prob of coll for particle i
		double MyCollFreq2 = mccGetMyCollFreqFunctional(mccSigmaCEXFunctional,
			vel[q],vel[q+1],vel[q+2],nt,CEX_a,CEX_b); //duplicate untill real cross sections are implemented

		//msg(STATUS,"is Rp = %f < MyCollFreq/maxfreqElectron = %f", Rp, ((MyCollFreq1+MyCollFreq2)/ *maxfreqIon));
		if (Rp<( (MyCollFreq1+MyCollFreq2)/maxfreqIon)){ // MyCollFreq1+MyCollFreq2 defines total coll freq.
			if(Rp < MyCollFreq1/maxfreqIon){ // if test slow, but enshures randomness...
				// elastic:
				//msg(STATUS,"elastic");
				errorcounter += 1;
				//point to new velocity
				vx = &vel[q];
				vy = &vel[q+1];
				vz = &vel[q+2];

				double testvel = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);

				scatterIon(vx,vy,vz,rng,mccVars,pop);

				double testvel2 = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);
				if(testvel<testvel2){
					msg(STATUS, "vel bef  = %f",testvel);
					msg(STATUS, "vel aftr  = %f",testvel2);
				}else{
					msg(STATUS, "NO ENERGY INCREASE");
				}
			}else{
				errorcounter += 1;
				errorcounter1 += 1;
				// Charge exchange:
				// flip ion and neutral. Neutral is new ion.
				//msg(STATUS,"ch-ex");

				vel[q] = gsl_ran_gaussian_ziggurat(rng,NvelThermal);  //Vx
				vel[q+1] = gsl_ran_gaussian_ziggurat(rng,NvelThermal);  //Vy
				vel[q+2] = gsl_ran_gaussian_ziggurat(rng,NvelThermal);  //Vz

			}
		}
		last_i = i;
	}
	// Special handling of last box to let every particle have posibillity to collide

	Rp = gsl_rng_uniform_pos(rng); //decides type of coll.
	//R = gsl_rng_uniform_pos(rng); // New random number per particle. maybe print?
	Rq = gsl_rng_uniform_pos(rng);
	q = ( (last_i) + floor(Rq*(((iStop)-last_i))) )*nDims; // particle q collides


	double MyCollFreq1 = mccGetMyCollFreqFunctional(
		mccSigmaIonElasticFunctional,vel[q],
			vel[q+1],vel[q+2],nt,ion_elastic_a,ion_elastic_b); //prob of coll for particle i
	double MyCollFreq2 = mccGetMyCollFreqFunctional(mccSigmaCEXFunctional,
			vel[q],vel[q+1],vel[q+2],nt,CEX_a,CEX_b); //duplicate untill real cross sections are implemented

	if (Rp<( (MyCollFreq1+MyCollFreq2)/maxfreqIon)){ // MyCollFreq1+MyCollFreq2 defines total coll freq.
		if(Rp < MyCollFreq1/maxfreqIon){ // if test slow, but enshures randomness...
			// elastic:
			//msg(STATUS,"elastic");
			errorcounter += 1;

			//point to new velocity
			vx = &vel[q];
			vy = &vel[q+1];
			vz = &vel[q+2];

			double testvel = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);

			scatterIon(vx,vy,vz,rng,mccVars,pop);

			double testvel2 = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);
			if(testvel<testvel2){
				msg(STATUS, "vel bef  = %f",testvel);
				msg(STATUS, "vel aftr  = %f",testvel2);
			}else{
				msg(STATUS, "NO ENERGY INCREASE");
			}

		}else{
			// Charge exchange:
			// flip ion and neutral. Neutral is new ion.
			errorcounter += 1;
			errorcounter1 += 1;
			//msg(STATUS,"ch-ex");
			vel[q] = gsl_ran_gaussian_ziggurat(rng,NvelThermal);  //Vx
			vel[q+1] = gsl_ran_gaussian_ziggurat(rng,NvelThermal);  //Vy
			vel[q+2] = gsl_ran_gaussian_ziggurat(rng,NvelThermal);  //Vz

		}
	}
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "counted  %i ION collisions on one MPI node, %i as CEX \n", errorcounter,errorcounter1);
	}
	//msg(STATUS,"counted  %i ION collisions on one MPI node, %i as CEX", errorcounter,errorcounter1);
}

void mccCollideElectronConstantFrq(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng,
	MpiInfo *mpiInfo){

	//  collfreq is constant
	//msg(STATUS,"colliding Electrons");

	mccGetPmaxElectronConstantFrq(ini,mccVars,pop,mpiInfo);

	int nDims = pop->nDims;
	double *vel = pop->vel;
	double Pmax = mccVars->pMaxElectron;

	long int q = 0;
	double* vx;
	double* vy;
	double* vz;
	double R = gsl_rng_uniform_pos(rng);
	//double vxMW, vyMW, vzMW;
	long int last_i = 0;

	//double ekinafter = 0;

	long int errorcounter = 0;

	long int iStart = pop->iStart[0];		// *nDims??
	long int iStop = pop->iStop[0];      // make shure theese are actually electrons!
	long int NparticleColl = (Pmax)*(pop->iStop[0]-pop->iStart[0]); // 1 dim ?, Number of particles too collide
	long int mccStepSize = floor((double)((iStop-iStart))\
	/ (double)(NparticleColl));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "colliding %i of %i electrons\n", (NparticleColl), ((iStop-iStart)));
	}
	//msg(STATUS,"colliding %i of %i electrons", NparticleColl, (iStop-iStart));

	long int mccStop = iStart + mccStepSize*NparticleColl; //enshure non out of bounds

	for(long int i=iStart;i<mccStop;i+=mccStepSize){  //check this statement #segfault long int p=pStart;p<pStart+Pmax*(pStop-pStart);p++)

		R = gsl_rng_uniform_pos(rng); // New random number per particle. maybe print?
		q = ((i + floor(R*mccStepSize))*nDims);

		vx = &vel[q];
		vy = &vel[q+1];
		vz = &vel[q+2];
		// need n,n+1,n+2 for x,y,j ?

		// Pmax gives the max possible colls. but we will "never" reach this
		// number, so we need to test for each particles probability.
		errorcounter += 1;
		double testvel = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);

		scatterElectron(vx,vy,vz,rng,mccVars,pop);

		double testvel2 = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);
		if(testvel<testvel2){
			msg(STATUS, "vel bef  = %f",testvel);
			msg(STATUS, "vel aftr  = %f",testvel2);
		}
		if(testvel2<1e-32){
			msg(STATUS, "vel bef  = %f",testvel);
			msg(STATUS, "vel aftr  = %f",testvel2);
		}
		//msg(STATUS, "vel aftr = %f,%f,%f",vel[q],vel[q+1],vel[q+2]);
		//msg(STATUS, "Ekin = %f, NewEkin = %f",Ekin,newEkin);
		last_i = i;
	}

	// Special handling of last box to let every particle have posibillity to collide

	R = gsl_rng_uniform_pos(rng); // New random number per particle. maybe print?
	q = ( (last_i) + floor(R*(((iStop)-last_i))) )*nDims; // particle q collides


	vx = &vel[q];
	vy = &vel[q+1];
	vz = &vel[q+2];
	// need n,n+1,n+2 for x,y,j ?


	// Pmax gives the max possible colls. but we will "never" reach this
	// number, so we need to test for each particles probability.
	errorcounter += 1;
	double testvel = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);

	scatterElectron(vx,vy,vz,rng,mccVars,pop);
	double testvel2 = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);
	if(testvel<testvel2){
		msg(STATUS, "vel bef  = %f",testvel);
		msg(STATUS, "vel aftr  = %f",testvel2);
	}
	if(testvel2<1e-32){
		msg(STATUS, "vel bef  = %f",testvel);
		msg(STATUS, "vel aftr  = %f",testvel2);
	}
	//msg(STATUS, "vel aftr = %f,%f,%f",vel[q],vel[q+1],vel[q+2]);
	//msg(STATUS, "Ekin = %f, NewEkin = %f",Ekin,newEkin);


	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "counted  %i Electron collisions on one MPI node \n",  errorcounter);
	}
	//msg(STATUS,"counted  %i Electron collisions on one MPI node", errorcounter);
}

void mccCollideIonConstantFrq(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng, MpiInfo *mpiInfo){

	// uses constant collfreq fixed number of colls per dt

	//double nt = iniGetDouble(ini,"collisions:numberDensityNeutrals"); //constant for now
	//double mccTimestep = iniGetDouble(ini,"time:timeStep");
	double NvelThermal = mccVars->NvelThermal;

	double collFrqIonElastic = mccVars->collFrqIonElastic;
	double collFrqIonCEX = mccVars->collFrqCex;

	mccGetPmaxIonConstantFrq(ini,mccVars,pop,mpiInfo);
	double Pmax = mccVars->pMaxIon;


	//int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *vel = pop->vel;
	//double velTh = NvelThermal;//(timeStep/stepSize)*NvelThermal; //neutrals
	//msg(STATUS,"neutral thermal vel used in Ion is %f", velTh);
	double Rp, Rq;
	long int q = 0;
	long int last_i = 0;

	long int errorcounter = 0; // count collisions
	long int errorcounter1 = 0; // count cex collisions

	double* vx;
	double* vy;
	double* vz;

	double MyCollFreq1 = collFrqIonElastic;
	double MyCollFreq2 = collFrqIonCEX;
	double maxfreqIon = MyCollFreq1+MyCollFreq2;

	long int iStart = pop->iStart[1];			// *nDims??
	long int iStop = pop->iStop[1];      // make shure theese are actually ions!
	long int NparticleColl = (Pmax)*(pop->iStop[1]-pop->iStart[1]); // 1 dim ?, Number of particles too collide
	long int mccStepSize = floor((double)((iStop-iStart))\
	/ (double)(NparticleColl));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "colliding %i of %i ions \n", (NparticleColl), ((iStop-iStart)));
	}
	///msg(STATUS,"colliding %i of %i ions", NparticleColl, (iStop-iStart));
	//msg(STATUS,"Particles per box to pick one from = %i", (mccStepSize));

	long int mccStop = iStart + mccStepSize*NparticleColl; //enshure non out of bounds

	//check energies
	//double ekin = 0; //test value
	//double ekinafter =0;

	for(long int i=iStart;i<mccStop;i+=mccStepSize){  //check this statement #segfault long int p=pStart;p<pStart+Pmax*(pStop-pStart);p++)

		Rp = gsl_rng_uniform_pos(rng); //decides type of coll.
		Rq = gsl_rng_uniform_pos(rng);
		q = (i + floor(Rq*mccStepSize))*nDims; // particle q collides
		//msg(STATUS,"mass = %f, i = %i, q = %i",pop->mass[1],i*nDims,q );

		//ekin = 0.5*(vxTran*vxTran+vyTran*vyTran+vzTran*vzTran)*pop->mass[1] ;
		//msg(STATUS, "vel befor = %f,%f,%f",vel[q],vel[q+1],vel[q+2]);

		//msg(STATUS,"is Rp = %f < MyCollFreq/maxfreqElectron = %f", Rp, ((MyCollFreq1+MyCollFreq2)/ *maxfreqIon));
		if(Rp < MyCollFreq1/ maxfreqIon){ // if test slow, but enshures randomness...
			// elastic:
			//msg(STATUS,"elastic");
			errorcounter += 1;

			//point to new velocity
			vx = &vel[q];
			vy = &vel[q+1];
			vz = &vel[q+2];

			double testvel = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);

			scatterIon(vx,vy,vz,rng,mccVars,pop);

			double testvel2 = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);
			// if(testvel<testvel2){
			// 	msg(STATUS, "vel bef  = %f",testvel);
			// 	msg(STATUS, "vel aftr  = %f",testvel2);
			// }else{
			// 	msg(STATUS, "NO ENERGY INCREASE");
			// }
			//msg(STATUS,"velocity changed by a factor %f",testvel2/testvel);
		}else{
			errorcounter += 1;
			errorcounter1 += 1;
			// Charge exchange:
			// flip ion and neutral. Neutral is new ion.
			//msg(STATUS,"ch-ex");
			vel[q] = gsl_ran_gaussian_ziggurat(rng,NvelThermal); //Vx
			vel[q+1] = gsl_ran_gaussian_ziggurat(rng,NvelThermal); //Vy
			vel[q+2] = gsl_ran_gaussian_ziggurat(rng,NvelThermal); //Vz
		}
		last_i = i;
	}



	// Special handling of last box to let every particle have posibillity to collide

	Rp = gsl_rng_uniform_pos(rng); //decides type of coll.
	Rq = gsl_rng_uniform_pos(rng);
	q = ( (last_i) + floor(Rq*(((iStop)-last_i))) )*nDims; // particle q collides
	// msg(STATUS,"istop = %i, total pop size = %i",(iStop),(iStop-iStart) );
	// msg(STATUS,"mccstep = %i ,last_q = %i, prev last_i = %i, last box = %i", mccStepSize,q,last_i,(iStop-last_i));

	if(Rp < MyCollFreq1/ maxfreqIon){ // if test slow, but enshures randomness...
		// elastic:
		//msg(STATUS,"elastic");
		errorcounter += 1;

		//point to new velocity
		vx = &vel[q];
		vy = &vel[q+1];
		vz = &vel[q+2];

		//double testvel = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);

		scatterIon(vx,vy,vz,rng,mccVars,pop);

		//double testvel2 = sqrt(vel[q]*vel[q]+ vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2]);
		// if(testvel<testvel2){
		// 	msg(STATUS, "vel bef  = %f",testvel);
		// 	msg(STATUS, "vel aftr  = %f",testvel2);
		// }else{
		// 	msg(STATUS, "NO ENERGY INCREASE");
		// }
	}else{
		// Charge exchange:
		// flip ion and neutral. Neutral is new ion.
		errorcounter += 1;
		errorcounter1 += 1;
		//msg(STATUS,"ch-ex");
		vel[q] = gsl_ran_gaussian_ziggurat(rng,NvelThermal); //Vx
		vel[q+1] = gsl_ran_gaussian_ziggurat(rng,NvelThermal); //Vy
		vel[q+2] = gsl_ran_gaussian_ziggurat(rng,NvelThermal); //Vz


	}

	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "counted  %i ION collisions on one MPI node, %i as CEX \n", errorcounter,errorcounter1);
	}
	//msg(STATUS,"counted  %i ION collisions on one MPI node, %i as CEX", errorcounter,errorcounter1);
}

/*************************************************
*		Method Handler
************************************************/

void mccCollideConstantCrossect(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng, MpiInfo *mpiInfo){

	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "\n Computing time-step \n");
	}
	mccCollideIonStatic(ini, pop, mccVars, rng, mpiInfo);
	mccCollideElectronStatic(ini, pop, mccVars, rng,mpiInfo);
}

funPtr constCrossect_set(dictionary *ini){

	mccSanity(ini,"mccMode",2);
	funPtr collissions;
	//point to function that calls both collision functions.. or more
	collissions = &mccCollideConstantCrossect;

	return collissions;
}

void mccCollideConstantFreq(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng, MpiInfo *mpiInfo){

	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "\n Computing time-step \n");
	}

	mccCollideIonConstantFrq(ini, pop,mccVars,rng,mpiInfo);
	mccCollideElectronConstantFrq(ini, pop,mccVars, rng,mpiInfo);
}

funPtr constFreq_set(dictionary *ini){

	mccSanity(ini,"mccMode",2);
	funPtr collissions;
	//point to function that calls both collision functions.. or more
	collissions = &mccCollideConstantFreq;

	return collissions;
}

void mccCollideFunctional(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng, MpiInfo *mpiInfo){

	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "\n Computing time-step \n");
	}
	mccCollideIonFunctional(ini,pop,mccVars,rng,mpiInfo);
	mccCollideElectronFunctional(ini, pop,mccVars, rng,mpiInfo);
}

funPtr functionalCrossect_set(dictionary *ini){

	mccSanity(ini,"mccMode",2);
	funPtr collissions;
	//point to function that calls both collision functions.. or more
	collissions = &mccCollideFunctional;

	return collissions;
}

void collissionsOff(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng, MpiInfo *mpiInfo){
	NULL;
}

funPtr collissionsOff_set(dictionary *ini){

	mccSanity(ini,"mccMode",2);
	funPtr collissions;
	//point to function that calls both collision functions.. or more
	collissions = &collissionsOff;

	return collissions;
}



/*************************************************
*		RUNS
************************************************/

funPtr mccMode_set(dictionary *ini){
	//test sanity here!
	mccSanity(ini,"mccMode",2);
	return mccMode;
}



void mccMode(dictionary *ini){

	msg(STATUS, "start mcc Test Mode");

	/*
	* SELECT METHODS
	*/
	void (*acc)()   = select(ini,"methods:acc",	puAcc3D1_set,
	puAcc3D1KE_set,
	puAccND1_set,
	puAccND1KE_set,
	puAccND0_set,
	puAccND0KE_set,
	puBoris3D1_set,
	puBoris3D1KE_set,
	puBoris3D1KETEST_set);

	void (*distr)() = select(ini,"methods:distr",	puDistr3D1split_set,
													puDistr3D1_set,
													puDistrND1_set,
													puDistrND0_set);

	void (*extractEmigrants)() = select(ini,"methods:migrate",	puExtractEmigrants3D_set,
	puExtractEmigrantsND_set);

	void (*solverInterface)()	= select(ini,	"methods:poisson",
												mgSolver_set
												//sSolver_set
												);

	void (*collide)() = select(ini,	"methods:mcc",
									collissionsOff_set,
									constCrossect_set,
									constFreq_set,
									functionalCrossect_set);

	//
	void (*solve)() = NULL;
	void *(*solverAlloc)() = NULL;
	void (*solverFree)() = NULL;
	solverInterface(&solve, &solverAlloc, &solverFree);


	/*
	* INITIALIZE PINC VARIABLES
	*
	*done in "main" for "regular mode"
	*/
	Units *units=uAlloc(ini);
	uNormalize(ini, units);
	// normalize mcc input variables
	//must be done after uNormalize, and before defining mcc vars
	mccNormalize(ini,units);

	MpiInfo *mpiInfo = gAllocMpi(ini);
	Population *pop = pAlloc(ini);
	Grid *phi = gAlloc(ini, SCALAR);
	Grid *E   = gAlloc(ini, VECTOR);
	Grid *rho = gAlloc(ini, SCALAR);
	Grid *rho_e = gAlloc(ini, SCALAR);
	Grid *rho_i = gAlloc(ini, SCALAR);
	void *solver = solverAlloc(ini, rho, phi);


	/*
	* mcc specific variables
	*/

	MccVars *mccVars=mccAlloc(ini,units);
	//double PmaxElectron = 0;
	//double maxfreqElectron = 0;
	//double PmaxIon = 0;
	//double maxfreqIon = 0;

	// using Boris algo
	int nSpecies = pop->nSpecies;
	double *S = (double*)malloc((3)*(nSpecies)*sizeof(double));
	double *T = (double*)malloc((3)*(nSpecies)*sizeof(double));

	// Creating a neighbourhood in the rho to handle migrants
	gCreateNeighborhood(ini, mpiInfo, rho);

	// Setting Boundary slices
	gSetBndSlices(phi, mpiInfo);

	//Set mgSolve
	//MgAlgo mgAlgo = getMgAlgo(ini);

	// Random number seeds
	gsl_rng *rngSync = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,mpiInfo->mpiRank+1); // Seed needs to be >=1

	// /*
	// * PREPARE FILES FOR WRITING
	// */
	// int rank = phi->rank;
	// double *denorm = malloc((rank-1)*sizeof(*denorm));
	// double *dimen = malloc((rank-1)*sizeof(*dimen));
	//
	// for(int d = 1; d < rank;d++) denorm[d-1] = 1.;
	// for(int d = 1; d < rank;d++) dimen[d-1] = 1.;

	/*
	 * PREPARE FILES FOR WRITING
	 */
	pOpenH5(ini, pop, units, "pop");
	gOpenH5(ini, rho, mpiInfo, units, units->chargeDensity, "rho");
	gOpenH5(ini, rho_e, mpiInfo, units, units->chargeDensity, "rho_e");
	gOpenH5(ini, rho_i, mpiInfo, units, units->chargeDensity, "rho_i");
	gOpenH5(ini, phi, mpiInfo, units, units->potential, "phi");
	gOpenH5(ini, E,   mpiInfo, units, units->eField, "E");
  // oOpenH5(ini, obj, mpiInfo, units, 1, "test");
  // oReadH5(obj, mpiInfo);

	hid_t history = xyOpenH5(ini,"history");
	pCreateEnergyDatasets(history,pop);
	hid_t temperature = xyOpenH5(ini,"temperature");
	pCreateTemperatureDatasets(temperature,pop);
	hid_t probe = arrOpenH5(ini,"probe");
	xyzProbeCreateDatasets(probe,phi,mpiInfo);

	// Add more time series to history if you want
	// xyCreateDataset(history,"/group/group/dataset");

	// free(denorm);
	// free(dimen);

	/*
	* INITIAL CONDITIONS
	*/

	//msg(STATUS, "constants = %f, %f", velThermal[0], velThermal[1]);
	// Initalize particles
	//pPosLattice(ini, pop, mpiInfo);
	pPosUniform(ini, pop, mpiInfo, rngSync);
	//pVelZero(pop);
	//double *thermalVelocity = iniGetDoubleArr(ini,"population:thermalVelocity",nSpecies);

	//pVelConstant(ini, pop,(thermalVelocity[0]),(thermalVelocity[1]) ); //constant values for vel.
	pVelMaxwell(ini, pop, rng);
	//double maxVel = iniGetDouble(ini,"population:maxVel");

	// Perturb particles
	//pPosPerturb(ini, pop, mpiInfo);

	// Migrate those out-of-bounds due to perturbation
	extractEmigrants(pop, mpiInfo);
	puMigrate(pop, mpiInfo, rho);

	/*
	* compute initial half-step
	*/

	// Get initial charge density
	distr(pop, rho); //two species
	gHaloOp(addSlice, rho, mpiInfo, FROMHALO);
	gHaloOp(addSlice, rho_e, mpiInfo, FROMHALO);
	gHaloOp(addSlice, rho_i, mpiInfo, FROMHALO);

	// Get initial E-field
	//solve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);
	solve(solver, rho, phi, mpiInfo);
	gFinDiff1st(phi, E);
	gHaloOp(setSlice, E, mpiInfo, TOHALO);
	gMul(E, -1.);

	// add External E
	//gZero(E);
	puAddEext(ini, pop, E);

	// Advance velocities half a step
	gMul(E, 0.5);
	puGet3DRotationParameters(ini, T, S, 0.5);
	acc(pop, E, T, S);
	gMul(E, 2.0);
	puGet3DRotationParameters(ini, T, S, 1.0);

	//Write initial h5 files
	//gWriteH5(E, mpiInfo, 0.0);
	//gWriteH5(rho, mpiInfo, 0.0);
	//gWriteH5(rho_e, mpiInfo, 0.0);
	//gWriteH5(rho_i, mpiInfo, 0.0);
	//gWriteH5(phi, mpiInfo, 0.0);
	pWriteH5(pop, mpiInfo, 0.0, 0.5,1);
	pWriteTemperature(temperature,pop,0.0,units,ini);
	pWriteEnergy(history,pop,0.0,units);
	xyzWriteProbe(probe, phi,0.0,mpiInfo);

	//msg(STATUS, "Pmax for Electrons is %f",PmaxElectron);
	//msg(STATUS, "Pmax for Ions is %f",PmaxIon);

	/*
	* TIME LOOP
	*/

	Timer *t = tAlloc(mpiInfo->mpiRank);

	// n should start at 1 since that's the timestep we have after the first
	// iteration (i.e. when storing H5-files).
	int nTimeSteps = iniGetInt(ini,"time:nTimeSteps");
	for(int n = 1; n <= nTimeSteps; n++){

		msg(STATUS,"Computing time-step %i of %i",n,nTimeSteps);
		//MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary

		// Check that no particle moves beyond a cell (mostly for debugging)
		//pVelAssertMax(pop,maxVel);

		tStart(t);


		// Move particles
		//msg(STATUS, "moving particles");
		puMove(pop);
		// oRayTrace(pop, obj);


		// Migrate particles (periodic boundaries)
		//msg(STATUS, "extracting emigrants");
		extractEmigrants(pop, mpiInfo);
		//msg(STATUS, "Migrating particles");
		puMigrate(pop, mpiInfo, rho);

		/*
		*   Collisions
		*/

		//mccCollideElectron_debug(ini,pop, &PmaxElectron, &maxfreqElectron, rng, mccTimestep, nt,DebyeLength); //race conditions?????
		//mccCollideIon_debug(ini, pop, &PmaxIon, &maxfreqIon, rng, mccTimestep, nt,DebyeLength);

		// mccGetPmaxElectronConstantFrq(ini, &PmaxElectron,
		// 	&collFrqElectronElastic,pop,mpiInfo);
		// mccGetPmaxIonConstantFrq(ini, &PmaxIon, &collFrqIonElastic,
		// 	&collFrqCEX, pop,mpiInfo);
		//
		// mccCollideIonConstantFrq(ini, pop,&PmaxIon, rng,NvelThermal,
		// 	&collFrqIonElastic,&collFrqCEX, mpiInfo);
		// mccCollideElectronConstantFrq(ini, pop,&PmaxElectron,&maxfreqElectron,
		// 	rng,NvelThermal,nt,mpiInfo,(1.0));

		//mccGetPmaxElectronStatic(ini,mccVars,pop,mpiInfo);
		//mccGetPmaxIonStatic(ini,mccVars, pop,mpiInfo);
		//
		// mccCollideIonStatic(ini, pop, mccVars, rng, mpiInfo);
		// mccCollideElectronStatic(ini, pop, mccVars, rng,mpiInfo,(1.0/10.0));

		collide(ini, pop, mccVars, rng,mpiInfo);

		// mccGetPmaxElectronFunctional(ini,
		// 		mccTimestep,nt,&PmaxElectron, &maxfreqElectron,pop,
		// 		mpiInfo,electron_a, electron_b);
		// mccGetPmaxIonFunctional(ini,
		// 		mccTimestep, nt, &PmaxIon, &maxfreqIon, pop,
		// 		mpiInfo,CEX_a,CEX_b, ion_elastic_a, ion_elastic_b);
		//
		// mccCollideElectronFunctional(ini, pop,
		// 	&PmaxElectron, &maxfreqElectron, rng,NvelThermal, nt, mpiInfo,
		// 	electron_a, electron_b,1.0);
		// mccCollideIonFunctional(ini, pop, &PmaxIon, &maxfreqIon,
		// 	rng, nt, NvelThermal, mpiInfo,CEX_a,CEX_b, ion_elastic_a,
		// 	ion_elastic_b);

		/*
		*   Collisions
		*/

		//MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary

		//MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary
		// Check that no particle resides out-of-bounds (just for debugging)
		//msg(STATUS, "checking particles out of bounds");
		//pPosAssertInLocalFrame(pop, rho);

		//MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary
		//-----------------------------
		// Compute charge density
		//msg(STATUS, "computing charge density");
		distr(pop, rho); //two species

		//MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary

		//msg(STATUS, "MPI exchanging density");
		gHaloOp(addSlice, rho, mpiInfo, FROMHALO);
		gHaloOp(addSlice, rho_e, mpiInfo, FROMHALO);
		gHaloOp(addSlice, rho_i, mpiInfo, FROMHALO);

		//---------------------------
		//msg(STATUS, "gAssertNeutralGrid rho");
		//gAssertNeutralGrid(rho, mpiInfo);

		// Compute electric potential phi
		//msg(STATUS, "Compute electric potential phi");
		//solve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);
		solve(solver, rho, phi, mpiInfo);
		//msg(STATUS, "gAssertNeutralGrid phi");
		gAssertNeutralGrid(phi, mpiInfo);

		// Compute E-field
		//msg(STATUS, "compute E field");
		gFinDiff1st(phi, E);
		//msg(STATUS, "MPI exchanging E");
		gHaloOp(setSlice, E, mpiInfo, TOHALO);
		//msg(STATUS, "gMul( E, -1)");
		gMul(E, -1.);

		//msg(STATUS, "gAssertNeutralGrid E");
		gAssertNeutralGrid(E, mpiInfo);
		// Apply external E
		// gAddTo(Ext);
		//gZero(E); ////temporary test
		puAddEext(ini, pop, E);

		// Accelerate particle and compute kinetic energy for step n
		//msg(STATUS, "Accelerate particles");

		acc(pop, E, T, S);
		tStop(t);

		// Sum energy for all species
		//msg(STATUS, "sum kinetic energy");
		pSumKinEnergy(pop); // kin_energy per species to pop
		//msg(STATUS, "compute pot energy");
		// Compute potential energy for step n
		gPotEnergy(rho,phi,pop);

		if( n%500 == 0){
			// Example of writing another dataset to history.xy.h5
			// xyWrite(history,"/group/group/dataset",(double)n,value,MPI_SUM);
			//Write h5 files
			//gWriteH5(E, mpiInfo, (double) n);
			//gWriteH5(rho, mpiInfo, (double) n);
			//gWriteH5(rho_e, mpiInfo, (double) n);
			//gWriteH5(rho_i, mpiInfo, (double) n);
			//gWriteH5(phi, mpiInfo, (double) n);
			pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5,1);
		}
		// if( n%10000 == 0 && n<20000){
		// 	pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5,1);
		// }
		// if(n==1){
		// 	pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5,1); //0.0001
		// 	//gWriteH5(E, mpiInfo, (double) n);
		// 	//gWriteH5(rho, mpiInfo, (double) n);
		// 	//gWriteH5(rho_e, mpiInfo, (double) n);
		// 	//gWriteH5(rho_i, mpiInfo, (double) n);
		// 	//gWriteH5(phi, mpiInfo, (double) n);
		// 	//gWriteH5(phi, mpiInfo, (double) n);
		// }

		//if(n == 30000){
			//pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5,1); //0.0001
		//}


		//if( n< 40000 && n%500 == 0){
		//}

		// if(n>nTimeSteps-1){
		// 	msg(STATUS, "writing over a given timestep to file");
		// 	//pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5,0.001);
		// 	//gWriteH5(rho, mpiInfo, (double) n);
		// 	//gWriteH5(rho_e, mpiInfo, (double) n);
		// 	//gWriteH5(rho_i, mpiInfo, (double) n);
		// 	//gWriteH5(phi, mpiInfo, (double) n);
		// 	//gWriteH5(E, mpiInfo, (double) n);
		// }

		//gWriteH5(phi, mpiInfo, (double) n);
		pWriteTemperature(temperature,pop,(double)n,units,ini);
		//fillGridIndexes(phi);
		//MPI_Barrier(MPI_COMM_WORLD);
		xyzWriteProbe(probe, phi,(double)n,mpiInfo);
		pWriteEnergy(history,pop,(double)n,units);
		//gWriteH5(phi, mpiInfo, (double) n);
		//free(nt);
		//free(mccTimestep);
		//free(frequency);
		//free(velThermal);
		//msg(STATUS, "   -    ");


	}
	//msg(STATUS, "Test returned %d", errorvar);

	if(mpiInfo->mpiRank==0) tMsg(t->total, "Time spent: ");
	MPI_Barrier(MPI_COMM_WORLD);

	/*
	* FINALIZE PINC VARIABLES
	*/
	gFreeMpi(mpiInfo);

	// Close h5 files
	pCloseH5(pop);
	gCloseH5(rho);
	gCloseH5(rho_e);
	gCloseH5(rho_i);
	gCloseH5(phi);
	gCloseH5(E);
	// oCloseH5(obj);
	xyCloseH5(history);
	xyCloseH5(temperature);
	arrCloseH5(probe);

	// Free memory
	mccFreeVars(mccVars);
	solverFree(solver);
	gFree(rho);
	gFree(rho_e);
	gFree(rho_i);
	gFree(phi);
	gFree(E);
	pFree(pop);
	free(S);
	free(T);
	// oFree(obj);

	gsl_rng_free(rngSync);
	gsl_rng_free(rng);

}


/// Old stuffz, u can delete if u are not me!

//
// funPtr mccMode2_set(dictionary *ini){
// 	//test sanity here!
// 	mccSanity(ini,"mccMode",2);
// 	return mccMode2;
// }
//
//
// void mccMode2(dictionary *ini){
//
//
// 	/*
// 	 * SELECT METHODS
// 	 */
// 	void (*acc)()   = select(ini,"methods:acc",	puAcc3D1_set,
// 												puAcc3D1KE_set,
// 												puAccND1_set,
// 												puAccND1KE_set,
// 												puAccND0_set,
// 												puAccND0KE_set,
// 												puBoris3D1_set,
// 												puBoris3D1KE_set,
// 												puBoris3D1KETEST_set);
//
// 	void (*distr)() = select(ini,"methods:distr",	puDistr3D1_set,
// 													puDistr3D1split_set,
// 													puDistrND1_set,
// 													puDistrND0_set);
//
//
// 	void (*extractEmigrants)() = select(ini,"methods:migrate",	puExtractEmigrants3D_set,
// 																puExtractEmigrantsND_set);
//
// 	//
// 	void (*solverInterface)()	= select(ini,	"methods:poisson",
// 												mgSolver_set
// 												//sSolver_set
// 												);
//
//
// 	//
// 	void (*solve)() = NULL;
// 	void *(*solverAlloc)() = NULL;
// 	void (*solverFree)() = NULL;
// 	solverInterface(&solve, &solverAlloc, &solverFree);
//
// 	/*
// 	 * INITIALIZE PINC VARIABLES
// 	 */
//
//
// 	MpiInfo *mpiInfo = gAllocMpi(ini);
// 	Population *pop = pAlloc(ini);
// 	Grid *E   = gAlloc(ini, VECTOR);
// 	Grid *rho = gAlloc(ini, SCALAR);
// 	Grid *res = gAlloc(ini, SCALAR);
// 	Grid *phi = gAlloc(ini, SCALAR);
// 	void *solver = solverAlloc(ini, rho, phi);
// 	/*
// 	* mcc specific variables
// 	*/
// 	// make struct???
// 	double PmaxElectron = 0;
// 	double maxfreqElectron = 0;
// 	double PmaxIon = 0;
// 	double maxfreqIon = 0;
// 	double *mass = pop->mass;
// 	int nSpecies = pop->nSpecies;
// 	double nt = iniGetDouble(ini,"collisions:numberDensityNeutrals"); //constant for now
// 	double mccTimestep = iniGetDouble(ini,"time:timeStep");
// 	//double frequency = iniGetDouble(ini,"collisions:collisionFrequency"); // for static Pmax...
// 	double *velThermal = iniGetDoubleArr(ini,"population:thermalVelocity",nSpecies);
//
//
// 	// using Boris algo
// 	//int nSpecies = pop->nSpecies;
// 	double *S = (double*)malloc((3)*(nSpecies)*sizeof(*S));
// 	double *T = (double*)malloc((3)*(nSpecies)*sizeof(*T));
//
// 	// Creating a neighbourhood in the rho to handle migrants
// 	gCreateNeighborhood(ini, mpiInfo, rho);
//
// 	// Setting Boundary slices
// 	gSetBndSlices(phi, mpiInfo);
//
// 	//Set mgSolve
// 	//MgAlgo mgAlgo = getMgAlgo(ini);
//
// 	// Random number seeds
// 	gsl_rng *rngSync = gsl_rng_alloc(gsl_rng_mt19937);
// 	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
// 	gsl_rng_set(rng,mpiInfo->mpiRank+1); // Seed needs to be >=1
//
//
// 	/*
// 	 * PREPARE FILES FOR WRITING
// 	 */
// 	int rank = phi->rank;
// 	double *denorm = malloc((rank-1)*sizeof(*denorm));
// 	double *dimen = malloc((rank-1)*sizeof(*dimen));
//
// 	for(int d = 1; d < rank;d++) denorm[d-1] = 1.;
// 	for(int d = 1; d < rank;d++) dimen[d-1] = 1.;
//
// 	pOpenH5(ini, pop, "pop");
// 	//gOpenH5(ini, rho, mpiInfo, denorm, dimen, "rho");
// 	gOpenH5(ini, phi, mpiInfo, denorm, dimen, "phi");
// 	//gOpenH5(ini, E,   mpiInfo, denorm, dimen, "E");
//   // oOpenH5(ini, obj, mpiInfo, denorm, dimen, "test");
//   // oReadH5(obj, mpiInfo);
//
// 	hid_t history = xyOpenH5(ini,"history");
// 	pCreateEnergyDatasets(history,pop);
//
// 	// Add more time series to history if you want
// 	// xyCreateDataset(history,"/group/group/dataset");
//
// 	free(denorm);
// 	free(dimen);
//
// 	/*
// 	 * INITIAL CONDITIONS
// 	 */
//
// 	// Initalize particles
// 	// pPosUniform(ini, pop, mpiInfo, rngSync);
// 	//pPosLattice(ini, pop, mpiInfo);
// 	//pVelZero(pop);
// 	// pVelMaxwell(ini, pop, rng);
// 	double maxVel = iniGetDouble(ini,"population:maxVel");
//
// 	// Manually initialize a single particle
// 	if(mpiInfo->mpiRank==0){
// 		double pos[3] = {8., 8., 8.};
// 		double vel[3] = {0.02, 0., 1.};
// 		pNew(pop, 0, pos, vel);
// 		double pos1[3] = {17., 17., 16.};
// 		double vel1[3] = {0.1, 0., 0.1};
// 		pNew(pop, 1, pos1, vel1); //second particle
// 	}
//
// 	// Perturb particles
// 	//pPosPerturb(ini, pop, mpiInfo);
//
// 	// Migrate those out-of-bounds due to perturbation
// 	extractEmigrants(pop, mpiInfo);
// 	puMigrate(pop, mpiInfo, rho);
//
// 	/*
// 	 * INITIALIZATION (E.g. half-step)
// 	 */
//
// 	// Get initial charge density
// 	distr(pop, rho);
// 	gHaloOp(addSlice, rho, mpiInfo, FROMHALO);
//
// 	// Get initial E-field
// 	solve(solver, rho, phi, mpiInfo);
// 	//solve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);
// 	gFinDiff1st(phi, E);
// 	gHaloOp(setSlice, E, mpiInfo, TOHALO);
// 	gMul(E, -1.);
//
//
// 	// Advance velocities half a step
//
//
// 	gMul(E, 0.5);
// 	puGet3DRotationParameters(ini, T, S, 0.5);
// 	// adScale(T, 3*nSpecies, 0.5);
// 	// adScale(S, 3*nSpecies, 0.5);
//
// 	gZero(E);
// 	puAddEext(ini, pop, E);
// 	acc(pop, E, T, S);
//
// 	gMul(E, 2.0);
// 	puGet3DRotationParameters(ini, T, S, 1.0);
// 	// adScale(T, 3*nSpecies, 2.0);
// 	// adScale(S, 3*nSpecies, 2.0);
//
// 	/*
// 	 * TIME LOOP
// 	 */
//
// 	Timer *t = tAlloc(rank);
//
// 	double x_min = 20;
// 	double x_max = 0;
// 	double y_min = 20;
// 	double y_max = 0;
//
// 	// n should start at 1 since that's the timestep we have after the first
// 	// iteration (i.e. when storing H5-files).
// 	PmaxElectron = 1.0;
// 	PmaxIon = 1.0; // first timestep
// 	int nTimeSteps = iniGetInt(ini,"time:nTimeSteps");
// 	for(int n = 1; n <= nTimeSteps; n++){
//
// 		msg(STATUS,"Computing time-step %i",n);
// 		MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary
//
// 		// Check that no particle moves beyond a cell (mostly for debugging)
// 		pVelAssertMax(pop,maxVel);
//
// 		tStart(t);
//
// 		// Move particles
// 		//adPrint(pop->pos, 3);
// 		//x_min = pop->pos[0]<x_min ? pop->pos[0] : x_min;
// 		//x_max = pop->pos[0]>x_max ? pop->pos[0] : x_max;
// 		//y_min = pop->pos[1]<y_min ? pop->pos[1] : y_min;
// 		//y_max = pop->pos[1]>y_max ? pop->pos[1] : y_max;
// 		puMove(pop);
// 		// oRayTrace(pop, obj);
//
// 		// Migrate particles (periodic boundaries)
// 		extractEmigrants(pop, mpiInfo);
// 		puMigrate(pop, mpiInfo, rho);
//
// 		mccGetPmaxElectron(ini,mass[0], velThermal[0], mccTimestep, nt,
// 			&PmaxElectron, &maxfreqElectron, pop,mpiInfo);
// 		mccGetPmaxIon(ini,mass[1], velThermal[1], mccTimestep, nt, &PmaxIon,
// 			&maxfreqIon, pop,mpiInfo);
//
// 		PmaxElectron = 0.0;
// 		PmaxIon = 0.0;
//
// 		if(n%4==0 ){
// 			PmaxIon = 1.0; // CEX coll
//
// 		}
//
//
// 		//MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary
// 		//mccCollideElectron(ini,pop, &PmaxElectron, &maxfreqElectron, rng, mccTimestep, nt,DebyeLength); //race conditions?????
// 		//MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary
// 		mccCollideIon_debug(ini, pop, &PmaxIon, &maxfreqIon, rng, mccTimestep, nt);
// 		//MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary
//
//
//
// 		// Check that no particle resides out-of-bounds (just for debugging)
// 		pPosAssertInLocalFrame(pop, rho);
// 		//msg(STATUS, "HEEEEEEEEEEERRRRRRRRRRREEEEEEEEEEEE");
// 		// Compute charge density
// 		distr(pop, rho);
// 		gHaloOp(addSlice, rho, mpiInfo, FROMHALO);
//
// 		// gAssertNeutralGrid(rho, mpiInfo);
//
// 		// Compute electric potential phi
// 		solve(solver, rho, phi, mpiInfo);
// 		//solve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);
//
// 		gAssertNeutralGrid(phi, mpiInfo);
//
// 		// Compute E-field
// 		gFinDiff1st(phi, E);
// 		gHaloOp(setSlice, E, mpiInfo, TOHALO);
// 		gMul(E, -1.);
//
// 		gAssertNeutralGrid(E, mpiInfo);
// 		// Apply external E
// 		// gAddTo(Ext);
// 		//puAddEext(ini, pop, E);
//
// 		// Accelerate particle and compute kinetic energy for step n
//
// 		gZero(E); // Turning off E-field for testing
// 		puAddEext(ini, pop, E);
// 		acc(pop, E, T, S);
//
// 		tStop(t);
//
// 		// Sum energy for all species
// 		pSumKinEnergy(pop);
//
// 		// Compute potential energy for step n
// 		gPotEnergy(rho,phi,pop);
//
// 		// Example of writing another dataset to history.xy.h5
// 		// xyWrite(history,"/group/group/dataset",(double)n,value,MPI_SUM);
//
// 		//Write h5 files
// 		//gWriteH5(E, mpiInfo, (double) n);
// 		//gWriteH5(rho, mpiInfo, (double) n);
// 		gWriteH5(phi, mpiInfo, (double) n);
// 		pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5);
// 		//pWriteEnergy(history,pop,(double)n);
//
// 	}
//
// 	msg(STATUS, "x in [%f, %f]", x_min, x_max);
// 	msg(STATUS, "y in [%f, %f]", y_min, y_max);
//
// 	if(mpiInfo->mpiRank==0) tMsg(t->total, "Time spent: ");
//
// 	/*
// 	 * FINALIZE PINC VARIABLES
// 	 */
// 	free(S);
// 	free(T);
//
// 	gFreeMpi(mpiInfo);
//
// 	// Close h5 files
// 	pCloseH5(pop);
// 	//gCloseH5(rho);
// 	gCloseH5(phi);
// 	//gCloseH5(E);
// 	// oCloseH5(obj);
// 	xyCloseH5(history);
//
// 	// Free memory
// 	solverFree(solver);
// 	gFree(rho);
// 	gFree(phi);
// 	gFree(E);
// 	pFree(pop);
// 	// oFree(obj);
//
// 	gsl_rng_free(rngSync);
// 	gsl_rng_free(rng);
//
// }
