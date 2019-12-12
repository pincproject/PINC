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

	//we use computational particles that contain many real particles.
	iniSetDouble(ini,"collisions:numberDensityNeutrals",nt);

	// in m/s
	double NvelThermal = iniGetDouble(ini,"collisions:thermalVelocityNeutrals");
	NvelThermal /= units->velocity;//(units->length/units->time);
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

	// "a" parameter is max crossect (m^2)
	// "b" parameter decides velocity to center about.
	// "b" given as 1/v^2

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
	}
	return MaxVelocity;
}

double mccGetMaxVelTran(const Population *pop, int species,const gsl_rng *rng, double NvelThermal){

	// iterate over pop and keep max velocity, uses
	double *vel = pop->vel;
	double MaxVelocity = 0;
	double NewVelocity = 0;
	int nDims = pop->nDims;
	double vxMW, vyMW, vzMW;

	long int iStart = pop->iStart[species];
	long int iStop  = pop->iStop[species];
	for(int i=iStart; i<iStop; i++){
		vxMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal); //maxwellian dist?
		vyMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal); //yes, gaussian in 3 dim
		vzMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal); // is maxwellian

		NewVelocity = sqrt((vel[i*nDims]-vxMW)*(vel[i*nDims]*vxMW)+(vel[i*nDims+1]-vyMW)\
			*(vel[i*nDims+1]-vyMW)+(vel[i*nDims+2]-vzMW)*(vel[i*nDims+2]-vzMW));
		if(NewVelocity>MaxVelocity){
			MaxVelocity=NewVelocity;
		}
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
	double *neutralDrift = iniGetDoubleArr(ini, "population:drift", 3*nSpecies);
	electronMassRatio = mass[0]*units->mass/units->weights[0];
	electronMassRatio /= iniGetDouble(ini,"collisions:realElectronMass");
	mccVars->neutralDrift = neutralDrift,
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

	//determines maximum collision probability

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
	Population *pop, MpiInfo *mpiInfo,const gsl_rng *rng){

	// Faster static version. uses static cross sections
	// to determine maximum collision probability

	double NvelThermal = mccVars->NvelThermal;
	double nt = mccVars->nt;
	double StaticSigmaCEX = mccVars->mccSigmaCEX;
	double StaticSigmaIonElastic = mccVars->mccSigmaIonElastic;
	double max_v = mccGetMaxVelTran(pop,1,rng,NvelThermal);
	//double max_v = mccGetMaxVel(pop,1);
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
	// to determine maximum collision probability

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
	double *drift = mccVars->neutralDrift;

	double angleChi = 0; //scattering angle in x, y
	double anglePhi = 0; //scattering angle in j
	double angleTheta = 0;
	double A = 0; //scaling Factor used for precomputing (speedup)
	double B = 0;
	double Ekin, newEkin;
	double R2 = gsl_rng_uniform_pos(rng);
	double R3 = gsl_rng_uniform_pos(rng);

	// let colls happen in frame comoving with drift
	double vx = *vx_point-drift[0];
	double vy = *vy_point-drift[1];
	double vz = *vz_point-drift[2];
	double vx_,vy_,vz_;

	// unit velocities
	double velchange = 0.0;// change of energy incident particle
	double velsquare = sqrt(vx*vx+vy*vy+vz*vz);
	if(velsquare<0.0000000000000001){ // safety, mostly for debugging
		velsquare = 0.0000000000000001;
	}

	Ekin = mccVars->energyConvFactor*0.5*(vx*vx + vy*vy + vz*vz)*mass[0]; //eV

	//make unit vector
	vx = vx/velsquare;
	vy = vy/velsquare;
	vz = vz/velsquare;

	// angles
	double argument = (2+Ekin-2*pow((1+Ekin),R2))/(Ekin);
	Ekin = Ekin/(mccVars->energyConvFactor); //energy in eV

	if(sqrt(argument*argument)>1.0){
		double newargument = sqrt(argument*argument)/argument;
		msg(WARNING,"old argument acos = %.32f new arg = %.64f, R = %.32f",
		argument, newargument,R2);
		argument = newargument;
	}
	angleChi =  acos(argument); // gives nan value if abs(argument) > 1
	newEkin = Ekin*(1-((2.0*mass[0])/(artificialLoss*mass[1]))
	*(1-(cos(angleChi))));

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

	// inc - scat vector relation
	//vector relation
	A = (sin(angleChi)*sin(anglePhi))/sin(angleTheta);
	B = (sin(angleChi)*cos(anglePhi))/sin(angleTheta);
	vx_ = (vx*cos(angleChi)+B*(vy*vy + vz*vz) ); //Vx
	vy_ = (vy*cos(angleChi)+A*vz-B*vx*vy); //Vy
	vz_ = (vz*cos(angleChi)-A*vy-B*vx*vz); //Vz

	// change in energy
	vx_ *= velchange;
	vy_ *= velchange;
	vz_ *= velchange;

	*vx_point=vx_+drift[0];
	*vy_point=vy_+drift[1];
	*vz_point=vz_+drift[2];

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
	double vxMW, double vyMW, double vzMW,
	const gsl_rng *rng, MccVars *mccVars, Population *pop){


	double *mass = pop->mass;
	double *drift = mccVars->neutralDrift;

	double angleChi = 0; //scattering angle in x, y
	double anglePhi = 0; //scattering angle in j
	double angleTheta = 0;
	double A = 0; //scaling Factor used for precomputing (speedup)
	double B = 0;
	double Ekin, newEkin;
	double velchange;

    // let colls happen in frame comoving with drift
	double vx = *vx_point;
	double vy = *vy_point;
	double vz = *vz_point;
	double vx_=0;
	double vy_=0;
	double vz_=0;
	double R1 = gsl_rng_uniform_pos(rng);
	double R2 = gsl_rng_uniform_pos(rng);


	//vxMW = vxMW;
	//vyMW = vyMW;
	//vzMW = vzMW;
	//printf("drift = %f, %f, %f \n",drift[3],drift[4],drift[5] );
	//printf("vxMW = %f, %f, %f \n",vxMW,vyMW,vzMW );
	//printf("vx = %f, %f, %f \n",vx,vy,vz );

	//The scattering happens in netral rest frame
	//transfer to neutral stationary frame
	double vxTran = vx-vxMW;
	double vyTran = vy-vyMW;
	double vzTran = vz-vzMW;

	double velsquare = sqrt(vxTran*vxTran+vyTran*vyTran+vzTran*vzTran);
	if(velsquare<0.00000000001){ // safety, mostly for debugging
		velsquare = 0.00000000001;
	}

	//make unit vector
	vx = vxTran/velsquare;
	vy = vyTran/velsquare;
	vz = vzTran/velsquare;

	// angles assuming neutral and ion have same mass.
	Ekin = 0.5*(vxTran*vxTran+vyTran*vyTran+vzTran*vzTran)*mass[1]; //PINC neutral frame
	angleChi =  acos(sqrt(1-R1)); //labframe goes from 0 to pi/2
	newEkin = Ekin*(cos(angleChi)*cos(angleChi));
	velchange = sqrt( (2.0*newEkin)/mass[1] );
	anglePhi = 2*PI*R2;

	if(vx>0.0){ // debugging, should not be necessary
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

	double ekinafter = 0.5*(vx_*vx_ + vy_*vy_ + vz_*vz_)*mass[1];

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
	double energydiff = Ekin-ekinafter;
	if(abs(energydiffexpr-energydiff)>1e-31){
		msg(WARNING,"too large energy error in collide Ion = %e",
		(abs(energydiffexpr-energydiff)));
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

	double nt = mccVars->nt;
	double mccSigmaElectronElastic = mccVars->mccSigmaElectronElastic;
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


	long int iStart = pop->iStart[0];
	long int iStop = pop->iStop[0];
	long int NparticleColl = (Pmax)*(pop->iStop[0]-pop->iStart[0]);
	long int mccStepSize = floor((double)((iStop-iStart))\
	/ (double)(NparticleColl));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "colliding %i of %i electrons\n",
		(NparticleColl), ((iStop-iStart)));
	}

	//enshure non out of bounds
	long int mccStop = iStart + mccStepSize*NparticleColl;

	for(long int i=iStart;i<mccStop;i+=mccStepSize){
		R = gsl_rng_uniform_pos(rng); // New random number per particle.
		Rp = gsl_rng_uniform_pos(rng); // separate rand num. for prob.
		q = ((i + floor(R*mccStepSize))*nDims);

		vx = &vel[q];
		vy = &vel[q+1];
		vz = &vel[q+2];

		double MyCollFreq = mccGetMyCollFreqStatic(mccSigmaElectronElastic,  \
			vel[q],vel[q+1],vel[q+2],nt);

		if (Rp<(MyCollFreq/maxfreqElectron)){
		// Pmax gives the max possible colls. but we will "never" reach this
		// number, so we need to test for each particles probability.
			errorcounter += 1;
			scatterElectron(vx,vy,vz,rng,mccVars,pop);
		}
		last_i = i;
	}
	// Special handling of last box to let every particle have posibillity
	// to collide

	R = gsl_rng_uniform_pos(rng); // New random number per particle.
	Rp = gsl_rng_uniform_pos(rng); // separate rand num. for prob.
	q = ( (last_i) + floor(R*(((iStop)-last_i))) )*nDims; // particle q collides
	vx = &vel[q];
	vy = &vel[q+1];
	vz = &vel[q+2];

	double MyCollFreq = mccGetMyCollFreqStatic(mccSigmaElectronElastic,  \
		vel[q],vel[q+1],vel[q+2],nt);

	if (Rp<(MyCollFreq/maxfreqElectron)){
		errorcounter += 1;
		scatterElectron(vx,vy,vz,rng,mccVars,pop);

	}
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision",
		"counted  %i Electron collisions on one MPI node \n",  errorcounter);
	}
}

void mccCollideIonStatic(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng, MpiInfo *mpiInfo){

	// uses static CROSS-Sections, collfreq is proportional to v

	mccGetPmaxIonStatic(ini,mccVars,pop,mpiInfo,rng);

	double nt = mccVars->nt;
	double NvelThermal = mccVars->NvelThermal;
	double mccSigmaCEX= mccVars->mccSigmaCEX;
	double mccSigmaIonElastic = mccVars->mccSigmaIonElastic;
	double maxfreqIon = mccVars->maxFreqIon;
	double Pmax = mccVars->pMaxIon;
	int nDims = pop->nDims;
	double *vel = pop->vel;
	double Rp, Rq;

	// pointers to pass to scatter function
	double* vx;
	double* vy;
	double* vz;

	long int q = 0;
	long int errorcounter = 0; // count collisions
	long int errorcounter1 = 0; // count cex collisions
	double vxMW, vyMW, vzMW;
	long int last_i = 0;

	long int iStart = pop->iStart[1];
	long int iStop = pop->iStop[1];
	long int NparticleColl = (Pmax)*(pop->iStop[1]-pop->iStart[1]);
	long int mccStepSize = (floor((double)((iStop-iStart))\
	/ (double)(NparticleColl)));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "colliding %i of %i ions \n", (NparticleColl), ((iStop-iStart)));
	}

	//enshure non out of bounds
	long int mccStop = iStart + mccStepSize*NparticleColl;
	for(long int i=iStart;i<mccStop;i+=mccStepSize){


		vxMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal);
		vyMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal);
		vzMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal);

		Rp = gsl_rng_uniform_pos(rng); //decides type of coll.
		Rq = gsl_rng_uniform_pos(rng);
		q = (i + floor(Rq*mccStepSize))*nDims; // particle q collides

		//transfer to neutral stationary frame
		double vxTran = vel[q]-vxMW;
		double vyTran = vel[q+1]-vyMW;
		double vzTran = vel[q+2]-vzMW;

		double MyCollFreq1 = mccGetMyCollFreqStatic(mccSigmaIonElastic,vxTran,
			vyTran,vzTran, nt);
		double MyCollFreq2 = mccGetMyCollFreqStatic(mccSigmaCEX,vxTran,
			vyTran,vzTran, nt);

		// MyCollFreq1+MyCollFreq2 defines total coll freq.
		if (Rp<( (MyCollFreq1+MyCollFreq2)/maxfreqIon)){
			if(Rp < MyCollFreq1/maxfreqIon){
				// elastic:
				errorcounter += 1;

				//point to new velocity
				vx = &vel[q];
				vy = &vel[q+1];
				vz = &vel[q+2];
				scatterIon(vx,vy,vz,vxMW,vyMW,vzMW,rng,mccVars,pop);
			}else{
				errorcounter += 1;
				errorcounter1 += 1;
				// Charge exchange:
				// flip ion and neutral. Neutral is new ion.
				vel[q] = vxMW; //Vx
				vel[q+1] = vyMW; //Vy
				vel[q+2] = vzMW; //Vz
			}
		}
		last_i = i;
	}

	// Special handling of last box to let every particle
	// have posibillity to collide

	Rp = gsl_rng_uniform_pos(rng); //decides type of coll.
	Rq = gsl_rng_uniform_pos(rng);
	q = ( (last_i) + floor(Rq*(((iStop)-last_i))) )*nDims; //particle q collides

	vxMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal);
	vyMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal);
	vzMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal);

	//transfer to neutral stationary frame
	double vxTran = vel[q]-vxMW;
	double vyTran = vel[q+1]-vyMW;
	double vzTran = vel[q+2]-vzMW;

	double MyCollFreq1 = mccGetMyCollFreqStatic(mccSigmaIonElastic,vxTran,
		vyTran,vzTran, nt);
	double MyCollFreq2 = mccGetMyCollFreqStatic(mccSigmaCEX,vxTran,
		vyTran,vzTran, nt); //prob of coll for particle i

	if (Rp<( (MyCollFreq1+MyCollFreq2)/maxfreqIon)){ // MyCollFreq1+MyCollFreq2 defines total coll freq.
		if(Rp < MyCollFreq1/maxfreqIon){ // if test slow, but enshures randomness...
			// elastic:
			errorcounter += 1;

			//point to new velocity
			vx = &vel[q];
			vy = &vel[q+1];
			vz = &vel[q+2];
			scatterIon(vx,vy,vz,vxMW,vyMW,vzMW,rng,mccVars,pop);
		}else{
			// Charge exchange:
			// flip ion and neutral. Neutral is new ion.
			errorcounter += 1;
			errorcounter1 += 1;
			vel[q] = vxMW;  //Vx
			vel[q+1] = vyMW;   //Vy
			vel[q+2] = vzMW;  //Vz

		}
	}
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision",
		"counted  %i ION collisions on one MPI node, %i as CEX \n"
		, errorcounter,errorcounter1);
		//printf("counted  %i Ion collisions on one MPI node \n",  errorcounter);
	}
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

	long int iStart = pop->iStart[0];
	long int iStop = pop->iStop[0];
	long int NparticleColl = (Pmax)*(pop->iStop[0]-pop->iStart[0]);
	long int mccStepSize = floor((double)((iStop-iStart))\
	/ (double)(NparticleColl));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "colliding %i of %i electrons\n",
		(NparticleColl), ((iStop-iStart)));
	}

	//enshure non out of bounds
	long int mccStop = iStart + mccStepSize*NparticleColl;
	for(long int i=iStart;i<mccStop;i+=mccStepSize){
		R = gsl_rng_uniform_pos(rng); // New random number per particle.
		Rp = gsl_rng_uniform_pos(rng); // separate rand num. for prob.
		q = ((i + floor(R*mccStepSize))*nDims);

		vx = &vel[q];
		vy = &vel[q+1];
		vz = &vel[q+2];

		double MyCollFreq = mccGetMyCollFreqFunctional(mccSigmaElectronElasticFunctional,  \
			vel[q],vel[q+1],vel[q+2],nt,electron_a,electron_b);

		if (Rp<(MyCollFreq/maxfreqElectron)){
			errorcounter += 1;
			scatterElectron(vx,vy,vz,rng,mccVars,pop);
		}
		if(mpiInfo->mpiRank==0){
			fMsg(ini, "collision", "counted  %i Electron collisions on one MPI node \n",  errorcounter);
		}
		last_i = i;
	}
	// Special handling of last box to let every particle have
	// posibillity to collide

	R = gsl_rng_uniform_pos(rng); // New random number per particle.
	Rp = gsl_rng_uniform_pos(rng); // separate rand num. for prob.
	q = ( (last_i) + floor(R*(((iStop)-last_i))) )*nDims; // particle q collides

	vx = &vel[q];
	vy = &vel[q+1];
	vz = &vel[q+2];

	double MyCollFreq = mccGetMyCollFreqFunctional(mccSigmaElectronElasticFunctional,  \
		vel[q],vel[q+1],vel[q+2],nt,electron_a,electron_b); //prob of coll for particle i

		if (Rp<(MyCollFreq/maxfreqElectron)){
			errorcounter += 1;
			scatterElectron(vx,vy,vz,rng,mccVars,pop);
		}
		if(mpiInfo->mpiRank==0){
			fMsg(ini, "collision", "counted  %i Electron collisions on one MPI node \n",  errorcounter);
		}
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
	double vxMW, vyMW, vzMW;

	long int errorcounter = 0; // count collisions
	long int errorcounter1 = 0; // count cex collisions

	// pointers to pass to scatter function
	double* vx;
	double* vy;
	double* vz;

	long int iStart = pop->iStart[1];
	long int iStop = pop->iStop[1];
	long int NparticleColl = (Pmax)*(pop->iStop[1]-pop->iStart[1]);
	long int mccStepSize = floor((double)((iStop-iStart))\
	/ (double)(NparticleColl));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "colliding %i of %i ions \n", (NparticleColl), ((iStop-iStart)));
	}

	//enshure non out of bounds
	long int mccStop = iStart + mccStepSize*NparticleColl;

	for(long int i=iStart;i<mccStop;i+=mccStepSize){

		vxMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal);
		vyMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal);
		vzMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal);

		//transfer to neutral stationary frame
		double vxTran = vel[q]-vxMW;
		double vyTran = vel[q+1]-vyMW;
		double vzTran = vel[q+2]-vzMW;

		Rp = gsl_rng_uniform_pos(rng); //decides type of coll.
		Rq = gsl_rng_uniform_pos(rng);
		q = (i + floor(Rq*mccStepSize))*nDims; // particle q collides

		double MyCollFreq1 = mccGetMyCollFreqFunctional(
			mccSigmaIonElasticFunctional,vxTran,
				vyTran,vzTran,nt,ion_elastic_a,ion_elastic_b);
		double MyCollFreq2 = mccGetMyCollFreqFunctional(mccSigmaCEXFunctional,
			vxTran,vyTran,vzTran,nt,CEX_a,CEX_b);

		if (Rp<( (MyCollFreq1+MyCollFreq2)/maxfreqIon)){
			if(Rp < MyCollFreq1/maxfreqIon){
				// elastic:
				errorcounter += 1;
				//point to new velocity
				vx = &vel[q];
				vy = &vel[q+1];
				vz = &vel[q+2];
				scatterIon(vx,vy,vz,vxMW,vyMW,vzMW,rng,mccVars,pop);
			}else{
				errorcounter += 1;
				errorcounter1 += 1;
				// Charge exchange:
				// flip ion and neutral. Neutral is new ion.
				vel[q] = vxMW;
				vel[q+1] = vyMW;
				vel[q+2] = vzMW;
			}
		}
		last_i = i;
	}
	// Special handling of last box to let every particle have posibillity to collide

	vxMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal);
	vyMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal);
	vzMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal);

	//transfer to neutral stationary frame
	double vxTran = vel[q]-vxMW;
	double vyTran = vel[q+1]-vyMW;
	double vzTran = vel[q+2]-vzMW;

	Rp = gsl_rng_uniform_pos(rng); //decides type of coll.
	Rq = gsl_rng_uniform_pos(rng);
	q = ( (last_i) + floor(Rq*(((iStop)-last_i))) )*nDims; //particle q collides

	double MyCollFreq1 = mccGetMyCollFreqFunctional(
		mccSigmaIonElasticFunctional,vxTran,
			vyTran,vzTran,nt,ion_elastic_a,ion_elastic_b);
	double MyCollFreq2 = mccGetMyCollFreqFunctional(mccSigmaCEXFunctional,
		vxTran,vyTran,vzTran,nt,CEX_a,CEX_b);

	if (Rp<( (MyCollFreq1+MyCollFreq2)/maxfreqIon)){
		if(Rp < MyCollFreq1/maxfreqIon){
			// elastic:
			//msg(STATUS,"elastic");
			errorcounter += 1;

			//point to new velocity
			vx = &vel[q];
			vy = &vel[q+1];
			vz = &vel[q+2];

			scatterIon(vx,vy,vz,vxMW,vyMW,vzMW,rng,mccVars,pop);
		}else{
			// Charge exchange:
			// flip ion and neutral. Neutral is new ion.
			errorcounter += 1;
			errorcounter1 += 1;
			vel[q] = vxMW;
			vel[q+1] = vyMW;
			vel[q+2] = vzMW;

		}
	}
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "counted  %i ION collisions on one MPI node, %i as CEX \n", errorcounter,errorcounter1);
	}
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
	long int last_i = 0;
	long int errorcounter = 0;

	long int iStart = pop->iStart[0];
	long int iStop = pop->iStop[0];
	long int NparticleColl = (Pmax)*(pop->iStop[0]-pop->iStart[0]);
	long int mccStepSize = floor((double)((iStop-iStart))\
	/ (double)(NparticleColl));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "colliding %i of %i electrons\n", (NparticleColl), ((iStop-iStart)));
	}
	//enshure non out of bounds
	long int mccStop = iStart + mccStepSize*NparticleColl;
	for(long int i=iStart;i<mccStop;i+=mccStepSize){

		R = gsl_rng_uniform_pos(rng); // New random number per particle. maybe print?
		q = ((i + floor(R*mccStepSize))*nDims);

		vx = &vel[q];
		vy = &vel[q+1];
		vz = &vel[q+2];

		errorcounter += 1;
		scatterElectron(vx,vy,vz,rng,mccVars,pop);
		last_i = i;
	}

	// Special handling of last box to let every particle have posibillity to collide

	R = gsl_rng_uniform_pos(rng); // New random number per particle. maybe print?
	q = ( (last_i) + floor(R*(((iStop)-last_i))) )*nDims; // particle q collides


	vx = &vel[q];
	vy = &vel[q+1];
	vz = &vel[q+2];

	errorcounter += 1;
	scatterElectron(vx,vy,vz,rng,mccVars,pop);

	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "counted  %i Electron collisions on one MPI node \n",  errorcounter);
	}
}

void mccCollideIonConstantFrq(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng, MpiInfo *mpiInfo){

	// uses constant collfreq fixed number of colls per dt

	double NvelThermal = mccVars->NvelThermal;
	double collFrqIonElastic = mccVars->collFrqIonElastic;
	double collFrqIonCEX = mccVars->collFrqCex;

	mccGetPmaxIonConstantFrq(ini,mccVars,pop,mpiInfo);
	double Pmax = mccVars->pMaxIon;
	int nDims = pop->nDims;
	double *vel = pop->vel;
	double *drift = mccVars->neutralDrift;

	double Rp, Rq;
	long int q = 0;
	long int last_i = 0;

	long int errorcounter = 0; // count collisions
	long int errorcounter1 = 0; // count cex collisions

	double vxMW, vyMW, vzMW;

	double* vx;
	double* vy;
	double* vz;

	double MyCollFreq1 = collFrqIonElastic;
	double MyCollFreq2 = collFrqIonCEX;
	double maxfreqIon = MyCollFreq1+MyCollFreq2;

	long int iStart = pop->iStart[1];
	long int iStop = pop->iStop[1];
	long int NparticleColl = (Pmax)*(pop->iStop[1]-pop->iStart[1]);
	long int mccStepSize = floor((double)((iStop-iStart))\
	/ (double)(NparticleColl));
	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "colliding %i of %i ions \n", (NparticleColl), ((iStop-iStart)));
	}

	//enshure non out of bounds
	long int mccStop = iStart + mccStepSize*NparticleColl;

	for(long int i=iStart;i<mccStop;i+=mccStepSize){

		vxMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal)+drift[3];
		vyMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal)+drift[4];
		vzMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal)+drift[5];

		//printf("vxMW = %f, %f, %f \n",vxMW,vyMW,vzMW );

		Rp = gsl_rng_uniform_pos(rng); //decides type of coll.
		Rq = gsl_rng_uniform_pos(rng);
		q = (i + floor(Rq*mccStepSize))*nDims; // particle q collides

		if(Rp < MyCollFreq1/ maxfreqIon){
			// elastic:
			errorcounter += 1;

			//point to new velocity
			vx = &vel[q];
			vy = &vel[q+1];
			vz = &vel[q+2];

			scatterIon(vx,vy,vz,vxMW,vyMW,vzMW,rng,mccVars,pop);
		}else{
			errorcounter += 1;
			errorcounter1 += 1;
			// Charge exchange:
			// flip ion and neutral. Neutral is new ion.
			vel[q] = vxMW;
			vel[q+1] = vyMW;
			vel[q+2] = vzMW;
		}
		last_i = i;
	}
	// Special handling of last box to let every particle have posibillity to collide

	Rp = gsl_rng_uniform_pos(rng); //decides type of coll.
	Rq = gsl_rng_uniform_pos(rng);
	q = ( (last_i) + floor(Rq*(((iStop)-last_i))) )*nDims; //particle q collides

	vxMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal)+drift[3];
	vyMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal)+drift[4];
	vzMW = gsl_ran_gaussian_ziggurat(rng,NvelThermal)+drift[5];

	if(Rp < MyCollFreq1/ maxfreqIon){
		// elastic:
		errorcounter += 1;
		//point to new velocity
		vx = &vel[q];
		vy = &vel[q+1];
		vz = &vel[q+2];

		scatterIon(vx,vy,vz,vxMW,vyMW,vzMW,rng,mccVars,pop);
	}else{
		// Charge exchange:
		// flip ion and neutral. Neutral is new ion.
		errorcounter += 1;
		errorcounter1 += 1;

		vel[q] = vxMW;
		vel[q+1] = vyMW;
		vel[q+2] = vzMW;
	}

	if(mpiInfo->mpiRank==0){
		fMsg(ini, "collision", "counted  %i ION collisions on one MPI node, %i as CEX \n", errorcounter,errorcounter1);
	}
}

/*************************************************
*		Method Handlers
************************************************/

// Functional pointer that points to a method handler function.
//method handlers are basically a collection of functions to call
// in a method.

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
	//does nothing to turn off all collisions.
	NULL;
}

funPtr collissionsOff_set(dictionary *ini){
	mccSanity(ini,"mccMode",2);
	funPtr collissions;
	collissions = &collissionsOff;

	return collissions;
}



/*************************************************
*		RUNS
************************************************/

// runs used in development.

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
	Population *pop = pAlloc(ini,mpiInfo);
	Grid *phi = gAlloc(ini, SCALAR,mpiInfo);
	Grid *E   = gAlloc(ini, VECTOR, mpiInfo);
	Grid *rho = gAlloc(ini, SCALAR, mpiInfo);
	Grid *rho_e = gAlloc(ini, SCALAR, mpiInfo);
	Grid *rho_i = gAlloc(ini, SCALAR, mpiInfo);
	void *solver = solverAlloc(ini, rho, phi);


	/*
	* mcc specific variables
	*/

	MccVars *mccVars=mccAlloc(ini,units);

	// using Boris algo
	int nSpecies = pop->nSpecies;
	double *S = (double*)malloc((3)*(nSpecies)*sizeof(double));
	double *T = (double*)malloc((3)*(nSpecies)*sizeof(double));

	// Creating a neighbourhood in the rho to handle migrants
	gCreateNeighborhood(ini, mpiInfo, rho);

	// Setting Boundary slices
	gSetBndSlices(ini, phi, mpiInfo);

	//Set mgSolve
	//MgAlgo mgAlgo = getMgAlgo(ini);

	// Random number seeds
	gsl_rng *rngSync = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,mpiInfo->mpiRank+1); // Seed needs to be >=1

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
	distr(pop, rho, rho_e, rho_i); //two species
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
	gWriteH5(E, mpiInfo, 0.0);
	gWriteH5(rho, mpiInfo, 0.0);
	//gWriteH5(rho_e, mpiInfo, 0.0);
	//gWriteH5(rho_i, mpiInfo, 0.0);
	gWriteH5(phi, mpiInfo, 0.0);
	//pWriteH5(pop, mpiInfo, 0.0, 0.5,1);
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

		// Check that no particle moves beyond a cell (mostly for debugging)
		//pVelAssertMax(pop,maxVel);

		tStart(t);
		// Move particles
		puMove(pop);
		// oRayTrace(pop, obj);

		// Migrate particles (periodic boundaries)
		extractEmigrants(pop, mpiInfo);
		puMigrate(pop, mpiInfo, rho);

		/*
		*   Collisions
		*/
		collide(ini, pop, mccVars, rng,mpiInfo);

		// Check that no particle resides out-of-bounds (just for debugging)
		//pPosAssertInLocalFrame(pop, rho);

		// Compute charge density
		distr(pop, rho); //two species

		gHaloOp(addSlice, rho, mpiInfo, FROMHALO);
		gHaloOp(addSlice, rho_e, mpiInfo, FROMHALO);
		gHaloOp(addSlice, rho_i, mpiInfo, FROMHALO);

		//gAssertNeutralGrid(rho, mpiInfo);

		// Compute electric potential phi
		//solve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);
		solve(solver, rho, phi, mpiInfo);

		// Compute E-field
		gFinDiff1st(phi, E);
		gHaloOp(setSlice, E, mpiInfo, TOHALO);
		gMul(E, -1.);

		//gAssertNeutralGrid(E, mpiInfo);
		// Apply external E
		// gAddTo(Ext);
		//gZero(E); ////temporary test
		puAddEext(ini, pop, E);

		// Accelerate particle and compute kinetic energy for step n
		acc(pop, E, T, S);
		tStop(t);

		// Sum energy for all species
		pSumKinEnergy(pop); // kin_energy per species to pop
		// Compute potential energy for step n
		gPotEnergy(rho,phi,pop);

		if( n%500 == 0){
			// Example of writing another dataset to history.xy.h5
			// xyWrite(history,"/group/group/dataset",(double)n,value,MPI_SUM);
			//Write h5 files
			gWriteH5(E, mpiInfo, (double) n);
			gWriteH5(rho, mpiInfo, (double) n);
			//gWriteH5(rho_e, mpiInfo, (double) n);
			//gWriteH5(rho_i, mpiInfo, (double) n);
			gWriteH5(phi, mpiInfo, (double) n);
			pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5);
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

		if(n == 80000){
			pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5); //0.0001
		}


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
		xyzWriteProbe(probe, phi,(double)n,mpiInfo);
		pWriteEnergy(history,pop,(double)n,units);
		//gWriteH5(phi, mpiInfo, (double) n);

	}

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



funPtr oCollMode_set(dictionary *ini){
	return oCollMode;
}

void oCollMode(dictionary *ini){

	/*
	 * SELECT METHODS
	 */
	void (*acc)()   			= select(ini,	"methods:acc",
												puAcc3D1_set,
												puAcc3D1KE_set,
												puAccND1_set,
												puAccND1KE_set,
												puAccND0_set,
												puAccND0KE_set,
                        puBoris3D1KETEST_set);

	void (*distr)() 			= select(ini,	"methods:distr",
												puDistr3D1split_set,
												puDistr3D1_set,
												puDistrND1_set,
												puDistrND0_set);


	void (*collide)() = select(ini,	"methods:mcc",
									collissionsOff_set,
									constCrossect_set,
									constFreq_set,
									functionalCrossect_set);


	void (*extractEmigrants)()	= select(ini,	"methods:migrate",
												puExtractEmigrants3D_set,
												puExtractEmigrantsND_set,
                        						puExtractEmigrants3DOpen_set);

	void (*solverInterface)()	= select(ini,	"methods:poisson",
												mgSolver_set
												//sSolver_set
											);

	void (*solve)() = NULL;
	void *(*solverAlloc)() = NULL;
	void (*solverFree)() = NULL;
	solverInterface(&solve, &solverAlloc, &solverFree);

	/*
	 * INITIALIZE PINC VARIABLES
	 */
	Units *units=uAlloc(ini);
	uNormalize(ini, units);
	mccNormalize(ini,units);

	MpiInfo *mpiInfo = gAllocMpi(ini);
	Population *pop = pAlloc(ini,mpiInfo);
	Grid *E   = gAlloc(ini, VECTOR,mpiInfo);
	Grid *rho = gAlloc(ini, SCALAR,mpiInfo);
	Grid *rho_e = gAlloc(ini, SCALAR, mpiInfo);
	Grid *rho_i = gAlloc(ini, SCALAR, mpiInfo);
    Grid *rhoObj = gAlloc(ini, SCALAR,mpiInfo);     // for capMatrix - objects
	Grid *phi = gAlloc(ini, SCALAR,mpiInfo);
	void *solver = solverAlloc(ini, rho, phi, mpiInfo);
	MccVars *mccVars=mccAlloc(ini,units);

    Object *obj = oAlloc(ini,mpiInfo,units);              // for capMatrix - objects
//TODO: look into multigrid E,rho,rhoObj

	// Creating a neighbourhood in the rho to handle migrants
	gCreateNeighborhood(ini, mpiInfo, rho);

  // Setting Boundary slices
  gSetBndSlices(ini, phi, mpiInfo);

	// Random number seeds
	gsl_rng *rngSync = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,mpiInfo->mpiRank+1); // Seed needs to be >=1

	/*
	 * PREPARE FILES FOR WRITING
	 */

	pOpenH5(ini, pop, units, "pop");
	//double denorm = units->potential;
	gOpenH5(ini, rho, mpiInfo, units, units->chargeDensity, "rho");
	gOpenH5(ini, rho_e, mpiInfo, units, units->chargeDensity, "rho_e");
	gOpenH5(ini, rho_i, mpiInfo, units, units->chargeDensity, "rho_i");
	gOpenH5(ini, phi, mpiInfo, units, units->potential, "phi");
	gOpenH5(ini, E,   mpiInfo, units, units->eField, "E");
  // oOpenH5(ini, obj, mpiInfo, units, 1, "test");
  // oReadH5(obj, mpiInfo);

    //msg(STATUS,"opening obj file");
		gOpenH5(ini, rhoObj, mpiInfo, units, units->chargeDensity, "rhoObj");        // for capMatrix - objects
		//oOpenH5(ini, obj, mpiInfo, units, units->chargeDensity, "object");          // for capMatrix - objects
		//oReadH5(obj->domain, mpiInfo, "Object");

    //msg(STATUS,"done");


		//Count the number of objects and fill the lookup tables.
    //msg(STATUS,"filling lookup table");
		//oFillLookupTables(obj,mpiInfo);

    //msg(STATUS,"finding surface nodes");
		// Find all the object nodes which are part of the object surface.
		//oFindObjectSurfaceNodes(obj, mpiInfo);


	hid_t history = xyOpenH5(ini,"history");
	pCreateEnergyDatasets(history,pop);

	// Add more time series to history if you want
	// xyCreateDataset(history,"/group/group/dataset");

	/*
	 * INITIAL CONDITIONS
	 */

    //Compute capacitance matrix
    //msg(STATUS, "com cap matrix");
    oComputeCapacitanceMatrix(obj, ini, mpiInfo);




	// Initalize particles
	pPosUniform(ini, pop, mpiInfo, rngSync);
	//pPosLattice(ini, pop, mpiInfo);
	//pVelZero(pop);
	pVelMaxwell(ini, pop, rng);
	double maxVel = iniGetDouble(ini,"population:maxVel");

	//add influx of new particles on boundary
    pPurgeGhost(pop, rho);
    pFillGhost(ini,pop,rng,mpiInfo);

	// Perturb particles
	//pPosPerturb(ini, pop, mpiInfo);

	// Migrate those out-of-bounds due to perturbation
	extractEmigrants(pop, mpiInfo);

	puMigrate(pop, mpiInfo, rho);



	/*
	 * INITIALIZATION (E.g. half-step)
	 */

    // Clean objects from any charge first.
    gZero(rhoObj);                                          // for capMatrix - objects
    oCollectObjectCharge(pop, rhoObj, obj, mpiInfo);        // for capMatrix - objects
    gZero(rhoObj);                                          // for capMatrix - objects


	// Get initial charge density
	distr(pop, rho,rho_e,rho_i);
	gHaloOp(addSlice, rho, mpiInfo, FROMHALO);
	gHaloOp(addSlice, rho_e, mpiInfo, FROMHALO);
	gHaloOp(addSlice, rho_i, mpiInfo, FROMHALO);
    //gWriteH5(rho, mpiInfo, (double) 0);

	// Get initial E-field

  //gBnd(phi, mpiInfo);
	solve(solver, rho, phi, mpiInfo);
    //gWriteH5(phi, mpiInfo, (double) 0);
    //pWriteH5(pop, mpiInfo, (double) 0, (double)0+0.5);

	gFinDiff1st(phi, E);
	gHaloOp(setSlice, E, mpiInfo, TOHALO);
	gMul(E, -1.);


  //Boris parameters
  int nSpecies = pop->nSpecies;
	double *S = (double*)malloc((3)*(nSpecies)*sizeof(double));
	double *T = (double*)malloc((3)*(nSpecies)*sizeof(double));

  // add External E
	//gZero(E); // for testing Boris
	//gAddTo(Ext); //needs grid definition of Eext
  puAddEext(ini, pop, E); // adds same value to whole grid

  gMul(E, 0.5);
	puGet3DRotationParameters(ini, T, S, 0.5);
	acc(pop, E, T, S);
	gMul(E, 2.0);
	puGet3DRotationParameters(ini, T, S, 1.0);

	/*
	 * TIME LOOP
	 */

	Timer *t = tAlloc(mpiInfo->mpiRank);

	// n should start at 1 since that's the timestep we have after the first
	// iteration (i.e. when storing H5-files).
	int nTimeSteps = iniGetInt(ini,"time:nTimeSteps");
	for(int n = 1; n <= nTimeSteps; n++){


		msg(STATUS,"Computing time-step %i",n);
        msg(STATUS, "Nr. of particles %i: ",(pop->iStop[0]- pop->iStart[0]));

		MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary

		// Check that no particle moves beyond a cell (mostly for debugging)
		pVelAssertMax(pop,maxVel);

		tStart(t);

		// Move particles
		// oRayTrace(pop, obj, deltaRho); <- do we need this still???
		puMove(pop); //puMove(pop, obj); Do not change functions such that PINC does
    // not work in other run modes!

		// Migrate particles (periodic boundaries)
		extractEmigrants(pop, mpiInfo);
		puMigrate(pop, mpiInfo, rho);

      //add influx of new particles on boundary
      pPurgeGhost(pop, rho);
      pFillGhost(ini,pop,rng,mpiInfo);

		// Check that no particle resides out-of-bounds (just for debugging)
		//pPosAssertInLocalFrame(pop, rho); //gives error with open boundary

        // Collect the charges on the objects.
        oCollectObjectCharge(pop, rhoObj, obj, mpiInfo);    // for capMatrix - objects


		/*
		*   Collisions
		*   Changes velocity component of some particles, not position.
		*/
		collide(ini, pop, mccVars, rng,mpiInfo);


		// Compute charge density
		distr(pop, rho,rho_e,rho_i);
		gHaloOp(addSlice, rho, mpiInfo, FROMHALO);
		gHaloOp(addSlice, rho_e, mpiInfo, FROMHALO);
		gHaloOp(addSlice, rho_i, mpiInfo, FROMHALO);


        // Keep writing Rho here.



        // Add object charge to rho.
        gAddTo(rho, rhoObj);
        //gHaloOp(addSlice, rho, mpiInfo, FROMHALO);
		//gAssertNeutralGrid(rho, mpiInfo);

        //gBnd(phi, mpiInfo);
        solve(solver, rho, phi, mpiInfo);                   // for capMatrix - objects




        // Second run with solver to account for charges
        oApplyCapacitanceMatrix(rho, phi, obj, mpiInfo);    // for capMatrix - objects


    //gBnd(phi, mpiInfo);
		solve(solver, rho, phi, mpiInfo);

		//gHaloOp(setSlice, phi, mpiInfo, TOHALO); // Needed by sSolve but not mgSolve

		// Compute E-field
		gFinDiff1st(phi, E);
		gHaloOp(setSlice, E, mpiInfo, TOHALO);
		gMul(E, -1.);

		//gAssertNeutralGrid(E, mpiInfo);
		// Apply external E
    //gZero(E);
    //gAddTo(Ext); //needs grid definition of Eext
    puAddEext(ini, pop, E); // adds same value to whole grid

		// Accelerate particle and compute kinetic energy for step n
		//acc(pop, E);
    acc(pop, E, T, S);

		tStop(t);

		// Sum energy for all species
		pSumKinEnergy(pop);

		// Compute potential energy for step n
		gPotEnergy(rho,phi,pop);

		// Example of writing another dataset to history.xy.h5
		// xyWrite(history,"/group/group/dataset",(double)n,value,MPI_SUM);

		if(n%1000 == 0 || n>29900){//50614
		//Write h5 files
		//gWriteH5(E, mpiInfo, (double) n);
			gWriteH5(rho, mpiInfo, (double) n);
			gWriteH5(rho_e, mpiInfo, (double) n);
			gWriteH5(rho_i, mpiInfo, (double) n);

			gWriteH5(phi, mpiInfo, (double) n);
			//pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5);
			gWriteH5(rhoObj, mpiInfo, (double) n);
		}

		pWriteEnergy(history,pop,(double)n,units);
	}

	if(mpiInfo->mpiRank==0) {
    tMsg(t->total, "Time spent: ");
}

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
    gCloseH5(rhoObj);       // for capMatrix - objects
    oCloseH5(obj);          // for capMatrix - objects
    // 11.10.19 segfault seems to link to oClose(), as calling this
    // alters the segfault.

	xyCloseH5(history);

  // Free memory
  // sFree(solver);
  // mgFreeSolver(solver);
  solverFree(solver);
  mccFreeVars(mccVars);
  gFree(rho);
  gFree(rho_e);
  gFree(rho_i);
  gFree(phi);
  free(S);
  free(T);

  gFree(E);
  pFree(pop);
  uFree(units);
    gFree(rhoObj);          // for capMatrix - objects
    oFree(obj);             // for capMatrix - objects




	gsl_rng_free(rngSync);
	gsl_rng_free(rng);


}
