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

/******************************************************************************
* 				Local functions
*****************************************************************************/

int mccTest(int one, int two){
	//bogus test to check if function calls and includes work
	return one+two;
}

/******************************************************************************
* 				Global functions
*****************************************************************************/


/******************************************************************************
* 				Local functions
*****************************************************************************/

static void mccSanity(dictionary *ini, const char* name, int nSpecies){
	// check error, see vahedi + surendra p 181
	// check for v_ion + v_neutral > maxVel (moves beyond a cell in one timestep)
	int countSpecies = iniGetInt(ini,"population:nSpecies");
	if(countSpecies!=nSpecies){
		msg(ERROR,"%s only supports 2 species",name);
	}
}

void mccReadcrossect(const char* filename){
	// read crossections from file?
	FILE* fp = fopen(filename, "r" );

	fclose(fp);
}
double mccGetMaxVel(const Population *pop, int species){

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
double mccGetMinVel(const Population *pop, int species){
	//degug function
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

double mccSigmaCEX(double eps, double DebyeLength){
	// need a function for calculating real cross sects from data...

	double a = 0.7*pow(10,-1); //minimum value of sigma
	double b = 0.5*pow(10,-2); //slew? rate
	double temp = a+(b/(sqrt(eps)));
	//temp = temp/(DebyeLength*DebyeLength); // convert to PINC dimensions
	return temp;
}
double mccSigmaIonElastic(double eps, double DebyeLength){
	// need a function for calculating real cross sects from data...

	double a = 0.3*pow(10,-1); //minimum value of sigma
	double b = 0.5*pow(10,-2); //slew? rate
	double temp = a+(b/(sqrt(eps)));
	//temp = temp/(DebyeLength*DebyeLength); // convert to PINC dimensions
	return temp;
}
double mccSigmaElectronElastic(double eps, double DebyeLength){
	// need a function for calculating real cross sects from data...
	double temp = 0.3*pow(10,-1); // hardcoded order of---
	//temp = temp/(DebyeLength*DebyeLength); // convert to PINC dimensions
	return temp;
}
void mccGetPmaxIon(double m, double thermalVel, double dt, double nt, double DebyeLength,double *Pmax, double *maxfreq, Population *pop){
	// determine max local number desity / max_x(n_t(x_i)) (for target species, neutrals)
	// Get å functional form of sigma_T (total crossect as funct. of energy)
	// determine the speed, v(eps_i) = sqrt(2*eps_i *m_s)
	// determine, max_eps(sigma_T *v)

	// for now lets do
	double max_v = mccGetMaxVel(pop,1); //2.71828*thermalVel; // e*thermalVel, needs to be max_velocity function
	double min_v = mccGetMinVel(pop,1);
	msg(STATUS,"maxVelocity Ion =  %f minVelocity Ion =  %f", max_v, min_v);
	*maxfreq = 0; //start at zero
	double eps = 0;

	double *vel = pop->vel;
	double NewFreq = 0;
	double NewVelocity = 0;
	int nDims = pop->nDims;

	long int iStart = pop->iStart[1]; //ions as specie 1
	long int iStop  = pop->iStop[1];
	for(int i=iStart; i<iStop; i++){
		NewVelocity = sqrt(vel[i*nDims]*vel[i*nDims]+vel[i*nDims+1]\
			*vel[i*nDims+1]+vel[i*nDims+2]*vel[i*nDims+2]);
		eps = 0.5*NewVelocity*NewVelocity*pop->mass[1];
		NewFreq = (mccSigmaCEX(eps, DebyeLength)+mccSigmaIonElastic(eps, DebyeLength))*NewVelocity*nt;
		if(NewFreq>*maxfreq){
			*maxfreq=NewFreq;
		}
	}

	//*maxfreq = (mccSigmaCEX(eps, DebyeLength)+mccSigmaIonElastic(eps, DebyeLength))*max_v*nt; //maxfreq = max_eps(sigma_T *v)*nt (or max frequency of colls?)
	msg(STATUS,"maxfreq Ion =  %f", *maxfreq );
	*Pmax = 1-exp(-(*maxfreq*dt));
	//*Pmax = 0.01;
	msg(STATUS,"getPmax Ion =  %f", *Pmax);
	//return Pmax, maxfreq;
}
void mccGetPmaxElectron(double m, double thermalVel, double dt, double nt, double DebyeLength,double *Pmax, double *maxfreq, Population *pop){
	// determine max local number desity / max_x(n_t(x_i)) (for target species, neutrals)
	// Get å functional form of sigma_T (total crossect as funct. of energy)
	// determine the speed, v(eps_i) = sqrt(2*eps_i *m_s)
	// determine, max_eps(sigma_T *v)

	// for now lets do
	double max_v = mccGetMaxVel(pop,0);//2.71828*thermalVel; // e*thermalVel, needs to be max_velocity
	double min_v = mccGetMinVel(pop,0);
	msg(STATUS,"maxVelocity Electron =  %f minVelocity Electron =  %f", max_v, min_v);
	double eps = 0.0;//0.5*max_v*max_v*pop->mass[0];

	double *vel = pop->vel;
	double NewFreq = 0;
	double NewVelocity = 0;
	int nDims = pop->nDims;

	long int iStart = pop->iStart[0]; //ions as specie 1
	long int iStop  = pop->iStop[0];
	for(int i=iStart; i<iStop; i++){
		NewVelocity = sqrt(vel[i*nDims]*vel[i*nDims]+vel[i*nDims+1]\
			*vel[i*nDims+1]+vel[i*nDims+2]*vel[i*nDims+2]);
		eps = 0.5*NewVelocity*NewVelocity*pop->mass[0];
		NewFreq = mccSigmaElectronElastic(eps, DebyeLength)*NewVelocity*nt;
		if(NewFreq>*maxfreq){
			*maxfreq=NewFreq;
		}
	}

	//*maxfreq = mccSigmaElectronElastic(eps, DebyeLength)*max_v*nt; //maxfreq = max_eps(sigma_T *v)*nt (or max frequency of colls?)
	msg(STATUS,"maxfreq electron =  %f", *maxfreq );
	*Pmax = 1-exp(-(*maxfreq*dt));
	//*Pmax = 0.01;
	msg(STATUS,"getPmax Electron =  %f", *Pmax);
	//return Pmax, maxfreq;
}
double mccGetMyCollFreq(double (*sigma)(double, double), double m, double vx,double vy,double vz, double dt, double nt, double DebyeLength){
	// collision freq given a cross section function
	double v = sqrt(vx*vx + vy*vy + vz*vz);
	double eps_i = 0.5*m*v*v;
	double sigma_T = sigma(eps_i, DebyeLength);
	double MyFreq = v*sigma_T*nt;
	//msg(STATUS,"Myfreq =  %f", MyFreq);
	return MyFreq;
}

double mccEnergyDiffIonElastic(double Ekin, double theta, double mass1, double mass2){
	double temp;
	//msg(STATUS, "Ekin = %.32f", Ekin);
	//msg(STATUS, "theta = %f", theta);
	//msg(STATUS, "mass1 = %f", mass1);
	//msg(STATUS, "mass2 = %f", mass2);
	temp = Ekin*( 1-((2*mass1*mass2)/((mass1+mass2)*(mass1+mass2)) )*(1.0-cos(theta)) );
	//temp = Ekin*(cos(2*theta)*cos(2*theta));
	//msg(STATUS, "temp = %.32f", temp);
	return temp;
}

void mccCollideElectron(Population *pop,  double *Pmax, double *maxfreqElectron, const gsl_rng *rng, double dt,double nt, double DebyeLength){

	msg(STATUS,"colliding Electrons");

	int nDims = pop->nDims;
	double *vel = pop->vel;
	double *mass = pop->mass;

	double R = gsl_rng_uniform_pos(rng);
	double Rp = 1.0;//gsl_rng_uniform(rng);

	double angleChi = 0; //scattering angle in x, y
	double anglePhi = 0; //scattering angle in j
	double angleTheta = 0;
	double A = 0; //scaling Factor used for precomputing (speedup)
	long int q = 0;
	double vx, vy, vz;

	long int errorcounter = 0;
	long int errorcounter3 = 0;
	long int debugindex = 0;
	double Ekin = 0;
	long int iStart = pop->iStart[0];		// *nDims??
	long int iStop = pop->iStop[0];      // make shure theese are actually electrons!
	long int NparticleColl = (*Pmax)*(pop->iStop[0]-pop->iStart[0]); // 1 dim ?, Number of particles too collide
	long int mccStepSize = floor((double)((iStop-iStart))\
	/ (double)(NparticleColl));
	//if ((mccStepSize-1)*(NparticleColl) > (iStop-iStart)){ //errorcheck, remove
	//	msg(WARNING,"particle collisions out of bounds in mccCollideElectron");
	//}

	msg(STATUS,"colliding %i of %i electrons", NparticleColl, (iStop-iStart));
	msg(STATUS,"Particles per box to pick one from = %i", (mccStepSize));

	long int mccStop = iStart + mccStepSize*NparticleColl; //enshure non out of bounds
	if (mccStop > iStop){ //errorcheck, remove
		msg(WARNING,"particle collisions out of bounds in mccCollideElectron");
		msg(WARNING,"mccStop = %i is bigger then iStop = %i",mccStop,iStop);
	}
	for(long int i=iStart;i<mccStop-nDims;i+=mccStepSize){  //check this statement #segfault long int p=pStart;p<pStart+Pmax*(pStop-pStart);p++)
		errorcounter3 += 1;
		//msg(STATUS,"errorcounter = %i", (errorcounter));
		if (errorcounter3 > NparticleColl){ //errorcheck, remove
			msg(ERROR,"errorcounter = %i should be the same or less as the \
			number of coll	particles = %i", (errorcounter),NparticleColl/nDims);
		//	msg(WARNING,"particle collisions out of bounds in mccCollideIon");
			//errorcounter = 0;
		}
		if (i > iStop){
			msg(STATUS,"iStop = %i and index=%i", iStop, i);
		}
		if (i < iStart){
			msg(STATUS,"iStop = %i and index=%i", iStop, i);
		}
		R = gsl_rng_uniform_pos(rng); // New random number per particle. maybe print?
		Rp = gsl_rng_uniform_pos(rng); // separate rand num. for prob.
		q = floor((i + (R*mccStepSize))*nDims);
		//q = i*nDims; //testing with first index
		if(q>iStop*nDims){
			msg(ERROR,"index out of range, q>iStop in collideElectron");
		}
		//msg(STATUS,"R=%f, Rp = %f, q = %i iStart = %i, iStop = %i", R,Rp,q,iStart,iStop );
		//if(q>iStart+5)
		//	msg(ERROR, "qsss");

		vx = vel[q];
		vy = vel[q+1];
		vz = vel[q+2];
		// need n,n+1,n+2 for x,y,j ?
		//msg(STATUS,"vx = %f, vy = %f, vz = %f", vx,vy,vz);

		double MyCollFreq = mccGetMyCollFreq(mccSigmaElectronElastic, mass[0],vx,vy,vz,dt,nt,DebyeLength); //prob of coll for particle i

		//msg(STATUS,"is Rp = %f < MyCollFreq/maxfreqElectron = %f", Rp, (MyCollFreq/ *maxfreqElectron));
		//msg(STATUS,"MyCollFreq %f, maxfreqElectron = %f", MyCollFreq, *maxfreqElectron);

		if (Rp<(MyCollFreq/ *maxfreqElectron)){
		// Pmax gives the max possible colls. but we will "never" reach this
		// number, so we need to test for each particles probability.
			errorcounter += 1;
			debugindex = i;
			if (( (MyCollFreq)/ *maxfreqElectron)>1.0000000001){
				//put in sanity later, need real cross sects
				msg(WARNING,"MyCollFreqElectron/ *maxfreqElectron > 1, MyCollFreq = %f, maxFreq = %f",MyCollFreq, *maxfreqElectron);
			}
			R = gsl_rng_uniform_pos(rng); // New random number per particle. maybe print?

			Ekin = 0.5*(vx*vx + vy*vy + vz*vz)*mass[0]; // values stored in population for speed??
			//msg(STATUS,"ekin = %f, R =  %f", Ekin, R);
			double argument = (2+Ekin-2*pow((1+Ekin),R))/(Ekin); // This one has the posibility to be greater than 1!!

			if(sqrt(argument*argument)>1.0){
				double newargument = sqrt(argument*argument)/argument;
				msg(WARNING,"old argument = %.32f new arg = %.64f", argument, newargument);
				argument = newargument;
			}
			angleChi =  acos(argument); // gives nan value if abs(argument) > 1

			//for debugging
			if(sqrt(argument*argument)>1.0){
				msg(ERROR,"argument of cos() is greater than 1. argument = %32f", argument);
			}
			if(isnan(angleChi)){
				msg(STATUS, "abs(argument) = %32f", sqrt(argument*argument));
				msg(ERROR,"angleChi == nan!!!! new velocity is cos(nan)  argument of acos() is %.32f!!", argument);
			}
			//msg(STATUS,"angleChi = %f, power is %.32f, argument is %f",angleChi, pow((1+Ekin),R),argument);
			anglePhi = 2*PI*R; // set different R?
			angleTheta = acos(vx);
			if(sin(angleTheta) == 0.0){
				msg(ERROR,"float division by zero in collideion", argument);
			}

			A = (sin(angleChi)*sin(anglePhi))/sin(angleTheta);
			vel[q] = vx*cos(angleChi)+A*(vy*vy + vz*vz); //Vx
			if(isnan(angleTheta)){
				msg(STATUS, "abs(argument) = %32f", sqrt(argument*argument));
				msg(ERROR,"angleTheta == nan!!!! new velocity is nan*(somestuffs)  = %.32f!!", vx*cos(angleChi)+A*(vy*vy + vz*vz));
			}
			if(isinf(A)){
				msg(ERROR,"A is inf in Collide Electron !!");
			}
			if(isnan(A)){
				msg(ERROR,"A is nan in Collide Electron !!");
			}
			vel[q+1] = vy*cos(angleChi)+A*vz-A*vx*vy; //Vy
			vel[q+2] = vz*cos(angleChi)-A*vy-A*vx*vz; //Vz
			//msg(STATUS,"newvx = %f, newvy = %f, newvz = %f", vel[q],vel[q+1],vel[q+2]);
		}

	}
	//msg(STATUS,"errorcounter = %i should be the same as number of\
	coll particles = %i", (errorcounter),NparticleColl);
	msg(STATUS,"%i Electron collisions actually performed and last index was %i", errorcounter, debugindex);
	msg(STATUS,"Done colliding Electrons");
	//free?
}
 //1590430

void mccCollideIon(const dictionary *ini, Population *pop, double *Pmax,
	double *maxfreqIon, const gsl_rng *rng, double timeStep, double nt,
	double DebyeLength){

	msg(STATUS,"colliding Ions");

	//int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *vel = pop->vel;
	double *mass = pop->mass;
	double NvelThermal = iniGetDouble(ini,"collisions:thermalVelocityNeutrals");
	//double timeStep = iniGetDouble(ini,"time:timeStep");
	double stepSize = iniGetDouble(ini,"grid:stepSize");
	double velTh = (timeStep/stepSize)*NvelThermal; //neutrals
	//free(velThermal);

	double R, R1, Rp;

	double angleChi = 0; //scattering angle in x, y
	double anglePhi = 0; //scattering angle in j
	double angleTheta = 0;
	double A = 0; //scaling Factor used for precomputing (speedup)
	long int q = 0;
	//double fraction = 0.7; //defines fraction of particles colliding w el/ch-ex
	double vxMW, vyMW, vzMW;
	long int errorcounter = 0;
	long int errorcounter1 = 0;
	long int errorcounter2 = 0;
	long int errorcounter3 = 0;
	long int errorcounter4 = 0;
	long int debugindex = 0;
	double Ekin = 0;
	double EkinAfter = 0;
	double EnergyDiff = 1;
	double AccumulatedEnergyDiff = 0;

	long int iStart = pop->iStart[1];			// *nDims??
	long int iStop = pop->iStop[1];      // make shure theese are actually ions!
	long int NparticleColl = (*Pmax)*(pop->iStop[1]-pop->iStart[1]); // 1 dim ?, Number of particles too collide
	long int mccStepSize = floor((double)((iStop-iStart))\
	/ (double)(NparticleColl));
	//if ((mccStepSize-1)*(NparticleColl)/nDims > (iStop-iStart)){ //errorcheck, remove
	//	msg(WARNING,"particle collisions out of bounds in mccCollideIon");
	//	msg(WARNING,"%i is bigger than array = %i",(mccStepSize-1)\
	//	*(NparticleColl),(iStop-iStart) );
	//}

	msg(STATUS,"colliding %i of %i ions", NparticleColl, (iStop-iStart));
	msg(STATUS,"Particles per box to pick one from = %i", (mccStepSize));

	long int mccStop = iStart + mccStepSize*NparticleColl; //enshure non out of bounds
	// but also that last index probably never collides...
	if (mccStop > iStop){ //errorcheck, remove
		msg(WARNING,"particle collisions out of bounds in mccCollideIon");
		msg(WARNING,"mccStop = %i is bigger then iStop = %i",mccStop,iStop);
	}
	for(long int i=iStart;i<mccStop;i+=mccStepSize){  //check this statement #segfault long int p=pStart;p<pStart+Pmax*(pStop-pStart);p++)
		errorcounter3 += 1;
		if (errorcounter3 > NparticleColl){ //errorcheck, remove
			msg(ERROR,"errorcounter = %i should be the same or less as the \
			number of coll	particles = %i", (errorcounter),NparticleColl/nDims);
		//	msg(WARNING,"particle collisions out of bounds in mccCollideIon");
			//errorcounter = 0;
		}
		if (i > iStop){
			msg(STATUS,"iStop = %i and index=%i", iStop, i);
		}
		if (i < iStart){
			msg(STATUS,"iStart = %i and index=%i", iStart, i);
		}
		Rp = gsl_rng_uniform_pos(rng); //decides type of coll.
		R = gsl_rng_uniform_pos(rng); // New random number per particle. maybe print?
		q = (i + (R*mccStepSize))*nDims;
		//msg(STATUS,"R=%g, p = %i, q = %i", R,p,q );
		if(q>iStop*nDims){
			msg(ERROR,"index out of range, q>iStop in collideIon");
		}
		vxMW = gsl_ran_gaussian_ziggurat(rng,velTh); //maxwellian dist?
		vyMW = gsl_ran_gaussian_ziggurat(rng,velTh);
		vzMW = gsl_ran_gaussian_ziggurat(rng,velTh);
		//msg(STATUS,"neutral velx = %f vely = %f velz = %f",vxMW,vyMW,vzMW );
		double vxTran = vel[q]-vxMW;   //simple transfer, should be picked from maxwellian
		double vyTran = vel[q+1]-vyMW;  //vel-velMaxwellian
		double vzTran = vel[q+2]-vzMW;  // this will for now break conservation laws???

		double MyCollFreq1 = mccGetMyCollFreq(mccSigmaIonElastic, mass[1],vxTran,vyTran,vzTran,timeStep,nt,DebyeLength); //prob of coll for particle i
		double MyCollFreq2 = mccGetMyCollFreq(mccSigmaCEX, mass[1],vxTran,vyTran,vzTran,timeStep,nt,DebyeLength); //duplicate untill real cross sections are implemented

		//msg(STATUS,"is Rp = %f < MyCollFreq/maxfreqElectron = %f", Rp, ((MyCollFreq1+MyCollFreq2)/ *maxfreqIon));
		if (Rp<( (MyCollFreq1+MyCollFreq2)/ *maxfreqIon)){ // MyCollFreq1+MyCollFreq2 defines total coll freq.
			errorcounter2 += 1;
			debugindex = i;
			if (( (MyCollFreq1+MyCollFreq2)/ *maxfreqIon)>1.0000000001){
				//put in sanity later, need real cross sects
				msg(WARNING,"(MyCollFreq1+MyCollFreq2)/ *maxfreqIon > 1");
				msg(WARNING,"(MyCollFreq1= %f MyCollFreq2 = %f *maxfreqIon =%f",MyCollFreq1,MyCollFreq2,*maxfreqIon);
			}
			if(Rp < MyCollFreq1/ *maxfreqIon){ // if test slow, but enshures randomness...
				// elastic:

				//msg(STATUS,"elastic");


				//msg(ERROR,"Elastic collision"); // temporary test
				Ekin = 0.5*(vel[q]*vel[q] + vel[q+1]*vel[q+1] + vel[q+2]*vel[q+2])*mass[1];
				R1 = gsl_rng_uniform_pos(rng);

				//The scattering happens in center of mass so we use THETA not chi
				double argument = 1-2*R; //chi is THETA

				if(sqrt(argument*argument)>1.0){
					double newargument = sqrt(argument*argument)/argument;
					msg(WARNING,"old argument = %.32f new arg = %.64f", argument, newargument);
					argument = newargument;
				}
				angleChi =  acos(argument); // gives nan value if abs(argument) > 1

				//for debugging
				if(sqrt(argument*argument)>1.0){
					msg(ERROR,"argument of cos() is greater than 1. argument = %32f", argument);
				}
				if(isnan(angleChi)){
					msg(STATUS, "abs(argument) = %32f", sqrt(argument*argument));
					msg(ERROR,"angleChi == nan!!!! new velocity is cos(nan)  argument of acos() is %.32f!!", argument);
				}

				//angleChi = acos(sqrt(1-R)); // for M_ion = M_neutral
				anglePhi = 2*PI*R1; // set different R?
				//angleTheta = 2*angleChi; //NOPE! not the same!
				angleTheta = acos(vxTran); //C_M frame
				EnergyDiff = mccEnergyDiffIonElastic(Ekin, angleChi, mass[1], mass[1]);//Acctually Ion and neutral mass
				//if(sin(angleTheta) == 0.0){
				//	msg(ERROR,"float division by zero in collideion", argument);
				//}
				A = (sin(angleChi)*sin(anglePhi))/sin(angleTheta);

				if(isinf(A)){
					msg(ERROR,"A is inf in Collide ion !!");
				}
				if(isnan(A)){
					msg(ERROR,"A is nan in Collide Ion !!");
				}
				vel[q] = vxTran*cos(angleChi)+A*(vyTran*vyTran + vzTran*vzTran) \
				+ vxMW; //Vx
				vel[q+1] = vyTran*cos(angleChi)+A*vzTran-A*vxTran*vyTran + vyMW; //Vy
				vel[q+2] = vzTran*cos(angleChi)-A*vyTran-A*vxTran*vzTran + vzMW; //Vz
				if(sqrt(vel[q]*vel[q] + vel[q+1]*vel[q+1]+vel[q+2]*vel[q+2])>sqrt(vxTran*vxTran+vyTran*vyTran+vzTran*vzTran)*10){
					msg(STATUS,"argument of acos() = %f", argument);
					msg(STATUS,"acos() = %f",angleChi );
					msg(STATUS,"A = (sin(angleChi)*sin(anglePhi))/sin(angleTheta) = %f",A );
					msg(STATUS,"vx = %f, vy = %f, vz = %f",vel[q],vel[q+1],vel[q+2]);
					errorcounter4 += 1;
				}
				//if(sqrt(A*A)<0.00000001){
				//	msg(STATUS,"argument of acos() = %f", argument);
				//	msg(STATUS,"acos() = %f",angleChi );
				//	msg(STATUS,"A = (sin(angleChi)*sin(anglePhi))/sin(angleTheta) = %f",A );
				//}
				EkinAfter = 0.5*(vel[q]*vel[q] + vel[q+1]*vel[q+1] + vel[q+2]*vel[q+2])*mass[1];
				AccumulatedEnergyDiff += (Ekin-EkinAfter);
				///msg(STATUS," ");
				///msg(STATUS, "---------- energydiff = %.32f", (Ekin-EkinAfter) );
				///msg(STATUS, "analytical energydiff = %.32f",EnergyDiff );
				///msg(STATUS, "deviation = %.32f",EnergyDiff-(Ekin-EkinAfter) );
				//msg(STATUS," ");
				errorcounter1 += 1;
			}else{
				// Charge exchange:
				// flip ion and neutral. Neutral is new ion.
				//msg(STATUS,"ch-ex");
				errorcounter += 1;
				vel[q] = vxMW; //Vx
				vel[q+1] = vyMW; //Vy
				vel[q+2] = vzMW; //Vz
			}
		}
	}
	msg(STATUS, "AccumulatedEnergyDiff = %.32f",AccumulatedEnergyDiff );
	msg(STATUS, "counted %i energy increases",errorcounter4);
	msg(STATUS,"%i ion collisions actually performed and last index was %i", errorcounter2, debugindex-iStart);
	msg(STATUS,"%i ion collisions as ch-ex and %i as elastic, the sum should be %i", errorcounter, errorcounter1,errorcounter2);

	//msg(STATUS,"errorcounter = %i should be the same as number of\
	coll particles = %i", (errorcounter),NparticleColl);

	//free?
	msg(STATUS,"Done colliding Ions");
}



/*************************************************
*		Inline functions
************************************************/



/*************************************************
*		DEFINITIONS
************************************************/



/*************************************************
* 		ALLOCATIONS
* 		DESTRUCTORS
************************************************/


/*******************************************************
*			VARIOUS COMPUTATIONS (RESIDUAL)
******************************************************/



/*************************************************
*		RUNS
************************************************/

funPtr mccTestMode_set(dictionary *ini){
	//test sanity here!
	mccSanity(ini,"mccTestMode",2);
	return mccTestMode;
}

void mccTestMode(dictionary *ini){
	int errorvar = 0;
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

	void (*distr)() = select(ini,"methods:distr",	puDistr3D1_set,
	puDistrND1_set,
	puDistrND0_set);

	void (*solve)() = select(ini,"methods:poisson", mgSolve_set);

	void (*extractEmigrants)() = select(ini,"methods:migrate",	puExtractEmigrants3D_set,
	puExtractEmigrantsND_set);



	/*
	* INITIALIZE PINC VARIABLES
	*
	*done in "main" for "regular mode"
	*/
	MpiInfo *mpiInfo = gAllocMpi(ini);
	Population *pop = pAlloc(ini);
	Grid *E   = gAlloc(ini, VECTOR);
	Grid *rho = gAlloc(ini, SCALAR);
	Grid *res = gAlloc(ini, SCALAR);
	Grid *phi = gAlloc(ini, SCALAR);
	Multigrid *mgRho = mgAlloc(ini, rho);
	Multigrid *mgRes = mgAlloc(ini, res);
	Multigrid *mgPhi = mgAlloc(ini, phi);

	/*
	* mcc specific variables
	*/
	// make struct???
	double PmaxElectron = 0;
	double maxfreqElectron = 0;
	double PmaxIon = 0;
	double maxfreqIon = 0;
	double *mass = pop->mass;
	int nSpecies = pop->nSpecies;
	double nt = iniGetDouble(ini,"collisions:numberDensityNeutrals"); //constant for now
	double DebyeLength = iniGetDouble(ini,"grid:debye"); // for converting crosssection dimensions
	double mccTimestep = iniGetDouble(ini,"time:timeStep");
	//double frequency = iniGetDouble(ini,"collisions:collisionFrequency"); // for static Pmax...
	double *velThermal = iniGetDoubleArr(ini,"population:thermalVelocity",nSpecies);

	// using Boris algo
	double *S = (double*)malloc((3)*(nSpecies)*sizeof(double));
	double *T = (double*)malloc((3)*(nSpecies)*sizeof(double));

	// Creating a neighbourhood in the rho to handle migrants
	gCreateNeighborhood(ini, mpiInfo, rho);

	// Setting Boundary slices
	gSetBndSlices(phi, mpiInfo);

	//Set mgSolve
	MgAlgo mgAlgo = getMgAlgo(ini);

	// Random number seeds
	gsl_rng *rngSync = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,mpiInfo->mpiRank+1); // Seed needs to be >=1

	/*
	* PREPARE FILES FOR WRITING
	*/
	int rank = phi->rank;
	double *denorm = malloc((rank-1)*sizeof(*denorm));
	double *dimen = malloc((rank-1)*sizeof(*dimen));

	for(int d = 1; d < rank;d++) denorm[d-1] = 1.;
	for(int d = 1; d < rank;d++) dimen[d-1] = 1.;

	pOpenH5(ini, pop, "pop");
	//gOpenH5(ini, rho, mpiInfo, denorm, dimen, "rho");
	gOpenH5(ini, phi, mpiInfo, denorm, dimen, "phi");
	//gOpenH5(ini, E,   mpiInfo, denorm, dimen, "E");
	// oOpenH5(ini, obj, mpiInfo, denorm, dimen, "test");
	// oReadH5(obj, mpiInfo);

	hid_t history = xyOpenH5(ini,"history");
	pCreateEnergyDatasets(history,pop);

	// Add more time series to history if you want
	// xyCreateDataset(history,"/group/group/dataset");

	free(denorm);
	free(dimen);

	/*
	* INITIAL CONDITIONS
	*/

	//msg(STATUS, "constants = %f, %f", velThermal[0], velThermal[1]);
	// Initalize particles
	pPosLattice(ini, pop, mpiInfo);
	//pPosUniform(ini, pop, mpiInfo, rngSync);
	//pVelZero(pop);
	pVelConstant(ini, pop, velThermal[0], velThermal[1]); //constant values for vel.
	//pVelMaxwell(ini, pop, rng);
	double maxVel = iniGetDouble(ini,"population:maxVel");

	// Perturb particles
	// pPosPerturb(ini, pop, mpiInfo);

	// compute Pmax once outside loop. More correct is every dt
	mccGetPmaxElectron(mass[0], velThermal[0], mccTimestep, nt, DebyeLength, &PmaxElectron, &maxfreqElectron, pop);
	mccGetPmaxIon(mass[1], velThermal[1], mccTimestep, nt, DebyeLength, &PmaxIon, &maxfreqIon, pop);
	//msg(STATUS, "freq is %f",frequency);

	// Migrate those out-of-bounds due to perturbation
	extractEmigrants(pop, mpiInfo);
	puMigrate(pop, mpiInfo, rho);

	/*
	* compute initial half-step
	*/

	// Get initial charge density
	distr(pop, rho);
	gHaloOp(addSlice, rho, mpiInfo, FROMHALO);

	// Get initial E-field
	solve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);
	gFinDiff1st(phi, E);
	gHaloOp(setSlice, E, mpiInfo, TOHALO);
	gMul(E, -1.);

	// add External E
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
	gWriteH5(phi, mpiInfo, 0.0);
	//pWriteH5(pop, mpiInfo, 0.0, 0.5);
	pWriteEnergy(history,pop,0.0);



	int collsOnOff;
	collsOnOff = iniGetInt(ini,"collisions:collisionsOnOff");
	if(collsOnOff == 0){
		PmaxElectron = 0;
		PmaxIon = 0;
		msg(STATUS, "----COLLISIONS ARE TURNED OFF----");
	}

	//msg(STATUS, "Pmax for Electrons is %f",PmaxElectron);
	//msg(STATUS, "Pmax for Ions is %f",PmaxIon);

	/*
	* TIME LOOP
	*/

	Timer *t = tAlloc(rank);

	// n should start at 1 since that's the timestep we have after the first
	// iteration (i.e. when storing H5-files).
	int nTimeSteps = iniGetInt(ini,"time:nTimeSteps");
	for(int n = 1; n <= nTimeSteps; n++){

		msg(STATUS,"Computing time-step %i",n);
		//MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary

		// Check that no particle moves beyond a cell (mostly for debugging)
		pVelAssertMax(pop,maxVel);

		tStart(t);


		// Move particles
		//msg(STATUS, "moving particles");
		puMove(pop);
		// oRayTrace(pop, obj);

		mccGetPmaxElectron(mass[0], velThermal[0], mccTimestep, nt, DebyeLength, &PmaxElectron, &maxfreqElectron, pop);
		mccGetPmaxIon(mass[1], velThermal[1], mccTimestep, nt, DebyeLength, &PmaxIon, &maxfreqIon, pop);

		//MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary
		mccCollideElectron(pop, &PmaxElectron, &maxfreqElectron, rng, mccTimestep, nt,DebyeLength); //race conditions?????
		//MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary
		mccCollideIon(ini, pop, &PmaxIon, &maxfreqIon, rng, mccTimestep, nt, DebyeLength);
		//MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary

		// Migrate particles (periodic boundaries)
		//msg(STATUS, "extracting emigrants");
		extractEmigrants(pop, mpiInfo);
		//msg(STATUS, "Migrating particles");
		puMigrate(pop, mpiInfo, rho);

		//MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary
		// Check that no particle resides out-of-bounds (just for debugging)
		//msg(STATUS, "checking particles out of bounds");
		pPosAssertInLocalFrame(pop, rho);

		//MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary
		//-----------------------------
		// Compute charge density
		//msg(STATUS, "computing charge density");
		distr(pop, rho);

		//MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary

		//msg(STATUS, "MPI exchanging density");
		gHaloOp(addSlice, rho, mpiInfo, FROMHALO);

		//---------------------------
		//msg(STATUS, "gAssertNeutralGrid");
		gAssertNeutralGrid(rho, mpiInfo);

		// Compute electric potential phi
		//msg(STATUS, "Compute electric potential phi");
		solve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);
		//msg(STATUS, "gAssertNeutralGrid");
		gAssertNeutralGrid(phi, mpiInfo);

		// Compute E-field
		//msg(STATUS, "compute E field");
		gFinDiff1st(phi, E);
		//msg(STATUS, "MPI exchanging E");
		gHaloOp(setSlice, E, mpiInfo, TOHALO);
		//msg(STATUS, "gMul( E, -1)");
		gMul(E, -1.);

		//msg(STATUS, "gAssertNeutralGrid");
		gAssertNeutralGrid(E, mpiInfo);
		// Apply external E
		// gAddTo(Ext);
		gZero(E);
		puAddEext(ini, pop, E);

		// Accelerate particle and compute kinetic energy for step n
		//msg(STATUS, "Accelerate particles");

		acc(pop, E, T, S);
		tStop(t);

		// Sum energy for all species
		//msg(STATUS, "sum kinetic energy");
		pSumKinEnergy(pop);
		//msg(STATUS, "compute pot energy");
		// Compute potential energy for step n
		gPotEnergy(rho,phi,pop);
		if(n>nTimeSteps-1000){
			msg(STATUS, "writing over a given timestep to file");
			// Example of writing another dataset to history.xy.h5
			// xyWrite(history,"/group/group/dataset",(double)n,value,MPI_SUM);
			//Write h5 files
			//gWriteH5(E, mpiInfo, (double) n);
			//gWriteH5(rho, mpiInfo, (double) n);
			gWriteH5(phi, mpiInfo, (double) n);
			//pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5);
		}
		if(n>nTimeSteps-50){
			pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5);
		}
//		if(n%1000==0 && n<48000){
//			msg(STATUS, "writing every 1000 timestep to file");
//			// Example of writing another dataset to history.xy.h5
//			// xyWrite(history,"/group/group/dataset",(double)n,value,MPI_SUM);
//			//Write h5 files
//			//gWriteH5(E, mpiInfo, (double) n);
//			//gWriteH5(rho, mpiInfo, (double) n);
//			gWriteH5(phi, mpiInfo, (double) n);
//			//pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5);
//		}
		//gWriteH5(phi, mpiInfo, (double) n);
		pWriteEnergy(history,pop,(double)n);

		//free(nt);
		//free(mccTimestep);
		//free(frequency);
		//free(velThermal);
		//msg(STATUS, "   -    ");


	}
	//msg(STATUS, "Test returned %d", errorvar);

	if(mpiInfo->mpiRank==0) tMsg(t->total, "Time spent: ");

	/*
	* FINALIZE PINC VARIABLES
	*/
	gFreeMpi(mpiInfo);

	// Close h5 files
	pCloseH5(pop);
	//gCloseH5(rho);
	gCloseH5(phi);
	//gCloseH5(E);
	// oCloseH5(obj);
	xyCloseH5(history);

	// Free memory
	mgFree(mgRho);
	mgFree(mgPhi);
	mgFree(mgRes);
	gFree(rho);
	gFree(phi);
	gFree(res);
	gFree(E);
	pFree(pop);
	free(S);
	free(T);
	// oFree(obj);

	gsl_rng_free(rngSync);
	gsl_rng_free(rng);

}
