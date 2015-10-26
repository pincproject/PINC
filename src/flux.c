/*DiP3D*/
/*Functions connected to flux calcul.*/
/*Author: Wojciech Jacek Miloch*/
/*University of Oslo, Norway*/
/*2009*/
/*Last revised 05.03.08*/

#include <math.h>
#include "const.h"
/**********Returns the complementary error function erfc(x)********/
/**********Fract. error everywhere less than 1.2 x 10 -7*****/
/***** Chebyshev fitting ***********/
/********From Numerical Recipes in C, page 221 ***********/
double erfcc(double x)
{
  double t,z,ans;  
  z=fabs(x);
  t=1.0/(1.0+0.5*z);
  ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));    
  return x >= 0.0 ? ans : 2.0-ans;
}

double cumfprim(double v0, int i) //distribution function for the flux
{
  return v0*exp(-((v0-vdriftx[i])*(v0-vdriftx[i]))/(2*vthx[i]*vthx[i]));
}
double zet1(int i) //statistical sum - partition function
{
  double PI=2*asin(1);
  double lb=-vdriftx[i]/(sqrt(2)*vthx[i]);
  return vthx[i]*vthx[i]*exp(-lb*lb)+(sqrt(2*PI)*vdriftx[i]*vthx[i]*erfcc(lb))/2;
}

double cumf(double v0, int i) //cumulative distr. function
{  
  double first;
  double second;
  double PI=2*asin(1);
  double lbound=-vdriftx[i]/(sqrt(2)*vthx[i]);
  double ub=(v0-vdriftx[i])/(sqrt(2)*vthx[i]);
  
  first=vthx[i]*vthx[i]*(exp(-lbound*lbound)-exp(-ub*ub));
  second=0.5*vdriftx[i]*vthx[i]*sqrt(2*PI)*(erfcc(lbound)-erfcc(ub));
  return first+second; 
}
void init_newpart() //initialize parameters for particle injection
{
  int i;
  for(i=0; i<S; i++) 
    {
      llb[i]=-vdriftx[i]/(sqrt(2)*vthx[i]);
      //calculate partition function - statistical sum
      lzet[i] = vthx[i]*vthx[i]*exp(-llb[i]*llb[i])+(sqrt(2*M_PI)*vdriftx[i]*vthx[i]*erfcc(llb[i]))/2;
      //calculate initial velocity for iterative procedure
      //it is v for which distr. function has maximum
      //to avoid instabilities
      lv0[i] =(sqrt(vdriftx[i]*vdriftx[i]+4*vthx[i]*vthx[i])+vdriftx[i])/2;

      //calculate right drift and parameters for the right side
      rvdriftx[i]=-vdriftx[i];
      rlb[i]=-rvdriftx[i]/(sqrt(2)*vthx[i]);
      //calculate partition function - statistical sum
      rzet[i] = vthx[i]*vthx[i]*exp(-rlb[i]*rlb[i])+(sqrt(2*M_PI)*rvdriftx[i]*vthx[i]*erfcc(rlb[i]))/2;
      //calculate initial velocity for iterative procedure
      rv0[i] =(sqrt(rvdriftx[i]*rvdriftx[i]+4*vthx[i]*vthx[i])+rvdriftx[i])/2;
    }
}

void calculate_flux(void)
{
  int i;
  double arger;  
  
  for(i=0; i<S; i++)
    {
		printf("flux %d %E %E\n",i, vdriftx[i], vthx[i]);
      arger=-(vdriftx[i]/(sqrt(2)*vthx[i]));
#ifdef MPI
      //flux per processor, no part on rank 0
      flux[i][0]=(dt*(dens[i])*Ly*Lz*((vthx[i]/sqrt_twopi)*exp(-(vdriftx[i]*vdriftx[i]/(2*vthx[i]*vthx[i]))) + vdriftx[i]*erfcc(arger)/2))/(numtasks-1); //left
      flux[i][1]=(dt*(dens[i])*Ly*Lz*((vthx[i]/sqrt_twopi)*exp(-(vdriftx[i]*vdriftx[i]/(2*vthx[i]*vthx[i]))) - vdriftx[i]*erfcc(arger)/2))/(numtasks-1); //right
      flux[i][2]=(dt*(dens[i])*Lx*Lz*(vthy[i]/sqrt_twopi))/(numtasks-1); //y up
      flux[i][3]=(dt*(dens[i])*Lx*Lz*(vthy[i]/sqrt_twopi))/(numtasks-1); //y down
	  flux[i][4]=(dt*(dens[i])*Lx*Ly*(vthz[i]/sqrt_twopi))/(numtasks-1); //bottom
      flux[i][5]=(dt*(dens[i])*Lx*Ly*(vthz[i]/sqrt_twopi))/(numtasks-1); //top
      //total flux per process
      totalflux[i]=(flux[i][0]+flux[i][1]+flux[i][2]+flux[i][3]+flux[i][4]+flux[i][5]);
      if(rank==0)
	  printf("SPECIE %d: Flux per process \n%f, %f, %f, %f, %f, %f\n",i, flux[i][0], flux[i][1], flux[i][2], flux[i][3], flux[i][4], flux[i][5]);
#else
      flux[i][0]=(dt*(dens[i])*Ly*Lz*((vthx[i]/sqrt_twopi)*exp(-(vdriftx[i]*vdriftx[i]/(2*vthx[i]*vthx[i]))) + vdriftx[i]*erfcc(arger)/2)); //left
      flux[i][1]=(dt*(dens[i])*Ly*Lz*((vthx[i]/sqrt_twopi)*exp(-(vdriftx[i]*vdriftx[i]/(2*vthx[i]*vthx[i]))) - vdriftx[i]*erfcc(arger)/2)); //right
      flux[i][2]=(dt*(dens[i])*Lx*Lz*(vthy[i]/sqrt_twopi)); //yup
      flux[i][3]=(dt*(dens[i])*Lx*Lz*(vthy[i]/sqrt_twopi)); //ydown
	  flux[i][4]=(dt*(dens[i])*Lx*Ly*(vthz[i]/sqrt_twopi)); //bottom
      flux[i][5]=(dt*(dens[i])*Lx*Ly*(vthz[i]/sqrt_twopi)); //top    
      //total flux per process
      totalflux[i]=(flux[i][0]+flux[i][1]+flux[i][2]+flux[i][3]+flux[i][4]+flux[i][5]);
      if(rank==0)
	     printf("SPECIE %d: Flux: %f, %f, %f, %f, %f, %f\n", i, flux[i][0], flux[i][1], flux[i][2], flux[i][3], flux[i][4], flux[i][5]);
	//	arger=arger;
		//arger=0;
	//	printf("check arger %E %E\n",  arger, erfcc(arger));

#endif
      //CHECK IF FLUX IS CORRECT
      
      //calculate mean velocity
      vmean[i]=((vthx[i]/sqrt_twopi)*exp(-(vdriftx[i]*vdriftx[i]/(2*vthx[i]*vthx[i]))) - vdriftx[i]*erfcc(arger)/2);
      //get the fractional part
      //here it is working fine but on ohter computer must be checked for negative values
	  fluxrest[i][0]=flux[i][0]-(int)flux[i][0];  
      fluxrest[i][1]=flux[i][1]-(int)flux[i][1];
      fluxrest[i][2]=flux[i][2]-(int)flux[i][2];  
      fluxrest[i][3]=flux[i][3]-(int)flux[i][3];   	
	  fluxrest[i][4]=flux[i][4]-(int)flux[i][4];  
      fluxrest[i][5]=flux[i][5]-(int)flux[i][5];   
      //printf("fluxrest \n%1.25f, %1.25f, %1.25f, %1.55f\n", fluxrest[i][1], fluxrest[i][2], fluxrest[i][3], fluxrest[i][0]);
      // getchar();
    }
}
