 /* DiP3D */
/* Particle movers */
/* Author: Wojciech Jacek Miloch */
/* University of Oslo */
/* 2008 */
#include "const.h"
#include <math.h>
#include <stdlib.h>
void accel(float factor)
{
  int i,ii;
  int j,k,l;
  double qmratio;
  double x,y,z,ax,ay,az;
  double wjkl,wj1kl,wj1k1l,wjk1l,wjkl1,wj1kl1,wj1k1l1,wjk1l1; //weights
  for(i=0;i<S;i++) //for each specie
    {
      qmratio=qm[i];
      for(ii=0;ii<npart[i];ii++) 
	{
	  //find forces from the closest grid points,and locate particle
	  j=spec[i].part[ii].x/dx;
	  k=spec[i].part[ii].y/dy;
	  l=spec[i].part[ii].z/dz;
	  x=spec[i].part[ii].x-j*dx;
	  y=spec[i].part[ii].y-k*dy;	
	  z=spec[i].part[ii].z-l*dz;	
	  /*find weights, Ex,Ey is multipl.with dt/dxdy*/
	  wjkl=(dx-x)*(dy-y)*(dz-z);
	  wj1kl=x*(dy-y)*(dz-z);
	  wj1k1l=x*y*(dz-z);
	  wjk1l=(dx-x)*y*(dz-z);
	  wjkl1=(dx-x)*(dy-y)*z;
	  wj1kl1=x*(dy-y)*z;
	  wj1k1l1=x*y*z;
	  wjk1l1=(dx-x)*y*z;
	  ax=factor*qmratio*(Fs[ix(0,j,k,l)]*wjkl+Fs[ix(0,j+1,k,l)]*wj1kl+Fs[ix(0,j+1,k+1,l)]*wj1k1l+Fs[ix(0,j,k+1,l)]*wjk1l+Fs[ix(0,j,k,l+1)]*wjkl1+Fs[ix(0,j+1,k,l+1)]*wj1kl1+Fs[ix(0,j+1,k+1,l+1)]*wj1k1l1+Fs[ix(0,j,k+1,l+1)]*wjk1l1);
	  ay=factor*qmratio*(Fs[ix(FsEy,j,k,l)]*wjkl+Fs[ix(FsEy,j+1,k,l)]*wj1kl+Fs[ix(FsEy,j+1,k+1,l)]*wj1k1l+Fs[ix(FsEy,j,k+1,l)]*wjk1l+Fs[ix(FsEy,j,k,l+1)]*wjkl1+Fs[ix(FsEy,j+1,k,l+1)]*wj1kl1+Fs[ix(FsEy,j+1,k+1,l+1)]*wj1k1l1+Fs[ix(FsEy,j,k+1,l+1)]*wjk1l1);
	  az=factor*qmratio*(Fs[ix(FsEz,j,k,l)]*wjkl+Fs[ix(FsEz,j+1,k,l)]*wj1kl+Fs[ix(FsEz,j+1,k+1,l)]*wj1k1l+Fs[ix(FsEz,j,k+1,l)]*wjk1l+Fs[ix(FsEz,j,k,l+1)]*wjkl1+Fs[ix(FsEz,j+1,k,l+1)]*wj1kl1+Fs[ix(FsEz,j+1,k+1,l+1)]*wj1k1l1+Fs[ix(FsEz,j,k+1,l+1)]*wjk1l1);
		  
	  /*accelerate!*/
	  spec[i].part[ii].vx+=ax;
	  spec[i].part[ii].vy+=ay;	
	  spec[i].part[ii].vz+=az;
	}
    }
}

/*** accel and move all the particles***/
void move(int t)
{
  FILE *assigned;
  int i,ii;
  int j,k,l;

   double x,y,z,ax,ay,az;
  double wjkl,wj1kl,wj1k1l,wjk1l,wjkl1,wj1kl1,wj1k1l1,wjk1l1; //weights
  
  double qmratio;
  double partxold, partyold, partzold;
  int vert,kk;
  double q;
  int offset;
  int kp1;
  double partvxold, partvyold,partvzold;
  
  int dno;
  double par;
  double licz, mian,odj;
  double t_tomove, cosphi, sinphi, tempx, tempy;
  double quad_a, quad_b, quad_c,quad_del, t_cross, t_cross1, t_cross2; 
  double xhitp, yhitp;

  //variables for suspected hitting
  double timehitmin_suspected;
  int part_suspected, dnohit_suspected;
  double deltat, deltapar; 
  double par_tmin, par_tmin_suspected;
  int dnoseghit,dnoseghit_suspected;
  deltat=1*dt;
  deltapar=1.0;
		
  for(i=0; i<S; i++)
    {
      printf("ACCEL: in move\n");
      offset=KE_off*i;
      lostpart[i]=0;
      qmratio=qm[i];     
      q=chargeandnorm[i];	          
      for(l=0; l<npart[i]; l++)
	{	 	  
	  //	  printf("before weight\t");
	  /*CALCULATE OLD velocity magnitude*/
	  spec[i].part[ii].kenergy=spec[i].part[ii].vx*spec[i].part[ii].vx+spec[i].part[ii].vy*spec[i].part[ii].vy+spec[i].part[ii].vz*spec[i].part[ii].vz;
	  
	  //find forces from the closest grid points,and locate particle:)
	 
	  j=spec[i].part[ii].x/dx;
	  k=spec[i].part[ii].y/dy;
	  l=spec[i].part[ii].z/dz;
	  x=spec[i].part[ii].x-j*dx;
	  y=spec[i].part[ii].y-k*dy;	
	  z=spec[i].part[ii].z-l*dz;	
	  //  printf("part %d no %d and jk: %d %d\n", i, l, j,k);

	  /*find weights, E is multiplied already by dt/dxdy*/
	  wjkl=(dx-x)*(dy-y)*(dz-z);
	  wj1kl=x*(dy-y)*(dz-z);
	  wj1k1l=x*y*(dz-z);
	  wjk1l=(dx-x)*y*(dz-z);
	  wjkl1=(dx-x)*(dy-y)*z;
	  wj1kl1=x*(dy-y)*z;
	  wj1k1l1=x*y*z;
	  wjk1l1=(dx-x)*y*z;
	  ax=qmratio*(Fs[ix(0,j,k,l)]*wjkl+Fs[ix(0,j+1,k,l)]*wj1kl+Fs[ix(0,j+1,k+1,l)]*wj1k1l+Fs[ix(0,j,k+1,l)]*wjk1l+Fs[ix(0,j,k,l+1)]*wjkl1+Fs[ix(0,j+1,k,l+1)]*wj1kl1+Fs[ix(0,j+1,k+1,l+1)]*wj1k1l1+Fs[ix(0,j,k+1,l+1)]*wjk1l1);
	  ay=qmratio*(Fs[ix(FsEy,j,k,l)]*wjkl+Fs[ix(FsEy,j+1,k,l)]*wj1kl+Fs[ix(FsEy,j+1,k+1,l)]*wj1k1l+Fs[ix(FsEy,j,k+1,l)]*wjk1l+Fs[ix(FsEy,j,k,l+1)]*wjkl1+Fs[ix(FsEy,j+1,k,l+1)]*wj1kl1+Fs[ix(FsEy,j+1,k+1,l+1)]*wj1k1l1+Fs[ix(FsEy,j,k+1,l+1)]*wjk1l1);
	  az=qmratio*(Fs[ix(FsEz,j,k,l)]*wjkl+Fs[ix(FsEz,j+1,k,l)]*wj1kl+Fs[ix(FsEz,j+1,k+1,l)]*wj1k1l+Fs[ix(FsEz,j,k+1,l)]*wjk1l+Fs[ix(FsEz,j,k,l+1)]*wjkl1+Fs[ix(FsEz,j+1,k,l+1)]*wj1kl1+Fs[ix(FsEz,j+1,k+1,l+1)]*wj1k1l1+Fs[ix(FsEz,j,k+1,l+1)]*wjk1l1);
	
	  
	 /*accelerate!*/
	  partvxold=spec[i].part[ii].vx;
	  partvyold=spec[i].part[ii].vy;
	  partvzold=spec[i].part[ii].vz;
	  spec[i].part[ii].vx+=ax;
	  spec[i].part[ii].vy+=ay;
	  spec[i].part[ii].vz+=az;
	  
	  /*CALCULATE KINETIC ENERGY only the velocity squared*/
	  spec[i].part[ii].kenergy+=spec[i].part[ii].vx*spec[i].part[ii].vx+spec[i].part[ii].vy*spec[i].part[ii].vy+spec[i].part[ii].vz*spec[i].part[ii].vz;	 
	  /*COLLECT KINETIC ENERGY to grids, need to be divided by dxdy and multiplied by 0.5*mass[ii] and 0.5 for average and norm factor normvel*/
	  KE[ix(offset,j,k,l)]+=spec[i].part[ii].kenergy*wjkl;
	  KE[ix(offset,j+1,k,l)]+=spec[i].part[ii].kenergy*wj1kl;
	  KE[ix(offset,j+1,k+1,l)]+=spec[i].part[ii].kenergy*wj1k1l;
	  KE[ix(offset,j,k+1,l)]+=spec[i].part[ii].kenergy*wjk1l;      
	  KE[ix(offset,j,k,l+1)]+=spec[i].part[ii].kenergy*wjkl1;
	  KE[ix(offset,j+1,k,l+1)]+=spec[i].part[ii].kenergy*wj1kl1;
	  KE[ix(offset,j+1,k+1,l+1)]+=spec[i].part[ii].kenergy*wj1k1l1;
	  KE[ix(offset,j,k+1,l+1)]+=spec[i].part[ii].kenergy*wjk1l1;  
	  //MOVE NOW
	  
	  vert=0;	  
	  partxold=spec[i].part[ii].x;
	  partyold=spec[i].part[ii].y;
	  partzold=spec[i].part[ii].z;
	  //move particle
	  spec[i].part[ii].x+=spec[i].part[ii].vx*dt;  
	  spec[i].part[ii].y+=spec[i].part[ii].vy*dt;
	  spec[i].part[ii].z+=spec[i].part[ii].vz*dt;

	  partnewx=spec[i].part[ii].x;
	  partnewy=spec[i].part[ii].y;
	  partnewz=spec[i].part[ii].z;
	  /*MOOVING DUST*/	  
	  //	  printf("m d \n");
	  /*check all dusts and all corner pairs*/
	  /*initialization*/
	  int dnohit,kkhit;
	  dnohit=dnohit_suspected=kkhit=-1;
	  double timehitmin=-10*dt;
	  timehitmin_suspected=timehitmin=-10*dt;
	  part_suspected=0;
	  par_tmin=par_tmin_suspected=10000; //infinity
	  dnoseghit=dnoseghit_suspected=-1;

	  /*check now*/
	  for(dno=0; dno<noofdusts; dno++)
	    {
	      if(dshape[dno]==0) //spherical dust
		{
		  //
		  //check if particle suspected
		  //  inside=0;
		  // if(initpartcheck(spec[i].part[ii].x, spec[i].part[ii].y, spec[i].part[ii].z, 0.000001*dx)%2 !=0) inside=1;
		  //check if particle close to the sphere
		  dist_centre=(partnewx-dustcxdx[dno])*(partnewx-dustcxdx[dno])+(partnewy-dustcydx[dno])*(partnewy-dustcydx[dno])+(partnewz-dustczdx[dno])*(partnewz-dustczdx[dno]);
		  if(dist_centre <= 4*dradiusdx[dno]*dradiusdx[dno]) //we are quite close and can do other checks, particle is suspected: (we have distance squared)
		    {
		      //calculate if the particle is within the sphere: ||Po+vdt  - Co|| <= dustradius -> find time.
		      //dust is static
		      quad_a=partvxold*partvxold+partvyold*partvyold+partvzold*partvzold;
		      quad_b=2*(partxold-dustcxdx[dno]+partyold-dustcydx[dno]+partzold-dustczdx[dno]);
		      quad_c=(partxold-dustcxdx[dno])*(partxold-dustcxdx[dno])+(partyold-dustcydx[dno])*(partyold-dustcydx[dno])+(partzold-dustczdx[dno])*(partzold-dustczdx[dno]-dradiusdx[dno]*dradiusdx[dno]);
		      t_cross=-2*deltat; //false time
		      
		      if(quad_a == 0) //linear
			if(quad_b != 0){	 
			  t_cross=-quad_c/quad_b;}
			else //quad_a=quad_b=0, quad_c!=0
			  if(quad_c==0)
			    t_cross=0;//0=0
			  else
			    continue;//1=0
		      else //quadratic equation
			{
			  quad_del=quad_b*quad_b-4.0*quad_a*quad_c;		     
			  if(quad_del< 0)
			    t_cross=-2*deltat; //no cross
			  else //find solution
			    {
			      t_cross1=(-quad_b+sqrt(quad_del))/(2.0*quad_a);
			      t_cross2=(-quad_b-sqrt(quad_del))/(2.0*quad_a);
			      //printf("quad_del %E time %E, %E\n", quad_del, t_cross1, t_cross2);
			      //take only positive time
			      double t_minus=0;		
			      if(t_cross1<0 && t_cross2<0)  //nocross
				{ 
				  t_minus = (t_cross1 > t_cross2) ? t_cross1 : t_cross2;
				  if(t_minus>-0.1*dt)
				    t_cross=t_minus; 
				}
			      //take only positive time
			      else
				if(t_cross1>=0 && t_cross2>=0) 
				  {
				    t_cross = (t_cross1 < t_cross2) ? t_cross1 : t_cross2; //take lower time -> 1st hit
				  }
				else //one is positive, one is negative
				  {				
				    t_cross = (t_cross1 > t_cross2) ? t_cross1 : t_cross2; //take positive time     
				    t_minus = (t_cross1 < t_cross2) ? t_cross1 : t_cross2;
				    if(t_minus>-0.1*dt)
				      t_cross=t_minus; 
				  }
			    }
			}		
		      
		      if((t_cross <= dt+deltat) && (t_cross >=-deltat)) 	//we hit the surface & need to find nearest points on the surface
			{
			  partxhit=partxold+partvxold*t_cross;
			  partyhit=partyold+partvyold*t_cross;
			  partzhit=partzold+partvzold*t_cross;
			  dist_sphere2=(partxhit-dpart[dno][0].x)*(partxhit-dpart[dno][0].x)+(partyhit-dpart[dno][0].y)*(partyhit-dpart[dno][0].y)+(partzhit-dpart[dno][0].z)*(partzhit-dpart[dno][0].z);
			  mindist[0]=dist_sphere2;
			  dist_sphere2=(partxhit-dpart[dno][1].x)*(partxhit-dpart[dno][1].x)+(partyhit-dpart[dno][1].y)*(partyhit-dpart[dno][1].y)+(partzhit-dpart[dno][1].z)*(partzhit-dpart[dno][1].z);
			  mindist[1]=dist_sphere2;
			  dist_sphere2=(partxhit-dpart[dno][2].x)*(partxhit-dpart[dno][2].x)+(partyhit-dpart[dno][2].y)*(partyhit-dpart[dno][2].y)+(partzhit-dpart[dno][2].z)*(partzhit-dpart[dno][2].z);
			  mindist[2]=dist_sphere2;
			  list[0]=0;
			  list[1]=1;
			  list[2]=2;
			  //sort the list...
			  do{
			    swapped=0;
			    for(ct=0; ct<2; ct++)
			      if(mindist[ct]>mindist[ct+1])
				{
				  tempdou=mindist[ct+11];
				  mindist[ct+1]=mindist[ct];
				  mindist[ct]=tempdou;
				  tempint=list[ct+1];
				  list[ct+1]=list[ct];
				  list[ct]=tempint;
				  swapped=1;
				}
			    }while(swapped==1);			
			    //search through all the points
			    for(ct=3; ct<dpartlast[dno]; ct++)
			      {
				dist_sphere2=(partxhit-dpart[dno][0].x)*(partxhit-dpart[dno][0].x)+(partyhit-dpart[dno][0].y)*(partyhit-dpart[dno][0].y)+(partzhit-dpart[dno][0].z)*(partzhit-dpart[dno][0].z);
				if(dist_sphere2<mindist[2])
				  if(dist_sphere2<mindist[1]){
				    if(dist_sphere2<mindist[0])
				      {
					list[2]=list[1];
					list[1]=list[0];
					mindist[2]=mindist[1];
					mindist[1]=mindist[0];
					mindist[0]=dist_sphere2;
					list[0]=ct;
				      }
				    else
				      {
					list[2]=list[1];
					mindist[2]=mindist[1];
					mindist[1]=dist_sphere2;
					list[1]=ct;
				      }}
				  else
				    {
				      mindist[2]=dist_sphere2;
				      list[2]=ct;
				    }
			      }
			    //now assign the charge to the nearest points from the list
			    //Area of the triangle *2 from the cross product of two vectors
			    vec_ax=dpart[dno][list[1]].x-dpart[dno][list[0]].x;
			    vec_ay=dpart[dno][list[1]].y-dpart[dno][list[0]].y;
			    vec_az=dpart[dno][list[1]].x-dpart[dno][list[0]].z;
			    vec_bx=dpart[dno][list[2]].x-dpart[dno][list[0]].x;
			    vec_by=dpart[dno][list[2]].y-dpart[dno][list[0]].y;
			    vec_bz=dpart[dno][list[2]].z-dpart[dno][list[0]].z;

			    localtrianglearea=sqrt(vec_ay*vec_bz-vec_az*vec_by)*(vec_ay*vec_bz-vec_az*vec_by)+(vec_ax*vec_bz-vec_az*vec_bz)*(vec_ax*vec_bz-vec_az*vec_bz)+(vec_ax*vec_by-vec_ay*vec_bx)*(vec_ax*vec_by-vec_ay*vec_bx);

			    //calculate partial weigths
			    vec_ax=dpart[dno][list[1]].x-partxhit;
			    vec_ay=dpart[dno][list[1]].y-partyhit;
			    vec_az=dpart[dno][list[1]].x-partzhit;
			    vec_bx=dpart[dno][list[2]].x-partxhit;
			    vec_by=dpart[dno][list[2]].y-partyhit;
			    vec_bz=dpart[dno][list[2]].z-partzhit;
			    weight0=sqrt(vec_ay*vec_bz-vec_az*vec_by)*(vec_ay*vec_bz-vec_az*vec_by)+(vec_ax*vec_bz-vec_az*vec_bz)*(vec_ax*vec_bz-vec_az*vec_bz)+(vec_ax*vec_by-vec_ay*vec_bx)*(vec_ax*vec_by-vec_ay*vec_bx);
			  
			    vec_ax=dpart[dno][list[0]].x-partxhit;
			    vec_ay=dpart[dno][list[0]].y-partyhit;
			    vec_az=dpart[dno][list[0]].x-partzhit;
			    vec_bx=dpart[dno][list[2]].x-partxhit;
			    vec_by=dpart[dno][list[2]].y-partyhit;
			    vec_bz=dpart[dno][list[2]].z-partzhit;
			    weight1=sqrt(vec_ay*vec_bz-vec_az*vec_by)*(vec_ay*vec_bz-vec_az*vec_by)+(vec_ax*vec_bz-vec_az*vec_bz)*(vec_ax*vec_bz-vec_az*vec_bz)+(vec_ax*vec_by-vec_ay*vec_bx)*(vec_ax*vec_by-vec_ay*vec_bx);

			    vec_ax=dpart[dno][list[1]].x-partxhit;
			    vec_ay=dpart[dno][list[1]].y-partyhit;
			    vec_az=dpart[dno][list[1]].x-partzhit;
			    vec_bx=dpart[dno][list[0]].x-partxhit;
			    vec_by=dpart[dno][list[0]].y-partyhit;
			    vec_bz=dpart[dno][list[0]].z-partzhit;
			    weight2=sqrt(vec_ay*vec_bz-vec_az*vec_by)*(vec_ay*vec_bz-vec_az*vec_by)+(vec_ax*vec_bz-vec_az*vec_bz)*(vec_ax*vec_bz-vec_az*vec_bz)+(vec_ax*vec_by-vec_ay*vec_bx)*(vec_ax*vec_by-vec_ay*vec_bx);


			    //assign charge to nearest grid points
			    dpart[dno][list[0]].q+=(weight0/localtrianglearea)*normalcharge[i];
			    dpart[dno][list[1]].q+=(weight1/localtrianglearea)*normalcharge[i];
			    dpart[dno][list[2]].q+=(weight2/localtrianglearea)*normalcharge[i];

			    /*particle lost*/
			    current[dnohit_suspected][i]++; //CORRECT THE CURRENT ARRAY
			    lostpart[i]++;  
			    lostlist[i][lostpart[i]-1]=l;	
			}			  
		    }
		}

	
	   
	      if(dshape[dno]==1) //other dust
		  ;
		  
		for(kk=0; kk<ncorners[dno]; kk++)
		{
		  quad_a=(daa[dno][kk]-(partvxold*daa1y[dno][kk]+partvyold*daa1x[dno][kk]));
		  quad_b=(dbb[dno][kk]-(partxold*daa1y[dno][kk]+partyold*daa1x[dno][kk]+partvxold*dbb1y[dno][kk]+partvyold*dbb1x[dno][kk]));
		  quad_c=(dcc[dno][kk]-(partxold*dcc1y[dno][kk]-partyold*dcc1x[dno][kk]));
		  //	           	printf("k a:%E\t b%E\t %E\n", quad_a, quad_b, quad_c);
		  t_cross=-2*deltat; //false time
		  
		  if(quad_a == 0) //linear
		    if(quad_b != 0){	 
		      t_cross=-quad_c/quad_b;}
		    else //quad_a=quad_b=0, quad_c!=0
		      if(quad_c==0)
			t_cross=0;//0=0
		      else
			continue;//1=0
		  else //quadratic equation
		    {
		      quad_del=quad_b*quad_b-4.0*quad_a*quad_c;
		      
		      if(quad_del< 0)
			t_cross=-2*deltat; //no cross
		      else //find solution
			{
			  t_cross1=(-quad_b+sqrt(quad_del))/(2.0*quad_a);
			  t_cross2=(-quad_b-sqrt(quad_del))/(2.0*quad_a);
			  //printf("quad_del %E time %E, %E\n", quad_del, t_cross1, t_cross2);
			  //take only positive time
			  double t_minus=0;		
			  if(t_cross1<0 && t_cross2<0)  //nocross
			    { 
			      t_minus = (t_cross1 > t_cross2) ? t_cross1 : t_cross2;
			      if(t_minus>-0.1*dt)
				t_cross=t_minus; 
			    }
			  //take only positive time
			  else
			    if(t_cross1>=0 && t_cross2>=0) 
			      {
				t_cross = (t_cross1 < t_cross2) ? t_cross1 : t_cross2; //take lower time -> 1st hit
			      }
			    else //one is positive, one is negative
			      {
				
				t_cross = (t_cross1 > t_cross2) ? t_cross1 : t_cross2; //take positive time     
				t_minus = (t_cross1 < t_cross2) ? t_cross1 : t_cross2;
				if(t_minus>-0.1*dt)
				  t_cross=t_minus; 
			      }
			}
		    }		  
		  //	printf("p %d time %E kk %d dno %d\n",  l, t_cross, kk, dno);
		  
       
		

		  if((t_cross <= dt+deltat) && (t_cross >=-deltat)) 
		    /* the particle is suspected to hit the dust */
		    /* now lets find the hitting parameter */
		    {		    
		      kp1=kk+1;
		      if(kk==(ncorners[dno]-1))
			kp1=0;
		      licz=mian=par=0;
		      /*I am finding the parameter p here*/
		      if((daa1y[dno][kk]*t_cross+dbb1y[dno][kk])==0);
		    
		      if(fabs(-daa1x[dno][kk]*t_cross-dbb1x[dno][kk]) > fabs(daa1y[dno][kk]*t_cross+dbb1y[dno][kk]))
			{
			  odj=dustxdxold[dno][kk]+dustvx[dno][kk]*t_cross;
			  licz=partxold+partvxold*t_cross-odj;
			  //mian=dustxdxold[dno][kp1]+dustvx[dno][kp1]*t_cross-odj;
			  /*this is doing the same*/
			  mian=	-daa1x[dno][kk]*t_cross-dbb1x[dno][kk];
			}
		      else
			{//only parameter for y
			  
			  odj=dustydyold[dno][kk]+dustvy[dno][kk]*t_cross;
			  licz=partyold+partvyold*t_cross-odj;
			  //  mian=dustydyold[dno][kp1]+dustvy[dno][kp1]*t_cross-odj;
			  
			  //  if(mian != daa1y[dno][kk]*t_cross+dbb1y[dno][kk]) {
			    //printf("t_cross %E and dt wall %d, %E wops %2.15E =?= %2.15E\n",t_cross,dt, kk, mian,daa1y[dno][kk]*t_cross+dbb1y[dno][kk]); 
			    
			    //printf("ncornrers %d, %d %d, %E %E %E %E\n",i, kk, kp1, dustydyold[dno][kk], dustvy[dno][kk], dustydyold[dno][kk+1], dustvy[dno][kk+1]);
			    //getchar();
			  //  }
			  mian=daa1y[dno][kk]*t_cross+dbb1y[dno][kk];
			  
			  //SHORT SAFETY TEST
			  if(mian==0) {			  
			    printf("sth wrong in dust crossing kk %d\nODJ %E LICZ %E MIAN+ODJ %E", kk, odj, licz,dustxdxold[dno][kp1]+dustvx[dno][kp1]*t_cross ); 
			    xhitp=partxold+partvxold*t_cross;
			    yhitp=partyold+partvyold*t_cross;
			    printf("XHIT and YHIT %E\t%E\n", xhitp*normx, yhitp*normx); //getchar();	  
			    printf(" a %E\t b%E\t,c%E\n",  quad_a, quad_b, quad_c);
			    printf("time for crossing: %E\n", t_cross);	
			    printf("kk %d daa %E\t dbb%E\t,dcc%E\n\n\ndcc1x,%E\t dcc1y %E\n", kk, daa[dno][kk], dbb[dno][kk], dcc[dno][kk], dcc1x[dno][kk], dcc1y[dno][kk]);
			    printf("kk %d daa1x: %E\t, daa1y: %E\n dbb1x %E, dbb1y %E\n Wall %d\n\n\n",kk, daa1x[dno][kk], daa1y[dno][kk], dbb1x[dno][kk], dbb1y[dno][kk], kk);
			    getchar();}
			}
		      
		      /*now I have the parameter calculated*/
		      par=licz/mian;
		      /*now I need to make decision rule*/
		      if(par>=-deltapar && par<=1.0 + deltapar)
			/*particle is close to the dust*/
			{
			  if((par>=-0.0) && (par<=1.0) && (t_cross<=dt)&& (t_cross>=0)) 
			    //Ready to assign particle :)
			    {
			      if((t_cross<timehitmin) || (timehitmin< -5*dt))
				{	dnohit=dno; 
				dnoseghit=kk;
				timehitmin=t_cross;
				par_tmin=par;			       
				//printf("hits %d, %E\n", l, t_cross);getchar();
			      }
			  }
			else
			  {
			    //Particle remain suspected
			    //add it to the list if not hit already
			    part_suspected=1;
			    if((t_cross<timehitmin_suspected) || (timehitmin_suspected < -5*dt))
			      {	
				dnohit_suspected=dno; 
				dnoseghit_suspected=kk;
				timehitmin_suspected=t_cross;
				par_tmin_suspected=par;
			      }	
			  }
			}
		    }
		}
	    }
	  	  

	  //LOOSE!
	  //do loose particle now
	  int licznik;

	  if(dnohit>=0)
	    { 
	      //	      printf("assigning normal\n");
	      if(dnoseghit<0)
		{
		  printf("SOMETHING WRONG in accel.c\n");
		  getchar();
		} 
	      //	      printf("part lost\n");
	      t_cross=timehitmin;
	      assigned=fopen("assign.dat", "a");
	      xhitp=partxold+partvxold*t_cross;
	      yhitp=partyold+partvyold*t_cross;
	      
	      fprintf(assigned, "%E\t%E\n", xhitp*normx, yhitp*normx); //getchar();
	      fclose(assigned);
	      
	      //Assign particle
	      double pvirtual;
	      int podloga, sufit;
	      double weight_low;
	      
       	      pvirtual=par_tmin*vipcorner[dnohit][dnoseghit];
	     
	      //number of the virtual particles on the given segment 
	      podloga=(int)(pvirtual); 	     
	      sufit=podloga+1;
	      weight_low=(pvirtual-podloga);
	      //  weight_high=(1-weight_low);
	      //find linear weighting
	      //assign now	    
	      // particles normalized
	      if(dnoseghit==0)
		{
		  licznik=0;
		}
	      else
		{
		  licznik=ccorner[dnohit][dnoseghit-1];
		}

	      dpart[dnohit][licznik+podloga].q+=normalcharge[i]*(1.0-weight_low)/ratio;
	      dpart[dnohit][licznik+sufit].q+=normalcharge[i]*weight_low/ratio;
	      		    		
	     
	      /*particle lost*/
	      current[dnohit][i]++; //CORRECT THE CURRENT ARRAY
	      //printf("getchere\n");  
	      lostpart[i]++;  
	      lostlist[i][lostpart[i]-1]=l;
	      //	printf("loastpart %d\n", lostpart[i]); getchar();
	      //    printf("particle is lost\n");
	    }
	  else
	    {
	      //check if particle is suspected at least
	      if(part_suspected)	      
		if(initpartcheck(spec[i].part[ii].x, spec[i].part[ii].y, spec[i].part[ii].z, 0.000001*dx) % 2 != 0)
		{
		  //  printf("Assigning suspected particle\n");
		  //Assign particle
		  double pvirtual;
		  int podloga, sufit;
		  double weight_low;
		  t_cross=timehitmin_suspected;
		  assigned=fopen("assign_susp.dat", "a");
		  xhitp=partxold+partvxold*t_cross;
		  yhitp=partyold+partvyold*t_cross;
		  fprintf(assigned, "%E\t%E\n", xhitp*normx, yhitp*normx); //getchar();
		  fclose(assigned);	 

		  pvirtual=par_tmin_suspected*vipcorner[dnohit_suspected][dnoseghit_suspected];
		  //number of the virtual particles on the given segment 
		  podloga=(int)(pvirtual); 	     
		  sufit=podloga+1;
		  weight_low=(pvirtual-podloga);
		  //  weight_high=(1-weight_low);
		  //find linear weighting
		  //assign now	    
		  // particles normalized
		  // printf("here\n");
		  if(par_tmin_suspected < 0)
		    podloga=sufit=0;
		  if(par_tmin_suspected > 1)
		    podloga=sufit=vipcorner[dnohit_suspected][dnoseghit_suspected];
		  
		  if(dnoseghit==0)
		    {
		      licznik=0;
		    }
		  else
		    {
		      licznik=ccorner[dnohit_suspected][dnoseghit_suspected-1];
		    }


		  dpart[dnohit_suspected][licznik+podloga].q+=normalcharge[i]*(1.0-weight_low)/ratio;
		  dpart[dnohit_suspected][licznik+sufit].q+=normalcharge[i]*weight_low/ratio;
		  
		  /*particle lost*/
		  current[dnohit_suspected][i]++; //CORRECT THE CURRENT ARRAY
		  lostpart[i]++;  
		  lostlist[i][lostpart[i]-1]=l;	
		  // printf("assigned\n");	  
		}	 	      	      
	    }
	  //	 	  printf("lost on dust %d\n", lostpart[i]);

	  //outer boundaries   
	  if((spec[i].part[ii].x >= Lx) || (spec[i].part[ii].x <= 0) || (spec[i].part[ii].y >= Ly) || (spec[i].part[ii].y <= 0) || (spec[i].part[ii].z >= Lz) || (spec[i].part[ii].z <= 0))
	    {      
	      lostpart[i]++;  	   
	      lostlist[i][lostpart[i]-1]=l; //create the list of lost particles
	    }
	  //	  printf("aft and dt %E\n", dt);
	  //  printf("part %d done\n", l);	
	}
      
      printf("finished now shifting SPECIE %d\n", i);
      // printf("end\n");
      /*shift particles*/
      int allpart;
      allpart=npart[i];
      if((numtasks>1 && rank==1)||(numtasks==1 && rank==0))
	{
	  printf("allpart %d, lostpart %ld\n", allpart, lostpart[i]);
	  // fprintf(history,"TIME %E, specie %d, allpart %d, lostpart %d\n",t*dt*normtime, i,  allpart, lostpart[i]);	  
	}
      j=0;
      while(j<lostpart[i])
	{
	  if(lostlist[i][lostpart[i]-1]!=allpart-1)
	    {
	      l=lostlist[i][j];
	      spec[i].part[ii]=spec[i].part[allpart-1]; //put the last to gap
	      allpart--; //decrease particles
	      j++;
	    }
	  else
	    {
	      allpart--; //Go to the last particle in the row
	      lostpart[i]--; //and forget last part which is lost
	    }
	}          
      npart[i]=allpart;   
    }

  //  printf("finished move(): current electr %d ions %d\n", current[0][0], current[0][1]);
}
