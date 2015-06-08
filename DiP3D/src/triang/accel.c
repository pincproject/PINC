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
  double vec_ax, vec_ay, vec_az, vec_bx, vec_by, vec_bz, weight0, weight1, weight2, totalweight, tempdou, mindist[3], localtrianglearea, dist_centre, dist_sphere2;
int ct, swapped, list[3], tempint;		
  double qmratio;
  double partxold, partyold, partzold;
  double partxhit, partyhit, partzhit;
  double partnewx, partnewy, partnewz;
    double qfactor=0.0;
  int vert,kk;
  double q;
  int offset;
  int kp1;
  double partvxold, partvyold,partvzold;
  double partvxnew, partvynew,partvznew;  
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
	
  //int ikl;
		  //	   for(ikl=0; ikl<FsMAX; ikl++)
		  //Fs[ikl]=0.0; 
	   //	printf("cleared %E", Fs[ix(0,64,0,0)]);
	   //	printf("cleared %E", Fs[ix(FsEy,64,0,0)]);
	   //	printf("cleared %E", Fs[ix(FsEz,64,36,58)]);	
	   //	printf("index %d\n", FsEz+6*ngy*ngz+36*ngz+58);	 
	//this was not the best idea   
  //create_linkedlist();	
	
	
	
  for(i=0; i<S; i++)
    {
      printf("ACCEL rank %d: in move, specie %d, npart %d \n", rank,i, npart[i]);
      offset=KE_off*i;
      lostpart[i]=0;
      qmratio=qm[i];     
      q=chargeandnorm[i];	
											
      for(ii=0; ii<npart[i]; ii++)
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
		
	
	//  if(j >= 64 || k>=64 || l>=64)
//{		  printf("part %d  and ijk: %d %d %d\n", ii, l,k,l);
//printf("particle no:  %d position %E vs %E", ii, spec[i].part[ii].z, Lz);	
//getchar();}
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
	  if((spec[i].part[ii].x >= Lx) || (spec[i].part[ii].x <= 0) || (spec[i].part[ii].y >= Ly) || (spec[i].part[ii].y <= 0) || (spec[i].part[ii].z >= Lz) || (spec[i].part[ii].z <= 0))
	    {      
	      printf("Time %d Somehow I am outside part ii %d: %E %E %E\n", t, ii, spec[i].part[ii].x, spec[i].part[ii].y, spec[i].part[ii].z); //create the list of lost particles
	      getchar();
	    }
	  //if(ax != 0 || ay !=0 || az!=0)
	  //	{	printf("accel nonzero %d %d %d \n ax %E ay %E az %E\n", j,k,l,ax, ay, az);	
	
	  //printf("%E\n", Fs[ix(FsEz,j,k,l)]);
	  //printf("%E\n", Fs[ix(FsEz,j,k,l+1)]);
	  //printf("%E\n", Fs[ix(FsEz,j,k+1,l+1)]);
	  //printf("%E\n", Fs[ix(FsEz,j+1,k+1,l+1)]);
	  //printf("%E\n", Fs[ix(FsEz,j+1,k,l+1)]);
	  //printf("%E\n", Fs[ix(FsEz,j+1,k+1,l)]);
	  //printf("%E\n", Fs[ix(FsEz,j+1,k,l)]);
	  //printf("%E\n", Fs[ix(FsEz,j,k+1,l)]);
	  //getchar();}
		  
	   //Check PP part
	  //if(ii==100)
	    //	    printf("rank %d specie %d part %d ijk %d %d %d, coord %E, %E, %E \n", rank, i, ii, j,k,l, spec[i].part[ii].x, spec[i].part[ii].y, spec[i].part[ii].z);
		//getchar();

	 int index, llist;	
	 int considerforhitting=0;
		double rho0000, rho0001, rho0011, rho0111, rho0100, rho0101, rho0110, rho0010;
		int dj,dk,dl, dii;
		double dpx,dpy,dpz;
		double dist_gp_p;
		double Cf,ddx,ddy,ddz;
		Cf=normPP*qmratio;
		
	//		printf("charge 2 mass ratio of the particle %E\n", qm[i]);
		
		//check if the gp is marked
		if(d_globallist[ix(0,j,k,l)]==1)
		{	
			
	 //  {printf("message %d %d %d\n", j,k,l);}
		//	printf("dii %d %E\n", 41, dpart[0][41].q);
			//getchar();
		// then do detailed check	
			
	  for(dno=0; dno<noofdusts; dno++)
	  {	
//      printf("checkingdust\n");
		  for(llist=0; llist<d_localmax[dno]; llist++)	
	   {
	//	   printf("checkingdust %d \n", i);

		   index=d_locallist[dno][llist];
		 //  printf("i %d  index %d %d %d %d check %d\n", i, index, j,k,l, ix(0,j,k,l));
		   // calculate indices back
		   dj=(int)(1.0*index/(ngy*ngz));
		   dk=(int)((index-dj*ngy*ngz)*1.0/ngy);
		   dl=index-dj*ngy*ngz-dk*ngy;
		 //  printf("dust index %d %d %d %d check %d\n", index, dj,dk,dl, ix(0,dj,dk,dl));  
		//if(((ix(0,j,k,l-1)<=index)&&(index<=ix(0,j,k,l+1)) || (ix(0,j,k-ngz,l-1)<=index)&&(index<=ix(0,j,k-ngz,l+1)) || (ix(0,j,k+ngz,l-1)<=index)&&(index<=ix(0,j,k+ngz,l+1)) || (ix(0,j-ngy*ngz,k,l)<=index)&&(index<=ix(0,j+ngy*ngz,k,l))))
	//checking if in the box
		   if((j-1<=dj) && (dj<=j+1) && (k-1<=dk) && (dk<=k+1) && (l-1<=dl) && (dl<=l+1))
		   {
		    considerforhitting=1; //plasma particle within 3x3 square from the dust
			
		//	getchar();
		//	printf("part %d befor ax:%E ay:%E az:%E magnitude: %E\n", ii, ax, ay,az, sqrt(ax*ax+ay*ay+az*az));
			double orax, oray, oraz;
			orax=ax;
			oray=ay;
			oraz=az;
			  // printf("dt %E part %d before ax:%E ay:%E az:%E magnitude: %E\n", dt, ii, ax, ay,az, sqrt(ax*ax+ay*ay+az*az));
   
		//corrections PP (MD part), change counters!!!
			//for all dpart at dno that belong 2 index
			for(dii=0; dii<dpartlast[dno]; dii++) 
			{
			 dj=dpart[dno][dii].x/dx; //take the integer
			 dk=dpart[dno][dii].y/dy;
			 dl=dpart[dno][dii].z/dz; 
			if(ix(0,j,k,l)==index) //only for those with indices
			 {
				 //q is a number of elementary charges, ratio calculation is hidden in Cf
				 
				 //MPI part 
				 
				// q=dpart[dno][dii].q*dt/dV;
			   q=dpartq[dno][dii]*dt/dV;
				 //TEST
				// q=1*dt/dV; 
				// q=1/dV;
				// q=1.0;
				 dpx=dpart[dno][dii].x-dj*dx;
				 dpy=dpart[dno][dii].y-dk*dy;
				 dpz=dpart[dno][dii].z-dl*dz;
				 
				 //dpart[dno][dii].x=dj*dx+0.5*dx;
				 //dpart[dno][dii].y=dk*dy+0.5*dx;
				 //dpart[dno][dii].z=dl*dz+0.5*dx;
				// dpx=0.5*dx;
				 //dpy=0.5*dy;
				 //dpz=0.5*dz;
				 
				 //spec[i].part[ii].x=
				 //spec[i].part[ii].y=
				 //spec[i].part[ii].z=

				 //rho tells the part of the charge being assigned to the gp
				 rho0000=q*(dx-dpx)*(dy-dpy)*(dz-dpz);
				 rho0100=q*dpx*(dy-dpy)*(dz-dpz);
				 rho0110=q*dpx*dpy*(dz-dpz);
				 rho0010=q*(dx-dpx)*dpy*(dz-dpz);
				 rho0001=q*(dx-dpx)*(dy-dpy)*dpz;
				 rho0101=q*dpx*(dy-dpy)*dpz;
				 rho0111=q*dpx*dpy*dpz;
				 rho0011=q*(dx-dpx)*dpy*dpz; 
			//	 printf("total weight %E and q %E\n", (rho0000+rho0001+rho0010+rho0100+rho0110+rho0111+rho0101+rho0011)/q/dV, q);
				// if(q!=0) {printf("dii %d %E\n", dii, dpart[dno][dii].q); getchar();}
				 //subtract: plasma particle - 8 corners
				 
				//now consider distance for the plasma particle 
				 
				 //source is the grid point
			
				 
				 //corner 0000
				 ddx=(spec[i].part[ii].x-dj*dx);
				 ddy=(spec[i].part[ii].y-dk*dy);
				 ddz=(spec[i].part[ii].z-dl*dz);
				 dist_gp_p=ddx*ddx+ddy*ddy+ddz*ddz;
				 dist_gp_p=sqrt(dist_gp_p);
				 ax-=Cf*rho0000*ddx/(dist_gp_p*dist_gp_p*dist_gp_p);
				 ay-=Cf*rho0000*ddy/(dist_gp_p*dist_gp_p*dist_gp_p);
				 az-=Cf*rho0000*ddz/(dist_gp_p*dist_gp_p*dist_gp_p);
				
				// printf("part %d minus ax:%E ay:%E az:%E magnitude: %E\n", ii, ax-orax, ay-oray,az-oraz, sqrt(ax*ax+ay*ay+az*az));

				 
				 //corner 001
				 ddx=(spec[i].part[ii].x-dj*dx);
				 ddy=(spec[i].part[ii].y-dk*dy);
				 ddz=(spec[i].part[ii].z-(dl+1)*dz);
				 dist_gp_p=ddx*ddx+ddy*ddy+ddz*ddz;
				 dist_gp_p=sqrt(dist_gp_p);
				 ax-=Cf*rho0001*ddx/(dist_gp_p*dist_gp_p*dist_gp_p);
				 ay-=Cf*rho0001*ddy/(dist_gp_p*dist_gp_p*dist_gp_p);
				 az-=Cf*rho0001*ddz/(dist_gp_p*dist_gp_p*dist_gp_p);
				 
				 //corner 010
				 ddx=(spec[i].part[ii].x-dj*dx);
				 ddy=(spec[i].part[ii].y-(dk+1)*dy);
				 ddz=(spec[i].part[ii].z-dl*dz);
				 dist_gp_p=ddx*ddx+ddy*ddy+ddz*ddz;
				 dist_gp_p=sqrt(dist_gp_p);
				 ax-=Cf*rho0010*ddx/(dist_gp_p*dist_gp_p*dist_gp_p);
				 ay-=Cf*rho0010*ddy/(dist_gp_p*dist_gp_p*dist_gp_p);
				 az-=Cf*rho0010*ddz/(dist_gp_p*dist_gp_p*dist_gp_p);
				 
				 //corner 011
				 ddx=(spec[i].part[ii].x-dj*dx);
				 ddy=(spec[i].part[ii].y-(dk+1)*dy);
				 ddz=(spec[i].part[ii].z-(dl+1)*dz);
				 dist_gp_p=ddx*ddx+ddy*ddy+ddz*ddz;
				 dist_gp_p=sqrt(dist_gp_p);
				 ax-=Cf*rho0011*ddx/(dist_gp_p*dist_gp_p*dist_gp_p);
				 ay-=Cf*rho0011*ddy/(dist_gp_p*dist_gp_p*dist_gp_p);
				 az-=Cf*rho0011*ddz/(dist_gp_p*dist_gp_p*dist_gp_p);
				 
				 
				 //corner 100
				 ddx=(spec[i].part[ii].x-(dj+1)*dx);
				 ddy=(spec[i].part[ii].y-dk*dy);
				 ddz=(spec[i].part[ii].z-dl*dz);
				 dist_gp_p=ddx*ddx+ddy*ddy+ddz*ddz;
				 dist_gp_p=sqrt(dist_gp_p);
				 ax-=Cf*rho0100*ddx/(dist_gp_p*dist_gp_p*dist_gp_p);
				 ay-=Cf*rho0100*ddy/(dist_gp_p*dist_gp_p*dist_gp_p);
				 az-=Cf*rho0100*ddz/(dist_gp_p*dist_gp_p*dist_gp_p);
				 
				 //corner 101
				 ddx=(spec[i].part[ii].x-(dj+1)*dx);
				 ddy=(spec[i].part[ii].y-dk*dy);
				 ddz=(spec[i].part[ii].z-(dl+1)*dz);
				 dist_gp_p=ddx*ddx+ddy*ddy+ddz*ddz;
				 dist_gp_p=sqrt(dist_gp_p);
				 ax-=Cf*rho0101*ddx/(dist_gp_p*dist_gp_p*dist_gp_p);
				 ay-=Cf*rho0101*ddy/(dist_gp_p*dist_gp_p*dist_gp_p);
				 az-=Cf*rho0101*ddz/(dist_gp_p*dist_gp_p*dist_gp_p);
				 
				 //corner 110
				 ddx=(spec[i].part[ii].x-(dj+1)*dx);
				 ddy=(spec[i].part[ii].y-(dk+1)*dy);
				 ddz=(spec[i].part[ii].z-dl*dz);
				 dist_gp_p=ddx*ddx+ddy*ddy+ddz*ddz;
				 dist_gp_p=sqrt(dist_gp_p);
				 ax-=Cf*rho0110*ddx/(dist_gp_p*dist_gp_p*dist_gp_p);
				 ay-=Cf*rho0110*ddy/(dist_gp_p*dist_gp_p*dist_gp_p);
				 az-=Cf*rho0110*ddz/(dist_gp_p*dist_gp_p*dist_gp_p);
				 
				 //corner 111
				 ddx=(spec[i].part[ii].x-(dj+1)*dx);
				 ddy=(spec[i].part[ii].y-(dk+1)*dy);
				 ddz=(spec[i].part[ii].z-(dl+1)*dz);
				 dist_gp_p=ddx*ddx+ddy*ddy+ddz*ddz;
				 dist_gp_p=sqrt(dist_gp_p);
				 ax-=Cf*rho0111*ddx/(dist_gp_p*dist_gp_p*dist_gp_p);
				 ay-=Cf*rho0111*ddy/(dist_gp_p*dist_gp_p*dist_gp_p);
				 az-=Cf*rho0111*ddz/(dist_gp_p*dist_gp_p*dist_gp_p);
				
			//	 printf("part %d minus ax:%E ay:%E az:%E magnitude: %E\n", ii, ax-orax, ay-oray,az-oraz, sqrt(ax*ax+ay*ay+az*az));
				 orax=ax;
				 oray=ay;
				 oraz=az;
				 
				 //add : plasma particle - dust particle				 
				 ddx=(spec[i].part[ii].x-dpart[dno][dii].x);
				 ddy=(spec[i].part[ii].y-dpart[dno][dii].y);
				 ddz=(spec[i].part[ii].z-dpart[dno][dii].z);
				 dist_gp_p=ddx*ddx+ddy*ddy+ddz*ddz;
				 dist_gp_p=sqrt(dist_gp_p);
				 ax+=Cf*dpartq[dno][dii]*ddx/(dist_gp_p*dist_gp_p*dist_gp_p);
				 ay+=Cf*dpartq[dno][dii]*ddy/(dist_gp_p*dist_gp_p*dist_gp_p);
				 az+=Cf*dpartq[dno][dii]*ddz/(dist_gp_p*dist_gp_p*dist_gp_p);

			//	 printf("part %d added ax:%E ay:%E az:%E magnitude: %E\n", ii, ax-orax, ay-oray,az-oraz, sqrt(ax*ax+ay*ay+az*az));
				 

			 }
				
				
				
				
			//	printf("part %d after ax:%E ay:%E az:%E magnitude: %E\n", ii, ax, ay,az, sqrt(ax*ax+ay*ay+az*az));

			}
			}
			
			
			
	   }
	  }
			
			
	} 
	
		
	//	ax=0;
	//	ay=0;
	//	az=0;
		
//	for(dno=0; dno<noofdusts; dno++)
//		for(dii=0; dii<dpartlast[dno]; dii++) 
//		{
			
			
//				ddx=(spec[i].part[ii].x-dpart[dno][dii].x);
//				ddy=(spec[i].part[ii].y-dpart[dno][dii].y);
//				ddz=(spec[i].part[ii].z-dpart[dno][dii].z);
//				dist_gp_p=ddx*ddx+ddy*ddy+ddz*ddz;
//				dist_gp_p=sqrt(dist_gp_p);
			
	//			ax+=Cf*dpart[dno][dii].q*ddx/(dist_gp_p*dist_gp_p*dist_gp_p);
	//			ay+=Cf*dpart[dno][dii].q*ddy/(dist_gp_p*dist_gp_p*dist_gp_p);
	//			az+=Cf*dpart[dno][dii].q*ddz/(dist_gp_p*dist_gp_p*dist_gp_p);
//			printf("Cf %E \n %E %E %E\n", Cf, ax, ay, az);
//		}

		
	 /*accelerate!*/
	  partvxold=spec[i].part[ii].vx;
	  partvyold=spec[i].part[ii].vy;
	  partvzold=spec[i].part[ii].vz;
	  spec[i].part[ii].vx+=ax;
	  spec[i].part[ii].vy+=ay;
	  spec[i].part[ii].vz+=az;
		partvxnew=spec[i].part[ii].vx;
		partvynew=spec[i].part[ii].vy;
		partvznew=spec[i].part[ii].vz;
		
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
	 // printf("move now\n");
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
		
	  //testowy=fopen("testowy.txt", "a");
	  //fprintf(testowy,"%E\t%E\t%E\t%E\t%E\t%E\t%E\n", t*dt*normtime, ax*normx/(normtime*normtime), ay*normx/(normtime*normtime), az*normx/(normtime*normtime), spec[i].part[ii].x, spec[i].part[ii].y*normx, spec[i].part[ii].z*normx);
	  //fclose(testowy);
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
	  if(considerforhitting==1)
	  {  
	  for(dno=0; dno<noofdusts; dno++)
	    {
	      if(dshape[dno]==0) //spherical dust
		{
	
	
		  	  //
		  //check if particle suspected
		  //  inside=0;
		  // if(initpartcheck(spec[i].part[ii].x, spec[i].part[ii].y, spec[i].part[ii].z, 0.000001*dx)%2 !=0) inside=1;
		  //check if particle close to the sphere
		// dustcxdx[dno]=dustcydx[dno]=dustczdx[dno]=0;
		// dradiusdx[dno]=sqrt(8);
//partnewx=2;	
//partnewy=2;
//partnewz=0;
  
		    dist_centre=(partnewx-dustcxdx[dno])*(partnewx-dustcxdx[dno])+(partnewy-dustcydx[dno])*(partnewy-dustcydx[dno])+(partnewz-dustczdx[dno])*(partnewz-dustczdx[dno]);
  //   printf("distcentre %E\n", dist_centre);	
		    if(dist_centre <= 4*dradiusdx[dno]*dradiusdx[dno]) //we are quite close and can do other checks, particle is suspected: (we have distance squared)
		    {
//		       if(i==1)
//			   printf("ACCEL: in move2\n");
		      //calculate if the particle is within the sphere: ||Po+vdt  - Co|| <= dustradius -> find time.
		      //dust is static
			//  partvxold=-100;
			//  partvyold=-100;
			//  partvzold=0;
			//  partxold=2.01;
			//  partyold=2.01;
			//  partzold=0.0;
			  
		      quad_a=partvxnew*partvxnew+partvynew*partvynew+partvznew*partvznew;
		      quad_b=2*((partxold-dustcxdx[dno])*partvxnew+(partyold-dustcydx[dno])*partvynew+(partzold-dustczdx[dno])*partvznew);
	 	      quad_c=(partxold-dustcxdx[dno])*(partxold-dustcxdx[dno])+(partyold-dustcydx[dno])*(partyold-dustcydx[dno])+(partzold-dustczdx[dno])*(partzold-dustczdx[dno])-dradiusdx[dno]*dradiusdx[dno];
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
		  //    printf("t_cross %E\n", t_cross/dt);
			  
		      if((t_cross <= dt+deltat) && (t_cross >=-deltat)) 	//we hit the surface & need to find nearest points on the surface
			{
			
			  partxhit=partxold+partvxnew*t_cross;
			  partyhit=partyold+partvynew*t_cross;
			  partzhit=partzold+partvznew*t_cross;
			  //printf("rank %d part %d hitting %E %E %E\n", rank, ii, partxhit, partyhit, partzhit); 
			  if(i==0)
			  {
			 //   printf("we hit the surface at t_cross/dt %E\n", t_cross/dt); 
			 // printf("particle hit %E %E %E", partxhit, partyhit, partzhit);
			 //getchar();
			}  
		    //  dist_sphere2=(partxhit-dustcxdx[dno])*(partxhit-dustcxdx[dno])+(partyhit-dustcydx[dno])*(partyhit-dustcydx[dno])+(partzhit-dustczdx[dno])*(partzhit-dustczdx[dno]);
			 // printf("\n relative distance from the centre %E\n", sqrt(dist_sphere2)/dradiusdx[dno]);
			
			  dist_sphere2=(partxhit-dpart[dno][0].x)*(partxhit-dpart[dno][0].x)+(partyhit-dpart[dno][0].y)*(partyhit-dpart[dno][0].y)+(partzhit-dpart[dno][0].z)*(partzhit-dpart[dno][0].z);
			  mindist[0]=dist_sphere2;
			 // printf("%E\n", dist_sphere2);
			  dist_sphere2=(partxhit-dpart[dno][1].x)*(partxhit-dpart[dno][1].x)+(partyhit-dpart[dno][1].y)*(partyhit-dpart[dno][1].y)+(partzhit-dpart[dno][1].z)*(partzhit-dpart[dno][1].z);
			  mindist[1]=dist_sphere2;
			 // printf("%E\n", dist_sphere2);
			  dist_sphere2=(partxhit-dpart[dno][2].x)*(partxhit-dpart[dno][2].x)+(partyhit-dpart[dno][2].y)*(partyhit-dpart[dno][2].y)+(partzhit-dpart[dno][2].z)*(partzhit-dpart[dno][2].z);
			  mindist[2]=dist_sphere2;
			 // printf("%E\n", dist_sphere2);
			  list[0]=0;
			  list[1]=1;
			  list[2]=2;
			  //sort the list...
			  //   printf("SORTING mindist %E %E %E\n", mindist[0], mindist[1], mindist[2]);

			  do{
			
			    swapped=0;
			    for(ct=0; ct<2; ct++)
			      if(mindist[ct]>mindist[ct+1])
				{
				  tempdou=mindist[ct+1];
				  mindist[ct+1]=mindist[ct];
				  mindist[ct]=tempdou;
				  tempint=list[ct+1];
				  list[ct+1]=list[ct];
				  list[ct]=tempint;
				  swapped=1;
				}
			    }while(swapped==1);		
			//	printf("postSORTING mindist %E %E %E\n", mindist[0], mindist[1], mindist[2]);	
			    //search through all the points
			    for(ct=3; ct<dpartlast[dno]; ct++)
			      {
				dist_sphere2=(partxhit-dpart[dno][ct].x)*(partxhit-dpart[dno][ct].x)+(partyhit-dpart[dno][ct].y)*(partyhit-dpart[dno][ct].y)+(partzhit-dpart[dno][ct].z)*(partzhit-dpart[dno][ct].z);
              //   printf("%E\n", dist_sphere2);
				//     printf("LOOP mindist %E %E %E\n", mindist[0], mindist[1], mindist[2]);

				if(dist_sphere2<mindist[2]){
				  if(dist_sphere2<mindist[1]){
				    if(dist_sphere2<mindist[0]){
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
				      }
				  }
				  else
				    {
				      mindist[2]=dist_sphere2;
				      list[2]=ct;
				    }
				}
			      }
				  //printf("mindist %E %E %E\n", mindist[0], mindist[1], mindist[2]);
			    //printf("Now assigning, list[0] %d, list[1] %d, list[2] %d\n hit: %E %E %E\n", list[0], list[1], list[2], partxhit, partyhit, partzhit); 
			    //printf("P1: %E %E %E\n", dpart[dno][list[0]].x,dpart[dno][list[0]].y,dpart[dno][list[0]].z);
			    //printf("P2: %E %E %E\n", dpart[dno][list[1]].x,dpart[dno][list[1]].y,dpart[dno][list[1]].z);
			    //printf("P3: %E %E %E\n", dpart[dno][list[2]].x,dpart[dno][list[2]].y,dpart[dno][list[2]].z);

/* dpart[dno][list[0]].x=dpart[dno][list[0]].y=dpart[dno][list[0]].z=0.0;
  dpart[dno][list[1]].x=dpart[dno][list[1]].y=0.0;
  dpart[dno][list[1]].z=2.0;
   dpart[dno][list[2]].x=dpart[dno][list[2]].z=0.0;
   dpart[dno][list[2]].y=1.0;
   partxhit=0.50;
   partyhit=0.5;
partzhit=1.0;*/

			    //now assign the charge to the nearest points from the list
			    //Area of the triangle *2 from the cross product of two vectors: area=module of the cross product
			  //  vec_ax=dpart[dno][list[1]].x-dpart[dno][list[0]].x;
			  //  vec_ay=dpart[dno][list[1]].y-dpart[dno][list[0]].y;
			  //  vec_az=dpart[dno][list[1]].z-dpart[dno][list[0]].z;
			  //  vec_bx=dpart[dno][list[2]].x-dpart[dno][list[0]].x;
			  //  vec_by=dpart[dno][list[2]].y-dpart[dno][list[0]].y;
			  //  vec_bz=dpart[dno][list[2]].z-dpart[dno][list[0]].z;

			  //  localtrianglearea=sqrt((vec_ay*vec_bz-vec_az*vec_by)*(vec_ay*vec_bz-vec_az*vec_by)+(vec_ax*vec_bz-vec_az*vec_bx)*(vec_ax*vec_bz-vec_az*vec_bx)+(vec_ax*vec_by-vec_ay*vec_bx)*(vec_ax*vec_by-vec_ay*vec_bx));
				//printf("localtriangle %E\n", localtrianglearea);
			    //calculate partial weigths
			    vec_ax=dpart[dno][list[1]].x-partxhit;
			    vec_ay=dpart[dno][list[1]].y-partyhit;
			    vec_az=dpart[dno][list[1]].z-partzhit;
			    vec_bx=dpart[dno][list[2]].x-partxhit;
			    vec_by=dpart[dno][list[2]].y-partyhit;
			    vec_bz=dpart[dno][list[2]].z-partzhit;
			    weight0=sqrt((vec_ay*vec_bz-vec_az*vec_by)*(vec_ay*vec_bz-vec_az*vec_by)+(vec_ax*vec_bz-vec_az*vec_bx)*(vec_ax*vec_bz-vec_az*vec_bx)+(vec_ax*vec_by-vec_ay*vec_bx)*(vec_ax*vec_by-vec_ay*vec_bx));

			    //=sqrt((vec_ay*vec_bz-vec_az*vec_by)*(vec_ay*vec_bz-vec_az*vec_by)+(vec_ax*vec_bz-vec_az*vec_bx)*(vec_ax*vec_bz-vec_az*vec_bx)+(vec_ax*vec_by-vec_ay*vec_bx)*(vec_ax*vec_by-vec_ay*vec_bx));
			  
			    vec_ax=dpart[dno][list[0]].x-partxhit;
			    vec_ay=dpart[dno][list[0]].y-partyhit;
			    vec_az=dpart[dno][list[0]].z-partzhit;
			    vec_bx=dpart[dno][list[2]].x-partxhit;
			    vec_by=dpart[dno][list[2]].y-partyhit;
			    vec_bz=dpart[dno][list[2]].z-partzhit;
			    weight1=sqrt((vec_ay*vec_bz-vec_az*vec_by)*(vec_ay*vec_bz-vec_az*vec_by)+(vec_ax*vec_bz-vec_az*vec_bx)*(vec_ax*vec_bz-vec_az*vec_bx)+(vec_ax*vec_by-vec_ay*vec_bx)*(vec_ax*vec_by-vec_ay*vec_bx));

			    // sqrt((vec_ay*vec_bz-vec_az*vec_by)*(vec_ay*vec_bz-vec_az*vec_by)+(vec_ax*vec_bz-vec_az*vec_bx)*(vec_ax*vec_bz-vec_az*vec_bx)+(vec_ax*vec_by-vec_ay*vec_bx)*(vec_ax*vec_by-vec_ay*vec_bx));

			    vec_ax=dpart[dno][list[1]].x-partxhit;
			    vec_ay=dpart[dno][list[1]].y-partyhit;
			    vec_az=dpart[dno][list[1]].z-partzhit;
			    vec_bx=dpart[dno][list[0]].x-partxhit;
			    vec_by=dpart[dno][list[0]].y-partyhit;
			    vec_bz=dpart[dno][list[0]].z-partzhit;
			    weight2=sqrt((vec_ay*vec_bz-vec_az*vec_by)*(vec_ay*vec_bz-vec_az*vec_by)+(vec_ax*vec_bz-vec_az*vec_bx)*(vec_ax*vec_bz-vec_az*vec_bx)+(vec_ax*vec_by-vec_ay*vec_bx)*(vec_ax*vec_by-vec_ay*vec_bx));
			    //sqrt((vec_ay*vec_bz-vec_az*vec_by)*(vec_ay*vec_bz-vec_az*vec_by)+(vec_ax*vec_bz-vec_az*vec_bx)*(vec_ax*vec_bz-vec_az*vec_bx)+(vec_ax*vec_by-vec_ay*vec_bx)*(vec_ax*vec_by-vec_ay*vec_bx));

totalweight=weight0+weight1+weight2;

			//    printf("partical weights calculated %E %E %E, test %E\n", (weight0/totalweight), (weight1/totalweight), (weight2/totalweight), (weight0+weight1+weight2)/totalweight);
			    //assign charge to nearest grid points
			 
			    qfactor=normalcharge[i]/ratio;			 
				//dpart[dno][list[0]].q+=qfactor*(weight0/totalweight);
				//*normalcharge[i];
			    //dpart[dno][list[1]].q+=qfactor*(weight1/totalweight);
				//*normalcharge[i];
			    //dpart[dno][list[2]].q+=qfactor*(weight2/totalweight);
				//*normalcharge[i];

				dpart[dno][list[0]].q+=qfactor*(weight0/totalweight);
				//*normalcharge[i];
			    dpart[dno][list[1]].q+=qfactor*(weight1/totalweight);
				//*normalcharge[i];
			    dpart[dno][list[2]].q+=qfactor*(weight2/totalweight);
				//*normalcharge[i];
				
				//CALCULATE DRAG FORCE
				//particle velocities, specie, 
				drag_force_direct(partvxnew, partvynew, partvznew, i, dno, partxhit, partyhit, partzhit);
			//	printf("in drag force specie %d part %d", i, ii);
			//	getchar();
				
//printf("assigning\n");
			    /*particle lost*/
			    current[dno][i]++; //CORRECT THE CURRENT ARRAY
			    lostpart[i]++;  
			    lostlist[i][lostpart[i]-1]=ii;	
//				printf("Now finished assigning \n"); 
			
			}			  
		    }
		}
	    }
		
	  }	
			  //outer boundaries   
	  if((spec[i].part[ii].x >= Lx) || (spec[i].part[ii].x <= 0) || (spec[i].part[ii].y >= Ly) || (spec[i].part[ii].y <= 0) || (spec[i].part[ii].z >= Lz) || (spec[i].part[ii].z <= 0))
	    {      
	      lostpart[i]++;  	   
	      lostlist[i][lostpart[i]-1]=ii; //create the list of lost particles
		  
		 if(spec[i].part[ii].x >= Lx) c0[i]++;
		 if(spec[i].part[ii].x <= 0) c1[i]++;
		 if(spec[i].part[ii].y >= Ly) c2[i]++;
		 if(spec[i].part[ii].y <= 0) c3[i]++;
		 if(spec[i].part[ii].z >= Lz) c4[i]++;
		 if(spec[i].part[ii].z <= 0) c5[i]++;
	    }
	  //	  printf("aft and dt %E\n", dt);
	  //  printf("part %d done\n", l);	
	}
      
      printf("ACCEL:  %d: Finished now shifting SPECIE %d\n", rank, i);
     
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
	      ii=lostlist[i][j];
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
	

	
	
	if(noofdusts>0)
   printf("finished move(): current electr %d ions %d average %d %d\n", current[0][0], current[0][1], curr_av[0][0], curr_av[0][1]);
else
	printf("finished move(): no dust\n");
	//for(i=0; i<dpartlast[0]; i++)
	 // printf("test po move rank %d i %d part %E\n", rank, i, dpart[0][i].q);
	//getchar();
    }
	
void create_linkedlist(void)
{
 int i,ii,j,k,l;
 long int index;
 for(ii=0; ii<S*llsize; ii++)
	 {
	   llmesh[ii]=-1;
	 }
	 
 for(i=0;i<S;i++) //for each specie
   {
	for(ii=0;ii<npart[i];ii++) 
	{
	  //locate particle on the coarse grid for linked list
	  j=spec[i].part[ii].x/lldx;
	  k=spec[i].part[ii].y/lldy;
	  l=spec[i].part[ii].z/lldz;
	  index=i*llsize*j*llngy*llngz+k*llngz+l; //offset + index

      spec[i].part[ii].llnext=llmesh[index];
	  llmesh[index]=ii;
	}
  }
}
