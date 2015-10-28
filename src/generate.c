/*DUSTY PLASMA program*/
/*particle generators*/
/*Author: Wojciech Jacek Miloch*/
/*University of Oslo, Norway*/
/*June 2006*/
#include <stdlib.h>
#include <math.h>
#include "const.h"
/*generate background particles in the chamber*/
void gen_bgnd(void)
{
  FILE *ffp;
  
  int i,j; /*counters*/
  double qx,qy,qz; /*coeff. for positions*/
  double lx,ly,lz,r2; /*coeff. for velocities*/
  int where;

  ffp=my_file_open("cos.txt", "w");
  for(i=0; i<S; i++) /*for each specie*/
    {  
	   init_primeroot(rank/(numtasks*1.0)); //ZAMIEN TO !!! DO GORY
      //  printf("primeroot: %3.60E", primeroot()); getchar();
     // init_primeroot(0.75);
      npart[i]=npartinit[i]; //already calcul in input.c
	//	if(i==0) npart[i]=1;
	//	if(i==1) npart[i]=0;
		
      for(j=0; j<npart[i]; j++) /*for each particle*/
	{
	  do //this loop is for not having particles in the object :)
	    {
	      where=1;
	      qx=primeroot(); /*locate the particle*/      
	      qy=primeroot();  
		  qz=primeroot(); 
	      spec[i].part[j].x=qx*Lx;
	      spec[i].part[j].y=qy*Ly;
		  spec[i].part[j].z=qz*Lz;
			
			
		 //  spec[i].part[j].x=dustxdx[0][1]-0.0000011*dx;
	     // spec[i].part[j].y=(dustydy[0][1]+dustydy[0][2])*0.5; 
		  
	      //for(z=0; z<(probesegments/4); z++)
		//	if((spec[i].part[j].x >= minx[2+4*z]) && (spec[i].part[j].x <= maxx[3+4*z] ) && (spec[i].part[j].y >=miny[0+4*z]) && (spec[i].part[j].y <= maxy[1+4*z]))
	      if(initpartcheck(spec[i].part[j].x, spec[i].part[j].y, spec[i].part[j].z, 0.000001*dx) % 2 != 0)	
		   {
		     where=0; 
		   }
	    }	  
	  while(where==0);  
	  /*Find the Maxwellian velocity*/	  
	  do
	    {
	      lx=primeroot()*2-1;
	      ly=primeroot()*2-1;	     
	      r2=lx*lx+ly*ly;
	    }
	  while(r2>=1);
	  
	  spec[i].part[j].vx = vdriftx[i]+vthx[i]*lx*sqrt(-2*log(r2)/r2);
	  spec[i].part[j].vy=vthy[i]*ly*sqrt(-2*log(r2)/r2);
	  
	   /*Find the Maxwellian velocity in z direction*/	  
	  do
	    {
	      lx=primeroot()*2-1;
	      lz=primeroot()*2-1;	     
	      r2=lx*lx+lz*lz;
	    }
	  while(r2>=1);
	  	  
	  spec[i].part[j].vz=vthz[i]*lz*sqrt(-2*log(r2)/r2);
		
		//TESTTT
	//	spec[i].part[j].x=0.5*Lx;
	//	spec[i].part[j].y=0.5*Ly;
	//	spec[i].part[j].z=0.65*Lz;
	//	spec[i].part[j].vx=0.50*vthx[1];
	//	spec[i].part[j].vy=0.0;
	//	spec[i].part[j].vz=0.0;

	   //if(spec[i].part[j].x <= dustxdx[0][1] && spec[i].part[j].x >= dustxdx[0][0] && spec[i].part[j].y >= dustydy[0][1] && spec[i].part[j].y <= dustydy[0][2])
	  //fprintf(ffp, "%E\t%E\t%E\n", spec[i].part[j].vx, spec[i].part[j].vy, spec[i].part[j].vz );	

	}   
    }
  fclose(ffp);
}

/*function to generate new particles*/
void newparticles(int timestep)
{
  int i,j;
  int last;
  double x,y,z,rsq;
  double iks1,r1,r2,r3;
  int side;
  double v, ub, vold, vtemp;
  int extra;
  //srand48(0);
  //init_primeroot(0);

    //	for(i=0; i<-1; i++) //empty TESTTT
  for(i=0; i<S; i++) /*for each specie do*/
	{ 
      allpart=npart[i];
      last=npart[i];
      /*check if enough space for new particles*/
      if(last+totalflux[0]>=NPART_MAX)
	{
	  printf("Problem, too many particles, rank %d\n. I am dumping them now\n", rank);
	  //  dump();
	  printf("I am terminating program and dumping particles.\nTry to increase N_PARTMAX and start again\n");
	  exit(1);
	}
      /*inject new particles*/
	  /*there are 6 sides now*/
	  
      side=0; //side 0 = left, positive x
      extrapart[i][side]+=fluxrest[i][side];
      extra=floor(extrapart[i][side]);	
      if(extra==1)
    	extrapart[i][side]-=extra; 
      //printf("in generate: last %d %E fluxside %E %E\n", last, primeroot(), flux[i][side], extrapart[i][side]);
      for(j=last; ((j-last)<(flux[i][side]+extra-1)); j++)
	{                                  
	  iks1=1.0-primeroot(); //substraction is to avoid 0.0 range (0,1]
	  r1=primeroot(); // [0,1]
	  r2=primeroot();
	  r3=primeroot();
		//printf("pos %E %E %E\n", r1, r2, r3);	
	  do
	    {	     //find maxwellian in y and z direction
	      z=primeroot()*2-1;
	      y=primeroot()*2-1;
	      rsq=z*z+y*y;
			//printf("vel %E %E\n", z, y);	
	    }
	  while(rsq>=1);	      
	 
	  iks1=iks1*lzet[i];
	  vold=lv0[i];
#ifdef BEAM
	  if(i!=2)
#endif
	    do
	      {
		ub=(vold-vdriftx[i])/(sqrt_two*vthx[i]);
		v=vold-((vthx[i]*vthx[i]*(exp(-llb[i]*llb[i])-exp(-ub*ub))+0.5*vdriftx[i]*vthx[i]*sqrt_twopi*(erfcc(llb[i])-erfcc(ub)))-iks1)/(vold*exp(-((vold-vdriftx[i])*(vold-vdriftx[i]))/(2*vthx[i]*vthx[i])));
		vtemp=v-vold;
		vold=v;  
	      }
	    while(fabs(vtemp)>0.00001*lv0[i]); 
#ifdef BEAM
	  else
	    v=vdriftx[i];
#endif	  
	  spec[i].part[j].vy=vthy[i]*y*sqrt(-2*log(rsq)/rsq); //maxwellian in y
	  spec[i].part[j].y=r2*Ly;
	  
	  spec[i].part[j].vx=v;	   //calculated in x
	  spec[i].part[j].x=r1*spec[i].part[j].vx*dt; //left

      spec[i].part[j].vz=vthz[i]*z*sqrt(-2*log(rsq)/rsq); //maxwellian in z
	  spec[i].part[j].z=r3*Lz;

	//if(spec[i].part[j].x >= dustxdx[0][0] && spec[i].part[j].x <= dustxdx[0][1])
      // 	if(spec[i].part[j].y >= dustydy[0][1] && spec[i].part[j].y <= dustydy[0][2])
//	{
//	printf("FLUX GENERATION Iam inside\n"); getchar();
//	}
// printf("side %d, %f %f %f\n%f %f %f\n",side, spec[i].part[j].x, spec[i].part[j].y, spec[i].part[j].z,spec[i].part[j].vx, spec[i].part[j].vy, spec[i].part[j].vz );
//	  getchar();
	  allpart++;
//	   if(j==124717) {printf("in generate side %d: %E %E %E",side, spec[i].part[j].x,spec[i].part[j].y,spec[i].part[j].z); getchar();}

	}
      last=allpart;

      fprintf(history, " %d\t",extra);

      
      side=1; //side 1 = right - flux in negative x direction
      extrapart[i][side]+=fluxrest[i][side];
      extra=floor(extrapart[i][side]);
      if(extra==1)
	extrapart[i][side]-=extra;      
      for(j=last; ((j-last)<(flux[i][side]+extra-1)); j++)
	{                                  
	  iks1=1-primeroot(); //substraction is to avoid 0.0 range (0,1]
	  r1=primeroot(); // [0,1]
	  r2=primeroot();
	  r3=primeroot();
	  do
	    {
	      z=primeroot()*2-1;
	      y=primeroot()*2-1;
	      rsq=z*z+y*y;
	    }
	  while(rsq>=1);	      
	  
       	  iks1=iks1*rzet[i];
	  //here drift must be made opposite - rvdrift
	  vold=rv0[i];
#ifdef BEAM
	  if(i!=2)
#endif
	    do
	      {
		ub=(vold-rvdriftx[i])/(sqrt_two*vthx[i]);
		v=vold-((vthx[i]*vthx[i]*(exp(-rlb[i]*rlb[i])-exp(-ub*ub))+0.5*rvdriftx[i]*vthx[i]*sqrt_twopi*(erfcc(rlb[i])-erfcc(ub)))-iks1)/(vold*exp(-((vold-rvdriftx[i])*(vold-rvdriftx[i]))/(2*vthx[i]*vthx[i])));
		vtemp=v-vold;
		vold=v;  
	      }
	    while(fabs(vtemp)>0.00001*rv0[i]); 
#ifdef BEAM
	  else
	    v=-vdriftx[i];
#endif
	  spec[i].part[j].vy=vthy[i]*y*sqrt(-2*log(rsq)/rsq); //maxwellian in y
	  spec[i].part[j].y=r2*Ly;		    
	  
	  spec[i].part[j].vz=vthz[i]*z*sqrt(-2*log(rsq)/rsq); //maxwellian in z
	  spec[i].part[j].z=r3*Lz;		    
	  
	  spec[i].part[j].vx=-v;
	  spec[i].part[j].x=Lx+r1*spec[i].part[j].vx*dt; //right
//	  printf("side %d, %f %f %f\n%f %f %f\n",side, spec[i].part[j].x, spec[i].part[j].y, spec[i].part[j].z,spec[i].part[j].vx, spec[i].part[j].vy, spec[i].part[j].vz );
//	  getchar();
// if(j==124717) {printf("in generate side %d: %E %E %E",side, spec[i].part[j].x,spec[i].part[j].y,spec[i].part[j].z); getchar();}

	  allpart++;
	}
      last=allpart;

   fprintf(history, " %d\t",extra);

      
      side=2; //side 2 = in positive y direction
      extrapart[i][side]+=fluxrest[i][side];
      extra=floor(extrapart[i][side]);
      if(extra==1)
	extrapart[i][side]-=extra;      
      for(j=last; ((j-last)<(flux[i][side]+extra-1)); j++)
	{
	  iks1=1.0-primeroot(); //substraction is to avoid 0.0 range (0,1]
	  r1=primeroot(); // [0,1]
	  r2=primeroot();	
	  r3=primeroot();  
	  do
	    {
	      x=primeroot()*2-1;
	      z=primeroot()*2-1;	      
	      rsq=x*x+z*z;
	    }
	  while(rsq>=1);
	  spec[i].part[j].vx=vdriftx[i]+vthx[i]*x*sqrt(-2*log(rsq)/rsq); //maxwellian in x
	  spec[i].part[j].x=r1*Lx;
	
	  spec[i].part[j].vz=vthx[i]*z*sqrt(-2*log(rsq)/rsq); //maxwellian in z
	  spec[i].part[j].z=r3*Lz;
	  
	  spec[i].part[j].vy=vthy[i]*sqrt(-2*log(iks1)); //positive if bottom, 
	  spec[i].part[j].y=r2*spec[i].part[j].vy*dt; //bottom	  
	  allpart++;
	  // if(j==124717) {printf("in generate side %d: %E %E %E",side, spec[i].part[j].x,spec[i].part[j].y,spec[i].part[j].z); getchar();}

	  
	}  
last=allpart;
	
    side=3; //side 3 = in negative y direction
      extrapart[i][side]+=fluxrest[i][side];
      extra=floor(extrapart[i][side]);
      if(extra==1)
	extrapart[i][side]-=extra;      
      for(j=last; ((j-last)<(flux[i][side]+extra-1)); j++)
	{
	  iks1=1.0-primeroot(); //substraction is to avoid 0.0 range (0,1]
	  r1=primeroot(); // [0,1]
	  r2=primeroot();	
	  r3=primeroot();  
	  do
	    {
	      x=primeroot()*2-1;
	      z=primeroot()*2-1;	      
	      rsq=x*x+z*z;
	    }
	  while(rsq>=1);
	  spec[i].part[j].vx=vdriftx[i]+vthx[i]*x*sqrt(-2*log(rsq)/rsq); //maxwellian in x
	  spec[i].part[j].x=r1*Lx;
	
	  spec[i].part[j].vz=vthz[i]*z*sqrt(-2*log(rsq)/rsq); //maxwellian in z
	  spec[i].part[j].z=r3*Lz;
	  
	  spec[i].part[j].vy=-vthy[i]*sqrt(-2*log(iks1)); //positive if bottom, 
	  spec[i].part[j].y=Ly+r2*spec[i].part[j].vy*dt; //top	  
	  allpart++;
	   //printf("%d allpart %d\n",j, allpart);

	//   if(j==124717) {printf("in generate side %d: %E %E %E",side, spec[i].part[j].x,spec[i].part[j].y,spec[i].part[j].z); getchar();}

	 

	}  
	last=allpart;		 					           
    fprintf(history, " %d\n",extra);

     side=4; //side 4 = through z in positive z
// printf("\n side %d last %d\n", side, last);
	       extrapart[i][side]+=fluxrest[i][side];
      extra=floor(extrapart[i][side]);
	if(extra==1)
	  extrapart[i][side]-=extra;
      for(j=last; ((j-last)<(flux[i][side]+extra-1)); j++)
	{
	  iks1=1.0-primeroot(); //substraction is to avoid 0.0 range (0,1]  
	  r1=primeroot(); // [0,1]
	  r2=primeroot(); 
	  r3=primeroot();  
	  do
	    {
	      x=primeroot()*2-1;
	      y=primeroot()*2-1;	 
	      rsq=x*x+y*y;
	    }
	  while(rsq>=1);
	  
	  spec[i].part[j].vx=vdriftx[i]+vthx[i]*x*sqrt(-2*log(rsq)/rsq); //maxwellian in x
	  spec[i].part[j].x=r1*Lx;
	  
      spec[i].part[j].vy=vthy[i]*y*sqrt(-2*log(rsq)/rsq); //maxwellian in y
	  spec[i].part[j].y=r2*Ly;
	  
	  spec[i].part[j].vz=vthz[i]*sqrt(-2*log(iks1)); //positive if bottom, 
	  spec[i].part[j].z=r3*spec[i].part[j].vz*dt; //bottom
	
	  allpart++;	
	   //printf("%d ",j);
	// 	 if(j==124717) {printf("in generate side %d: %E %E %E",side, spec[i].part[j].x,spec[i].part[j].y,spec[i].part[j].z); getchar();}


	}
      last=allpart;
      fprintf(history, "extra %d part: %d\t",i,extra);


   side=5; //side 5 =  through z in negative z direction
      extrapart[i][side]+=fluxrest[i][side];
      extra=floor(extrapart[i][side]);
      if(extra==1)
	  extrapart[i][side]-=extra;      
      for(j=last; ((j-last)<(flux[i][side]+extra-1)); j++)
	{
	  iks1=1.0-primeroot(); //substraction is to avoid 0.0 range (0,1]
	  r1=primeroot(); // [0,1]
	  r2=primeroot();	
	  r3=primeroot();  
	  do
	    {
	      x=primeroot()*2-1;
	      y=primeroot()*2-1;	      
	      rsq=x*x+y*y;
	    }
	  while(rsq>=1);
	  spec[i].part[j].vx=vdriftx[i]+vthx[i]*x*sqrt(-2*log(rsq)/rsq); //maxwellian in x
	  spec[i].part[j].x=r1*Lx;
	
	  spec[i].part[j].vy=vthy[i]*y*sqrt(-2*log(rsq)/rsq); //maxwellian in y
	  spec[i].part[j].y=r2*Ly;
	  
	  spec[i].part[j].vz=-vthz[i]*sqrt(-2*log(iks1)); //positive if bottom, 
	  spec[i].part[j].z=Lz+r3*spec[i].part[j].vz*dt; //top	  
	  allpart++; 
	   //printf("side %d, specie %d, %f %f %f\n%f %f %f\n",side, i, spec[i].part[j].x, spec[i].part[j].y, spec[i].part[j].z,spec[i].part[j].vx, spec[i].part[j].vy, spec[i].part[j].vz );
	  //getchar();
	  //printf("%d ",j);
	  // if(j==124717) {printf("in generate side %d: %E %E %E",side, spec[i].part[j].x,spec[i].part[j].y,spec[i].part[j].z); getchar();}

	}  
	last=allpart;	
	      printf("in generate: after last %d\n", last);
	 	
	  				           
    fprintf(history, " %d\n",extra);

      lostpart[i]=0; 
      npart[i]=allpart;      
    }
	
}

/****PRIME ROOT FRACTION NUMBER GENERATOR****/
double primeroot(void)
{
  double r;
  int draw=floor((drand48()*BUCKETSIZE));
  //draw the number, generate new one, and put the number in an epty space
  r=primerootbucket[draw];
  primerootno=primerootno+sqrt_two-floor((primerootno+sqrt_two));
  //NOTE we could also use floor(x) here -> x>=0
  primerootbucket[draw]=primerootno;
  return r;
}
/*****PRIME ROOT FRACTION NUMBER INITIALIZER*****/
void init_primeroot(double seed)
{
  int i;
  double x;
  srand(rank);
  srand48(rank);
  x=seed;
  for(i=0; i<BUCKETSIZE; i++)
    {      
      x=x+sqrt_two-(int)(x+sqrt_two); //prime root fraction
      primerootbucket[i]=x;
    }
	primerootno=0.1;
}



/*************************************/
int initpartcheck(double px, double py, double pz, double delta)
{
  int crossed=0;
  int maybein=0;
  double distance;
  double a,b,x1,y1;
  int dno,vert,kk,solution,kp1;
    
  //as for now no dust at the boundaries py!=0, px!=0
  vert=0;
  if(px!=0) //this is always true
    {
	  //PARAMETERS FOR "POLYGONAL" OR "POLYHEDRAL" DUST GRAINS
      a=py/px;
      b=0;
      /*We found linear function for particle*/
      /*now we check if it intersects with probe segment*/	  
      delta=0.000001*dx; //necessary due to computational errors
      
	  //PARAMETERS FOR "SPHERICAL" GRAINS, 
	  
      /*check all dusts*/
  for(dno=0; dno<noofdusts; dno++)
  {
    if(dshape[dno]==0) //spherical dust grains
    	{
		 distance=(px-dustcxdx[dno])*(px-dustcxdx[dno])+(py-dustcydx[dno])*(py-dustcydx[dno])+(pz-dustczdx[dno])*(pz-dustczdx[dno]);
		 distance=sqrt(distance); 
		 distance-=delta; //due to computational errors
		 //	 printf("distance %E dradiusdx[] %E dno %d noofdusts %d\n", distance, dradiusdx[dno], dno, noofdusts);
		 if(distance<=dradiusdx[dno])
		   return 1; 
		 else
           return 2;
        }
  if(dshape[dno]==1) //polyhedral dust grains
   {
   /*
	for(kk=0; kk<ncorners[dno]; kk++)
	  {	
	    //   printf("IN CHECK corner %d v %d ax %E b %E\n", kk, dustv[dno][kk], dusta[dno][kk], dustbdy[dno][kk]); getchar();   

	    //  printf(" delta %E IN CHECK corner %d v %d ax %E b %E and %E %E and a %E b %E\n",delta, kk, dustv[dno][kk], dusta[dno][kk], dustbdy[dno][kk], dustxdx[dno][kk], dustydy[dno][kk], a, b); getchar();

	    solution=1;
	    if(vert==1)
	      if((dustv[dno][kk]==1)&&(dustxdx[dno][kk]!=px))
		solution=0;
	      else
		if(dustv[dno][kk]==1)
		  {
		    y1=py; 
		    x1=px;
		  } 
		else
		  {   
		    x1=px;
		    y1=dusta[dno][kk]*px+dustbdy[dno][kk]; 
		  }
	    else
	      {
		if(dustv[dno][kk]==1)
		  {
		    y1=a*dustxdx[dno][kk]+b;
		    x1=dustxdx[dno][kk];
		  }
		else
		  if((a-dusta[dno][kk])!=0) //they will cross somewhere
		    {
		      x1=(dustbdy[dno][kk]-b)/(a-dusta[dno][kk]);
		      y1=a*x1+b;
		      //printf("cross\n");
		      //  getchar();	    
		    }
		  else 
		    if((dustbdy[dno][kk]-b)!=0)
		      solution=0;
	      }

	    //check if they crossed
	    //CHANGE!!!
	    if(solution==1)
		  {
		    //    printf("px py %E %E x1, y1, %E, %E\n", px,py,x1, y1); getchar();
		    // printf("solution for wall%d x1 %f y1 %f for a %f and b %f\n", kk, x1, y1, a, b);
		    if(kk<ncorners[dno]-1)
		      kp1=kk+1;
		    else 
		      kp1=0;
			  
		    if(((x1+delta>=dustxdx[dno][kk])&&(x1-delta<=dustxdx[dno][kp1])&&(x1+delta>=0)&&(x1-delta<=px))||((x1+delta>=dustxdx[dno][kp1])&&(x1-delta<=dustxdx[dno][kk])&&(x1+delta>=0)&&(x1-delta<=px)))
		      {
    			if(((y1+delta>=dustydy[dno][kk])&&(y1-delta<=dustydy[dno][kp1])&&(y1+delta>=0)&&(y1-delta<=py))||((y1+delta>=dustydy[dno][kp1])&&(y1-delta<=dustydy[dno][kk])&&(y1+delta>=0)&&(y1-delta<=py)))
			  {
	   
			  crossed++;//particle hit the probe;
			  // if((crossed % 2)==0) //we may go out too early!!
			  if((x1-delta <= px) && (x1+delta >= px) && (y1 - delta <= py) && (y1+delta >= py)) //we may be still inside!! 		
			    {			       maybein=1; }
			  }
		      }  
			  

		  }
	  }
    }	
  //SOLVED THE PROBE
 
  if((crossed % 2)==0)
      if(maybein) 
	crossed--;

  return crossed;
  */
  }
  }
}
}



int initpartcheck_restart(int dno, double px, double py, double pz, double delta)
{
  int crossed=0;
  int maybein=0;
  double distance;
  double a,b,x1,y1;
  int vert,kk,solution,kp1;
    
  //as for now no dust at the boundaries py!=0, px!=0
  vert=0;
  if(px!=0) //this is always true
    {
      //PARAMETERS FOR "POLYGONAL" OR "POLYHEDRAL" DUST GRAINS
      a=py/px;
      b=0;
      /*We found linear function for particle*/
      /*now we check if it intersects with probe segment*/	  
      delta=0.000001*dx; //necessary due to computational errors
      
	  //PARAMETERS FOR "SPHERICAL" GRAINS, 
      
      /*check all dusts*/
      //  for(dno=0; dno<noofdusts; dno++)
      
      if(dshape[dno]==0) //spherical dust grains
    	{
	  distance=(px-dustcxdx[dno])*(px-dustcxdx[dno])+(py-dustcydx[dno])*(py-dustcydx[dno])+(pz-dustczdx[dno])*(pz-dustczdx[dno]);
	  distance=sqrt(distance); 
	  distance-=delta; //due to computational errors
	  //	  printf("delta %E distance %E dradiusdx[] %E dno %d noofdusts %d\n", delta, distance, dradiusdx[dno], dno, noofdusts);
	  if(distance<=dradiusdx[dno])
	    return 1; 
	  else
	    return 2;
        }
      if(dshape[dno]==1) //polyhedral dust grains
	{
	  
	}
      
    }
}
/*************************************/

/*************************************/

