 /* DUSTY PLASMA */
/* Various calculation for dust grains*/
/* Author: Wojciech Jacek Miloch */
/* University of Oslo */
/* 2007 */
#include "const.h"

/***dust movement***/
void d_move(int t)
{
  int i,ii,j,k,l;
  double transx,transy, transz;
  double tempx, tempy, tempz;
  double sinbeta, cosbeta, r;
  double x,y,z,x1,y1,z1;
  double wjkl,wj1kl,wj1k1l,wjk1l,wjkl1,wj1kl1,wj1k1l1,wjk1l1; //weights
  double ftransX, ftransY, ftransZ, torque;

#ifdef MPI
 double rftransX, rftransY, rtorque;
#endif
  double fpar;
  double forceX, forceY, forceZ;
  double cosphi,sinphi;
  int jp1;

  double normforce=normmass*normx/(normtime*normtime);
  double normtorque=normforce*normx;
  //  printf("in dust move\n noofdusts %d \n", noofdusts); getchar();

  for(i=0; i<noofdusts; i++)
    {
      ftransX=ftransY=ftransZ=torque=0.0;
	
  
      /*FIND FORCEs*/
      /*take all particles assigned to dust[i]-> weight force as usual*/
     

	   for(ii=0; ii<dpartlast[i]; ii++) 
	{  
	  //no particles on rank 0 if MPI 
//	  printf("dpart: x %E y %E z %E\n", dpart[i][ii].x, dpart[i][ii].y, dpart[i][ii].z);
	 
	  x=dpart[i][ii].x;
	  y=dpart[i][ii].y;
	  z=dpart[i][ii].z;
	  
	  //find forces from the closest grid points,and locate particle:)
	  j=dpart[i][ii].x/dx;
	  k=dpart[i][ii].y/dy;
	  l=dpart[i][ii].z/dz;
	  x1=dpart[i][ii].x-j*dx;
	  y1=dpart[i][ii].y-k*dy;
	  z1=dpart[i][ii].z-l*dz;
	  /*find weights, E is multiplied already by dt/dxdy*/
	  wjkl=(dx-x1)*(dy-y1)*(dz-z1);
	  wj1kl=x1*(dy-y1)*(dz-z1);
	  wj1k1l=x1*y1*(dz-z1);
	  wjk1l=(dx-x1)*y1*(dz-z1);
	  wjkl1=(dx-x1)*(dy-y1)*z1;
	  wj1kl1=x1*(dy-y1)*z1;
	  wj1k1l1=x1*y1*z1;
	  wjk1l1=(dx-x1)*y1*z1;	  	  
	  forceX=dpart[i][ii].q*(Fs[ix(0,j,k,l)]*wjkl+Fs[ix(0,j+1,k,l)]*wj1kl+Fs[ix(0,j+1,k+1,l)]*wj1k1l+Fs[ix(0,j,k+1,l)]*wjk1l+Fs[ix(0,j,k,l+1)]*wjkl1+Fs[ix(0,j+1,k,l+1)]*wj1kl1+Fs[ix(0,j+1,k+1,l+1)]*wj1k1l1+Fs[ix(0,j,k+1,l+1)]*wjk1l1);
	  forceY=dpart[i][ii].q*(Fs[ix(FsEy,j,k,l)]*wjkl+Fs[ix(FsEy,j+1,k,l)]*wj1kl+Fs[ix(FsEy,j+1,k+1,l)]*wj1k1l+Fs[ix(FsEy,j,k+1,l)]*wjk1l+Fs[ix(FsEy,j,k,l+1)]*wjkl1+Fs[ix(FsEy,j+1,k,l+1)]*wj1kl1+Fs[ix(FsEy,j+1,k+1,l+1)]*wj1k1l1+Fs[ix(FsEy,j,k+1,l+1)]*wjk1l1);
	  forceZ=dpart[i][ii].q*(Fs[ix(FsEz,j,k,l)]*wjkl+Fs[ix(FsEz,j+1,k,l)]*wj1kl+Fs[ix(FsEz,j+1,k+1,l)]*wj1k1l+Fs[ix(FsEz,j,k+1,l)]*wjk1l+Fs[ix(FsEz,j,k,l+1)]*wjkl1+Fs[ix(FsEz,j+1,k,l+1)]*wj1kl1+Fs[ix(FsEz,j+1,k+1,l+1)]*wj1k1l1+Fs[ix(FsEz,j,k+1,l+1)]*wjk1l1);
		  
	  
//printf("charge %E\n",charge[dpart[i][ii].spec]);getchar();
	 
	  /*find perp and paral component*/
 	//  r=sqrt((x-dmass_centr_x[i])*(x-dmass_centr_x[i])+(y-dmass_centr_y[i])*(y-dmass_centr_y[i]));
	//  cosbeta=(y-dmass_centr_y[i])/r;
	//  sinbeta=(x-dmass_centr_x[i])/r;
	  /*if force par is negative then it is directed towards dust*/
	  /*if force perp is posit. then it is directed CW for Force along x direction and ACW for force along y direction */ 
	//  fpar=forceY*cosbeta+forceX*sinbeta; /*if negative -> push dust*/
	//  ftransX+=fpar*sinbeta;
	//  ftransY+=fpar*cosbeta;

	  //  fperp=forceY*sinbeta-forceX*cosbeta; /*if positive -> ACW*/	
	  /*find torque from this particle r x F */
	  // torque+=r*fperp;
	 // torque+=forceY*(x-dmass_centr_x[i])-forceX*(y-dmass_centr_y[i]);

	}
#ifdef MPI
      /*CHECK MPI HERE!*/
      /* sum the forces from different ranks*/
      /*check if it works*/
      MPI_Barrier(MPI_COMM_WORLD);
     // MPI_Reduce(&ftransX,&rftransX, 1, MPI_DOUBLE,MPI_SUM,0, MPI_COMM_WORLD);
     // MPI_Reduce(&ftransY,&rftransY, 1, MPI_DOUBLE,MPI_SUM,0, MPI_COMM_WORLD);
     // MPI_Reduce(&torque,&rtorque,1, MPI_DOUBLE,MPI_SUM,0, MPI_COMM_WORLD);
     // MPI_Barrier(MPI_COMM_WORLD);
      if(rank==0)
	{
	//  ftransX=rftransX;
	//  ftransY=rftransY;
	//  torque=rtorque;
	}
//	getchar();
      MPI_Barrier(MPI_COMM_WORLD);
     // MPI_Bcast(&ftransX,1,MPI_DOUBLE,0,MPI_COMM_WORLD);  
     // MPI_Bcast(&ftransY,1,MPI_DOUBLE,0,MPI_COMM_WORLD);  
     // MPI_Bcast(&torque,1,MPI_DOUBLE,0,MPI_COMM_WORLD);  
#endif
   //   force_chk=fopen("./data/force_chk.dat","a");
    //  fprintf(force_chk,"%d\t%E\t%E\t%E\t%E\n", i, t*dt*normtime, ftransX*normforce, ftransY*normforce, torque*normtorque);
    //  fclose(force_chk);
        ftransY=0.0;
        ftransX=0.0;
        ftransZ=0.0;	
		torque=0.0;
      //printf("DUST %d: moving Fx:%E Fy: %E Fz: %E  T:%E\n", i, ftransX,ftransY,ftransZ,torque); 
 //     getchar(); 
      /*MOVE DUST translational movement fprce already multipl with dt*/
	  //printf("mass %E\n", dmass[i]);
      dustaccx[i]=ftransX/dmass[i];
      dustaccy[i]=ftransY/dmass[i];
	  dustaccz[i]=ftransZ/dmass[i]; 
      //printf("mass %E\n", dmass[i]);
	  dustvxc[i]+=dustaccx[i];
      dustvyc[i]+=dustaccy[i];
      dustvzc[i]+=dustaccz[i];
	  transx=dustvxc[i]*dt;
      transy=dustvyc[i]*dt;
      transz=dustvzc[i]*dt;
	  /*move coordinates*/
      dmass_centr_x[i]+=transx;
      dmass_centr_y[i]+=transy;
      dmass_centr_z[i]+=transz;
	  //printf("befffore\n");
	  for(j=0; j<ncorners[i]; j++)
	   {
	    dustxdxold[i][j]=dustxdx[i][j];
	    dustydyold[i][j]=dustydy[i][j];
		dustzdzold[i][j]=dustzdz[i][j];
	  dustxdx[i][j]+=transx;
	  dustydy[i][j]+=transy;
	  dustzdz[i][j]+=transz;
	  dustx[i][j]+=transx/dx;
	  dusty[i][j]+=transy/dy;
	  dustz[i][j]+=transz/dz;
	  // printf("after moving corner %d %E %E\n", j, dustx[i][j], dusty[i][j]);
	  //  getchar();
      //printf("by if photons\n");	
	  /*photoelectrons unit*/
	  if(photons)
	    {
//	      dustxnormv[i][j]+=transx;
//	      dustynormv[i][j]+=transy;
//		  dustznormv[i][j]+=transz;
	    }
	   }
      /*angular movement torqe already multipl with dt*/
      duste[i]=torque/dmomI[i];
      dustomega[i]+=duste[i];
      cosphi=cos(dustomega[i]*dt);
      sinphi=sin(dustomega[i]*dt);
      //   printf("cosphi %E\t sinphi %E dt %E\n", cosphi, sinphi, dt);
   
      for(j=0; j<ncorners[i]; j++)
	{	
	  //	  printf("before rotatin corner %d centr %E %E and corn %E %E\n", j, dmass_centr_x[i]/dx, dmass_centr_y[i]/dy, dustx[i][j], dusty[i][j]); //getchar();
	  /*rotate*/
//	  tempx=cosphi*(dustxdx[i][j]-dmass_centr_x[i])-sinphi*(dustydy[i][j]-dmass_centr_y[i]);
//	  tempy=sinphi*(dustxdx[i][j]-dmass_centr_x[i])+cosphi*(dustydy[i][j]-dmass_centr_y[i]);

//	  dustxdx[i][j]=(tempx+dmass_centr_x[i]);
//	  dustydy[i][j]=(tempy+dmass_centr_y[i]);
//	  dustx[i][j]=dustxdx[i][j]/dx;
//	  dusty[i][j]=dustydy[i][j]/dy;
	  /*and now find vx and vy total for each corner*/
//	  dustvx[i][j]=(dustxdx[i][j]-dustxdxold[i][j])/dt;
//	  dustvy[i][j]=(dustydy[i][j]-dustydyold[i][j])/dt;

	  /*photoelectrons unit*/
	  if(photons)
	    {
//	      tempx=cosphi*(dustxnormv[i][j])-sinphi*(dustynormv[i][j]);
//	      tempy=sinphi*(dustxnormv[i][j])+cosphi*(dustynormv[i][j]);
	    }
	//   printf("after rotatin corner %d centr %E %E and corn %1.15E %1.15E\n", j, dmass_centr_x[i], dmass_centr_y[i], dustxdx[i][j]*normx, dustydy[i][j]*normx); 	  
	  // printf("after rotatin corner not norm %d centr %E %E and corn %1.15E %1.15E\n", j, dmass_centr_x[i]/dx, dmass_centr_y[i]/dy, dustxdx[i][j], dustydy[i][j]); 
	}



      //    getchar();
      /*now calculate coefficients for dust hitting*/  
      for(j=0; j<ncorners[i]; j++)
	{
	  jp1=j+1;
	  if(j==(ncorners[i]-1))
	    {jp1=0; 	    }
	
//	  daa[i][j]=dustvx[i][j]*dustvy[i][jp1]-dustvy[i][j]*dustvx[i][jp1];
//	  dbb[i][j]=dustxdxold[i][j]*dustvy[i][jp1]+dustydyold[i][jp1]*dustvx[i][j]-dustydyold[i][j]*dustvx[i][jp1]-dustxdxold[i][jp1]*dustvy[i][j];
//	  dcc[i][j]=dustxdxold[i][j]*dustydyold[i][jp1]-dustydyold[i][j]*dustxdxold[i][jp1];
//	  daa1y[i][j]=dustvy[i][jp1]-dustvy[i][j];
//	  daa1x[i][j]=dustvx[i][j]-dustvx[i][jp1];
//	  dbb1x[i][j]=dustxdxold[i][j]-dustxdxold[i][jp1]; 
//	  dbb1y[i][j]=dustydyold[i][jp1]-dustydyold[i][j];
//	  dcc1x[i][j]=dustxdxold[i][jp1]-dustxdxold[i][j]; 
//	  dcc1y[i][j]=dustydyold[i][jp1]-dustydyold[i][j];

	  // printf("IN DUSTG kk %d daa %E\t dbb%E\t,dcc%E\n\n\ndcc1x,%E\t dcc1y %E\n", j, daa[i][j], dbb[i][j], dcc[i][j], dcc1x[i][j], dcc1y[i][j]);
	  //	  printf("j %d daa1x: %E\t, daa1y: %E\n dbb1x %E, dbb1y %E\n Wall %d",j, daa1x[i][j], daa1y[i][j], dbb1x[i][j], dbb1y[i][j], j);
	  //printf("jeszcze raz %d\n", j);
			  //	  getchar();
	}

      /*MOVE also particles on the dust surface*/
      for(j=0; j<dpartlast[i]; j++)
	{
//	  dpart[i][j].x+=transx;
//	  dpart[i][j].y+=transy;
	  /*rotate*/
//	  tempx=cosphi*(dpart[i][j].x-dmass_centr_x[i])-sinphi*(dpart[i][j].y-dmass_centr_y[i]);
//	  tempy=sinphi*(dpart[i][j].x-dmass_centr_x[i])+cosphi*(dpart[i][j].y-dmass_centr_y[i]);
//	  dpart[i][j].x=tempx+dmass_centr_x[i];
//	  dpart[i][j].y=tempy+dmass_centr_y[i];
   
	  if(dtype[i]==D_COND)
	    {
//	      csq[i][j].x1+=transx;
//	      csq[i][j].y1+=transy;/*rotate*/
//	      tempx=cosphi*(csq[i][j].x1-dmass_centr_x[i])-sinphi*(csq[i][j].y1-dmass_centr_y[i]);
//	      tempy=sinphi*(csq[i][j].x1-dmass_centr_x[i])+cosphi*(csq[i][j].y1-dmass_centr_y[i]);
//	      csq[i][j].x1=tempx+dmass_centr_x[i];
//	      csq[i][j].y1=tempy+dmass_centr_y[i];
	      
//	      csq[i][j].x2+=transx;
//	      csq[i][j].y2+=transy;/*rotate*/
//	      tempx=cosphi*(csq[i][j].x2-dmass_centr_x[i])-sinphi*(csq[i][j].y2-dmass_centr_y[i]);
///	      tempy=sinphi*(csq[i][j].x2-dmass_centr_x[i])+cosphi*(csq[i][j].y2-dmass_centr_y[i]);
//	      csq[i][j].x2=tempx+dmass_centr_x[i];
//	      csq[i][j].y2=tempy+dmass_centr_y[i];
	      
//	      csq[i][j].x3+=transx;
//	      csq[i][j].y3+=transy;/*rotate*/
//	      tempx=cosphi*(csq[i][j].x3-dmass_centr_x[i])-sinphi*(csq[i][j].y3-dmass_centr_y[i]);
//	      tempy=sinphi*(csq[i][j].x3-dmass_centr_x[i])+cosphi*(csq[i][j].y3-dmass_centr_y[i]);
//	      csq[i][j].x3=tempx+dmass_centr_x[i];
//	      csq[i][j].y3=tempy+dmass_centr_y[i];
	      
//	      csq[i][j].x4+=transx;
//	      csq[i][j].y4+=transy;/*rotate*/
//	      tempx=cosphi*(csq[i][j].x4-dmass_centr_x[i])-sinphi*(csq[i][j].y4-dmass_centr_y[i]);
//	      tempy=sinphi*(csq[i][j].x4-dmass_centr_x[i])+cosphi*(csq[i][j].y4-dmass_centr_y[i]);
//	      csq[i][j].x4=tempx+dmass_centr_x[i];
//	      csq[i][j].y4=tempy+dmass_centr_y[i];
	    }
	}
//      findabv(i);
    }
}

/*******************calculate all static parameters for dust grains***/
void calculate_staticparameters(int arc, char *arv[])
{
  int i,j;
  printf("I am assigning memory for dust particles\n");   
  memoryduststatic(noofdusts);
//  memorydpart(noofdusts, DPSEGM);
 
 
 
 // printf("Triangulating\n");  
//  d_polygon(arc, arv);  /***check it here***/
  printf("Calculating centre of mass and momentum of inertia\n");
  d_centreofmass_and_momI(); 
  
  /*find the r2 distance between vertex and centerofmass*/
  /*I do not use that!!!!*/
//  for(i=0; i<noofdusts; i++)
 //   for(j=0; j<ncorners[i]; j++)
  //    {
//	dr2v2[i][j]=(dx*dustx[i][j]-dmass_centr_x[i])*(dx*dustx[i][j]-dmass_centr_x[i])+(dy*dusty[i][j]-dmass_centr_y[i])*(dy*dusty[i][j]-dmass_centr_y[i]);
 //     }
  printf("Calculated static parameters\n");
}

/***************************memory dpart*****************/
void memorydpart(int no, int dmax)
{
  // int i;
  /*I initialize dust particles here*/
  dpart=(d_particle **)malloc(no*sizeof(d_particle *));
  drho=(d_rho **)malloc(no*sizeof(d_rho *));

#ifdef MPI
  rdpart=(d_particle **)malloc(no*sizeof(d_particle *));
#endif
  /*I change now to virtual particles for insulators*/
  /*
  for(i=0; i<no; i++)
    if(dtype[i]!=D_COND)
     {
      dpart[i]=(d_particle *)malloc(dmax*sizeof(d_particle));
      dpartlast[i]=0; //assign zero to each point.var.
      dpartmax[i]=dmax; //this will show when to reallocate memory 
    }
  if (dpartlast == NULL || dpartmax == NULL || dpart==NULL)
    { 
      printf("Error: Out of Memory rank:%d\n", rank);
      exit(1);
    }
  */
}
/*****************memory static dust parameters*************/
 void memoryduststatic(int no)
  {
  int i;
  int nholes=0; /*need to be changed if holes included*/
  /***allocate mem for static elements***/  
  dmass_centr_x=(double *)malloc(no*sizeof(double));
  dmass_centr_y=(double *)malloc(no*sizeof(double));
  dr2v2=(double **)malloc(no*sizeof(double *));
  dmomI=(double *)malloc(no*sizeof(double));
  dmass=(double *)malloc(no*sizeof(double));
  nooftriangles=(int *)malloc(no*sizeof(int));
  dtrian=(dtriangle **)malloc(no*sizeof(dtriangle *));
  for(i=0; i<no; i++)
    {
      nooftriangles[i]=(ncorners[i]-2)+(2*nholes);
      //      printf("ncorners %d nholes %d\n", ncorners[i], nholes);
      dtrian[i]=(dtriangle *)malloc(nooftriangles[i]*sizeof(dtriangle));
      dr2v2[i]=(double *)malloc(ncorners[i]*sizeof(double));
    }
  if (dtrian == NULL || dr2v2 == NULL)
    { 
      printf("Error: Out of Memory rank:%d\n", rank);
      exit(1);
    }



	}








/*this function needs to be called only once -> assuming that dust shape/mass doesnot change*/
void d_centreofmass_and_momI()
{
  /*find the centre of mass for each dust grain in units of grid spacing*/
  int i,j;
  int pt0,pt1,pt2;
  for(i=0; i<noofdusts; i++)
    {
      dmass_centr_x[i]=0;
      dmass_centr_y[i]=0;
      dmass[i]=0.0;
      dmomI[i]=0.0;

      for(j=0; j<nooftriangles[i]; j++)
	{
	  /*find area of triangles -> mass*/
	  pt0=dtrian[i][j].pt1;
	  pt1=dtrian[i][j].pt2;
	  pt2=dtrian[i][j].pt3;
	  dtrian[i][j].area=0.5*(((dustx[i][pt1]-dustx[i][pt0])*(dusty[i][pt2]-dusty[i][pt0]))-((dustx[i][pt2]-dustx[i][pt0])*(dusty[i][pt1]-dusty[i][pt0])));
	  /* mass in real units * dx*dy */
	  printf("dx dy %E %E\n", dx, dy);
	  dtrian[i][j].mass=dtrian[i][j].area*dustrho[i]*dx*dy;
	  //mormalize it we have rho -> kg/m3 
	  dtrian[i][j].mass*=normx*normx/mass[0];
	  //normmass;
	  printf("DRHO for %d is %E\n",i, dustrho[i]); //getchar();
	  /*find centre of mass*/
	  /*find median of a triangle*/
	  dtrian[i][j].tcx=(dustx[i][pt0]+dustx[i][pt1]+dustx[i][pt2])/3.0;
	  dtrian[i][j].tcy=(dusty[i][pt0]+dusty[i][pt1]+dusty[i][pt2])/3.0;
	  printf("triangle center %E %E\n",  dtrian[i][j].tcx, dtrian[i][j].tcy);
	  dmass[i]+=dtrian[i][j].mass;
	  dmass_centr_x[i]+=dtrian[i][j].mass*dtrian[i][j].tcx;
	  dmass_centr_y[i]+=dtrian[i][j].mass*dtrian[i][j].tcy;
	  /*build moment of inertia -> axis in centre of mass*/
	  /*mom of inertia in real coordinates*/
	  dmomI[i]+=dtrian[i][j].mass*(dx*dx*dtrian[i][j].tcx*dtrian[i][j].tcx+dtrian[i][j].tcy*dtrian[i][j].tcy*dy*dy);
	}
      /*get the coordimates for centre of mass in real coordinates*/
      dmass_centr_x[i]=dmass_centr_x[i]*dx/dmass[i];
      dmass_centr_y[i]=dmass_centr_y[i]*dy/dmass[i];    
      printf("Dust %d: centre of mass %E %E and MASS %E MOM %E\n",i,dmass_centr_x[i]/dx,dmass_centr_y[i]/dy, dmass[i], dmomI[i] );
      
      printf("Dust %d: centre of mass %E %E\n",i,dmass_centr_x[i]/dx+dustcx[i],dmass_centr_y[i]/dy+dustcy[i] );
      //      getchar();
    }
}

void redistribute(void)
{
  int i,ii;
  double total, redist;
  for(i=0;i<noofdusts;i++)
    {
      total=0.0;
      redist=0.0;
      for(ii=0; ii<dpartlast[i]; ii++) 
	{
	  total+=dpart[i][ii].q;
	}
		//TESTTT
		//total=-normalcharge[0]/ratio;
		//printf("charge of the dust %E\n", total);
      redist=total/dpartlast[i];
      for(ii=0; ii<dpartlast[i]; ii++) 
	{
	  dpart[i][ii].q=redist;
	}
    }
}



/****First order, linear weighting, to calculate charge dens on grid points****/
void weightingdust1(int ko)
{
	
	int i,ii; //particles
	int j,k,l,index; //grid
	double q;
	double x,y,z;	
	
  // for MD list and tests, clearing the d_globallist for each dust
	for(i=0; i<ngx*ngy*ngz; i++)
	   d_globallist[i]=0;	
	
  //dpart[0][0].q=-(10000/ratio)*(1/(dens[0]*dV*normx*normx*normx));
//*(-1)*chargeandnorm[0]/(dV);
 // dpart[0][0].q=1000*normalcharge[0];
  printf("weighting dust\n");
  redistribute();
  //printf("CONDUCTING DUST!\n");



  for(i=0;i<noofdusts;i++)
  {
	  
// for MD list and tests
	  d_localmax[i]=0;
	  
//	  printf("number of dust cells %d\n", dpartlast[i]); getchar();
	  
    if(dtype[i]==D_COND)
      if(ko==0)
        chargeoncond(i);
    for(ii=0; ii<dpartlast[i]; ii++) 
      {
	//I HOPE ALL IS OK HERE
	//minus beacause of ... chargeandnorm of 0 is negative...
		  //dV  is to facilitate the weithting process
	q=dpart[i][ii].q*(-1)*chargeandnorm[0]/(dV);
	//  printf("Q for ii %d is %E, check %d is %E\n", ii, q, dpart[i][ii].spec, chargeandnorm[dpart[i][ii].spec]);
	  j=dpart[i][ii].x/dx; //take the integer
	  k=dpart[i][ii].y/dy;
	  l=dpart[i][ii].z/dz;  
	  x=dpart[i][ii].x-j*dx;
	  y=dpart[i][ii].y-k*dy;
	  z=dpart[i][ii].z-l*dz;
	  rho[ix(0,j,k,l)]+=q*(dx-x)*(dy-y)*(dz-z);
	  rho[ix(0,j+1,k,l)]+=q*x*(dy-y)*(dz-z);
	  rho[ix(0,j+1,k+1,l)]+=q*x*y*(dz-z);
	  rho[ix(0,j,k+1,l)]+=q*(dx-x)*y*(dz-z);
	  rho[ix(0,j,k,l+1)]+=q*(dx-x)*(dy-y)*z;
	  rho[ix(0,j+1,k,l+1)]+=q*x*(dy-y)*z;
	  rho[ix(0,j+1,k+1,l+1)]+=q*x*y*z;
	  rho[ix(0,j,k+1,l+1)]+=q*(dx-x)*y*z; 
	//  printf("Weighting dust %d %d %d\n", j,k,l);
	//	  	printf("parameters irho ngrid %E %d of %d \n", irho[6][1][1][1], ii, dpartlast[i]);	
	//for MD part and tests
	
	  int llist, checkk;
	  checkk=-1;
	  index=ix(0,j,k,l);
		
	  d_globallist[index]=1;
	  if(j-1 > 0)
    	 d_globallist[ix(0,j-1,k,l)]=1;	  
	  if(j+1 < ngx)	
		 d_globallist[ix(0,j+1,k,l)]=1;	
	  if(k-1 > 0)
		 d_globallist[ix(0,j,k-1,l)]=1;	
	  if(k+1 < ngy)
		 d_globallist[ix(0,j,k+1,l)]=1;	
	  if(l-1 > 0)
		 d_globallist[ix(0,j,k,l-1)]=1;	
	  if(l+1 < ngz)
		 d_globallist[ix(0,j,k,l+1)]=1;	
		  
		 // printf("before in dust weightilng j %d k %d l %d\n", j, k, l);
		//  getchar();
      for(llist=0; llist<d_localmax[i]; llist++)
	   if(d_locallist[i][llist]==index)
	   {
		   checkk=1;
	   }
	  if(checkk==-1)
		{
		//assume max number of occupied cells by one dust = 50
			d_localmax[i]++;	
			if(d_localmax[i]>LIST_SIZE) 
			{
			printf("LIST_SIZE too small, program will crash in weithting dust\n");
			exit(1);
			}
	        d_locallist[i][d_localmax[i]-1]=index;	
		//	printf("in dust weightilng j %d k %d l %d max %d\n", j,k,l, d_localmax[i]);
		//	getchar();
		}
	  }   
  } 
  printf("rank %d finished weighting\n", rank);
}	

/***********check cond ***********************/
void checkcond(void)
{
/*
  int i,ii,z;
  int j,k;
  double pot,x,y;
  double wjk,wj1k, wjk1, wj1k1;
  
  for(i=0; i<noofdusts; i++)
    {
      for(ii=0; ii<dpartlast[i]; ii++)
	{
	  pot=0.0;
	  x=dpart[i][ii].x;
	  y=dpart[i][ii].y;
	  
	  j=dpart[i][ii].x/dx;
	  k=dpart[i][ii].x/dy;
	  x=dpart[i][ii].x-j*dx;
	  y=dpart[i][ii].x-k*dy;
	  wjk=(dx-x)*(dy-y);
	  wj1k=x*(dy-y);
	  wj1k1=x*y;
	  wjk1=(dx-x)*y;
	  
	  pot=(gp[j][k].phi*wjk+gp[j+1][k].phi*wj1k+gp[j][k+1].phi*wjk1+gp[j+1][k+1].phi*wj1k1)/dxdy;
	  z=0;
	  if(gp[j][k].marker>1)
	    {pot+=gp[j][k].phi;z++;}
	  if(gp[j+1][k].marker>1)
	    {pot+=gp[j+1][k].phi;z++;}
	  if(gp[j][k+1].marker>1)
	    {pot+=gp[j][k+1].phi;z++;}
	  if(gp[j+1][k+1].marker>1)
	    {pot+=gp[j+1][k+1].phi;z++;}
	  //if(z>0);	
	  pot/=(1.0*z);
	  
	  if(pot>0  && dphifl[i]<0)
	    if(dpart[i][ii].spec==1)
	      dpart[i][ii].spec=0;
	    else
	      dpart[i][ii].spec=1;
	  
	  if(pot<0 && dphifl[i]>0)
	    if(dpart[i][ii].spec==1)
	      dpart[i][ii].spec=0;
	    else
	      dpart[i][ii].spec=1;
	  
	  
	  if(dphifl[i]<0)
	    {
	      if(fabs(pot) > fabs(dphifl[i]))
		dpart[i][ii].q*=0.75;
	      else
		dpart[i][ii].q*=1.25;
	    }
	  else
	    {
	      if(fabs(pot) < fabs(dphifl[i]))
		dpart[i][ii].q*=0.75;
	      else
		dpart[i][ii].q*=1.25;
	    }
	  
	  
	  //printf("i %d and pot %E and z %d %E %E spec %d\n", ii, dpart[i][ii].q, z, pot, dphifl[i], dpart[i][ii].spec);
	}
    }
	*/
}


/***********find charge on conductor**********/
void chargeoncond(int i)
{
/*
int ii,j,k,bz,kk;
double x,y;
double wjk,wj1k, wjk1, wj1k1;
double phi1,phi2,phi3,phi4;
double q;

//h2=delta*delta;
bz=kk=0;
for(ii=0; ii<dpartlast[i]; ii++)
{   
   if(csq[i][ii].in1==1)
       phi1=dphifl[i];
	else{
	 x=csq[i][ii].x1;
     y=csq[i][ii].y1;
     //interpolate potential from the closest grid points,and locate particle:)
     j=csq[i][ii].x1/dx;
     k=csq[i][ii].y1/dy;
     x=csq[i][ii].x1-j*dx;
     y=csq[i][ii].y1-k*dy;
     wjk=(dx-x)*(dy-y);
     wj1k=x*(dy-y);
     wj1k1=x*y;
     wjk1=(dx-x)*y;
     phi1=(gp[j][k].phi*wjk+gp[j+1][k].phi*wj1k+gp[j][k+1].phi*wjk1+gp[j+1][k+1].phi*wj1k1)/dxdy;
	} 
   
   if(csq[i][ii].in2==1)
     phi2=dphifl[i];
   else{
     x=csq[i][ii].x2;
     y=csq[i][ii].y2;
     //interpolate potential from the closest grid points,and locate particle:)
     j=csq[i][ii].x2/dx;
     k=csq[i][ii].y2/dy;
     x=csq[i][ii].x2-j*dx;
     y=csq[i][ii].y2-k*dy;
     wjk=(dx-x)*(dy-y);
     wj1k=x*(dy-y);
     wj1k1=x*y;
     wjk1=(dx-x)*y;
     phi2=(gp[j][k].phi*wjk+gp[j+1][k].phi*wj1k+gp[j][k+1].phi*wjk1+gp[j+1][k+1].phi*wj1k1)/dxdy;
   }
   
   if(csq[i][ii].in3==1)
     phi3=dphifl[i];
   else{
     x=csq[i][ii].x3;
     y=csq[i][ii].y3;
     //interpolate potential from the closest grid points,and locate particle:)
     j=csq[i][ii].x3/dx;
     k=csq[i][ii].y3/dy;
     x=csq[i][ii].x3-j*dx;
     y=csq[i][ii].y3-k*dy;
     wjk=(dx-x)*(dy-y);
     wj1k=x*(dy-y);
     wj1k1=x*y;
     wjk1=(dx-x)*y;
     phi3=(gp[j][k].phi*wjk+gp[j+1][k].phi*wj1k+gp[j][k+1].phi*wjk1+gp[j+1][k+1].phi*wj1k1)/dxdy;
   }
   
   if(csq[i][ii].in4==1)
     phi4=dphifl[i];
   else{
     x=csq[i][ii].x4;
     y=csq[i][ii].y4;
     //interpolate potential from the closest grid points,and locate particle:)
     j=csq[i][ii].x4/dx;
     k=csq[i][ii].y4/dy;
     x=csq[i][ii].x4-j*dx;
     y=csq[i][ii].y4-k*dy;
     wjk=(dx-x)*(dy-y);
     wj1k=x*(dy-y);
     wj1k1=x*y;
     wjk1=(dx-x)*y;
     phi4=(gp[j][k].phi*wjk+gp[j+1][k].phi*wj1k+gp[j][k+1].phi*wjk1+gp[j+1][k+1].phi*wj1k1)/dxdy;
   }
   
   q=-((phi1+phi2+phi3+phi4)-4*dphifl[i]);
   
   
   if((ii==(dpartlast[i]-1)) || (ii==ccorner[i][kk])) //we have the last part before the corner
     {
       q/=2.0;
     }
   if(ii==ccorner[i][kk+1]-1) //we have the last part before the corner
     {
       if((kk+1)<ncorners[i])
	 kk++;
       q/=2.0;
     }

   printf("i %d, x %E ,y %E q %E,phi1 %E, %E, %E, %E and dphi %E\n", ii, dpart[i][ii].x, dpart[i][ii].y, q, phi1, phi2, phi3, phi4, dphifl[i]);
   
   if(q<0)
     dpart[i][ii].spec=0;
   else
     dpart[i][ii].spec=1;
   
   dpart[i][ii].q=fabs(q)/(numtasks); //on each node also on rank0 a little bit ;)
   
   //normalcharge[dpart[i][ii].spec]; do not use normcharge, because the charge is already normalized, i do have number of charges
}
*/
}

/****generate virtual particles on all particles****/
/****20 particles per grid length!*/
void virtpart(void)
{
  //check
  FILE *checkvp;
  checkvp=fopen("checkvp.dat", "w");

  int i,j,k,jp1;
  double circum=0.0,length;
  int jj,ii;
  ddelta=dx/20;
  rdrholast=0;
  
  unitvec=(vectorst **)malloc(noofdusts*sizeof(vectorst *));
  orthvec=(vectorst **)malloc(noofdusts*sizeof(vectorst *));
  drholast=(int *)malloc(noofdusts*sizeof(int));
  ccorner=(int **)malloc(noofdusts*sizeof(int *));
  printf("hopsssa noofdusts %d\n", noofdusts); //getchar();
  vipcorner=(int **)malloc(noofdusts*sizeof(int *));
  


  for(i=0; i<noofdusts; i++)
    { 
      drholast[i]=0;
      ccorner[i]=(int *)malloc(ncorners[i]*sizeof(int));	
      vipcorner[i]=(int *)malloc(ncorners[i]*sizeof(int));
      
      for(j=0; j<ncorners[i]; j++)
	{ 
	  /*periodicity*/
	  jp1=j+1;
	  if(j==ncorners[i]-1) jp1=0;
	  length=sqrt((dustx[i][jp1]-dustx[i][j])*(dustx[i][jp1]-dustx[i][j]) + (dusty[i][jp1]-dusty[i][j])*(dusty[i][jp1]-dusty[i][j]));
	  printf("length %E and %E %E\n", length, dustx[i][jp1],dusty[i][j]);
	 
	  circum+=length;
	  ccorner[i][j]=length*20+1;
	  printf("noof particles per segm %d\n", ccorner[i][j]);
	  vipcorner[i][j]=ccorner[i][j];
	  //  printf("ccorner assigned %d %d", j, ccorner[i][j]);
	  drholast[i]+=ccorner[i][j];
	  //find how many virtual particles you produce (20 per one grid length
	  // getchar();
	} 
      /*zero because of periodicity*/
      //      ccorner[i][0]=0;
      for(j=1; j<ncorners[i]; j++)
	{
	  ccorner[i][j]+=ccorner[i][j-1];
	  printf("ccorner assigned new %d %d\n", j, ccorner[i][j]);
	  //  getchar();
	}
      
      rdrholast+=drholast[i]; 
      drho[i]=(d_rho *)malloc(drholast[i]*sizeof(d_rho));
      if(drho[i]==NULL)
	{ 
	  printf("Error: Out of Memory rank:%d\n", rank);
	  exit(1);
	}
      
      //Do it for all dusts now :(
      //      if(dtype[i]==D_COND)
      {
	dpartlast[i]=dpartmax[i]=drholast[i];
	//and assign memory for virtual dust particles
	dpart[i]=(d_particle *)malloc(dpartlast[i]*sizeof(d_particle));

#ifdef MPI
	rdpart[i]=(d_particle *)malloc(dpartlast[i]*sizeof(d_particle));
#endif
	if(dpart[i]==NULL)
	  { 
	    printf("Error: Out of Memory rank:%d\n", rank);
	    exit(1);
	  }
      }
      printf("DUST %d circumference: %E\n", i, circum); 
      //create particles NOW
      
      unitvec[i]=(vectorst *)malloc(dpartlast[i]*sizeof(vectorst));
      orthvec[i]=(vectorst *)malloc(dpartlast[i]*sizeof(vectorst));
      /*
	jj=0;
	ii=1;
	drho[i][0].x=dustx[i][0]*dx;
	drho[i][0].y=dusty[i][0]*dy;
	if(dtype[i]==D_COND)
	{  
	dpart[i][0].x=dustx[i][0]*dx;
	dpart[i][0].y=dusty[i][0]*dy;
	}
	unitvec[i][0].x=unitvecseg[i][0].x;
	orthvec[i][0].x=unitvecseg[i][0].x;
	unitvec[i][0].y=unitvecseg[i][0].y;
	orthvec[i][0].y=unitvecseg[i][0].y;
      */
      jj=0;
      ii=0;
      unitvec[i][0].x=unitvecseg[i][0].x;
      orthvec[i][0].x=orthvecseg[i][0].x;
      unitvec[i][0].y=unitvecseg[i][0].y;
      orthvec[i][0].y=orthvecseg[i][0].y;
      dpart[i][0].x=dustx[i][0]*dx;
      dpart[i][0].y=dusty[i][0]*dy;
      drho[i][0].x=dustx[i][0]*dx;
      drho[i][0].y=dusty[i][0]*dy;
      for(k=1; k<drholast[i]; k++)
	{
	  //	  printf("k, last %d jj %d %d ccorner jp1 %d\n", drholast[i],k, jj, ccorner[i][jj+1]);
	  if(k==ccorner[i][jj])
	    {
	      jj++;
	      ii++;
	      drho[i][k].x=dustx[i][ii]*dx;
	      drho[i][k].y=dusty[i][ii]*dy;
	      //   if(dtype[i]==D_COND)
	      {  
		dpart[i][k].x=dustx[i][ii]*dx;
		dpart[i][k].y=dusty[i][ii]*dy;
	      }
	      unitvec[i][k].x=unitvecseg[i][jj].x;
	      orthvec[i][k].x=orthvecseg[i][jj].x;
	      unitvec[i][k].y=unitvecseg[i][jj].y;
	      orthvec[i][k].y=orthvecseg[i][jj].y;
	      printf("unit at %d: x %E, y %E and orht x %E, y %E\n", k,unitvec[i][k].x, unitvec[i][k].y, orthvec[i][k].x, orthvec[i][k].y);
	      //getchar();	   
	      //ii++;	    
	    }
	  else
	    {
	      drho[i][k].x=drho[i][k-1].x+unitvecseg[i][jj].x*ddelta;	
	      drho[i][k].y=drho[i][k-1].y+unitvecseg[i][jj].y*ddelta;	
	      //if(dtype[i]==D_COND)
	      {  	
		dpart[i][k].x=dpart[i][k-1].x+unitvecseg[i][jj].x*ddelta;	
		dpart[i][k].y=dpart[i][k-1].y+unitvecseg[i][jj].y*ddelta;	
		//printf("%d x %E,y %E\n",k, dpart[i][k].x, dpart[i][k].y);
	      }
	      unitvec[i][k].x=unitvecseg[i][jj].x;
	      orthvec[i][k].x=orthvecseg[i][jj].x;
	      unitvec[i][k].y=unitvecseg[i][jj].y;
	      orthvec[i][k].y=orthvecseg[i][jj].y;
	      //printf(unit (%E,%E) orth (%E,%E)\n", )
	    }	
	  fprintf(checkvp, "%E\t%E\n", dpart[i][k].x, dpart[i][k].y); 
	}
    }
  printf("WOJTEK\n");
  rdrho=(double *)malloc(rdrholast*sizeof(double));
  tordrho=(double *)malloc(rdrholast*sizeof(double));
  printf("Finished virtpart()\n");
  // getchar();
  for(i=0; i<noofdusts; i++)
    for(j=0; j<ncorners[i]; j++)
      {
       	printf("Dust %d nvirt per segm %d  is: %d vert %d\n", i, j,vipcorner[i][j], dustv[i][j]);
      }
  //  getchar();
  fclose(checkvp);

}

/*****orth and normal vectors on dust grain segments***/
/*****all dusts****/
void ortnormvec(void)
{
/*
  int i,j,jp1;
  double length;
  
  unitvecseg=(vectorst **)malloc(noofdusts*sizeof(vectorst *));
  orthvecseg=(vectorst **)malloc(noofdusts*sizeof(vectorst *));
  

  if(photons)
    {
      dustxnormv=(double **)malloc(noofdusts*sizeof(double *));
      dustynormv=(double **)malloc(noofdusts*sizeof(double *));
      for(i=0; i<noofdusts; i++)   
	{
	  dustxnormv[i]=(double *)malloc(ncorners[i]*sizeof(double));
	  dustynormv[i]=(double *)malloc(ncorners[i]*sizeof(double));	
	}
    }

  for(i=0; i<noofdusts; i++)
    {
      //allocate memory
      unitvecseg[i]=(vectorst *)malloc(ncorners[i]*sizeof(vectorst));
      orthvecseg[i]=(vectorst *)malloc(ncorners[i]*sizeof(vectorst));
      
      //calculate circumference of the particle in grid spacing units
      for(j=0; j<ncorners[i]; j++)
	{ 
	  jp1=j+1;
	  if(j==(ncorners[i]-1)) jp1=0;
	  
	  length=sqrt((dustx[i][jp1]-dustx[i][j])*(dustx[i][jp1]-dustx[i][j]) + (dusty[i][jp1]-dusty[i][j])*(dusty[i][jp1]-dusty[i][j]));
	  unitvecseg[i][j].x=(dustx[i][jp1]-dustx[i][j])/length;
	  unitvecseg[i][j].y=(dusty[i][jp1]-dusty[i][j])/length;
	  printf("Dust %d Unitvector %d is %E %Elength %E (%E, %E), and (%E, %E)\n",i,j, unitvecseg[i][j].x,unitvecseg[i][j].y, length, dustx[i][jp1], dusty[i][jp1], dustx[i][j], dusty[i][j]);
	  
	  
	  if(dustv[i][j]==0)
	    {
	      orthvecseg[i][j].x=dusta[i][j]*sqrt(1/(dusta[i][j]*dusta[i][j]+1));
	      orthvecseg[i][j].y=-sqrt(1/(dusta[i][j]*dusta[i][j]+1));
	    }
	  else
	    {
	      orthvecseg[i][j].x=1;
	      orthvecseg[i][j].y=0;
	    }
	  if(photons)
	    {
	      dustxnormv[i][j]=orthvecseg[i][j].x;
	      dustynormv[i][j]=orthvecseg[i][j].y;
	      if(initpartcheck(0.1*dustxnormv[i][j]+0.5*(dustxdx[i][j]+dustxdx[i][jp1]), 0.1*dustynormv[i][j]+0.5*(dustydy[i][j]+dustydy[i][jp1]), ddelta*0.0000000001) % 2 != 0)
		{
		  dustxnormv[i][j]*=-1;
		  dustynormv[i][j]*=-1; 
		  printf("changed\n");
		}
	      printf("\n\n normal vector for segment %d is x %E y %E\n", j, dustxnormv[i][j], dustynormv[i][j]);
	      
	    }
	}
    }
	*/

  /*separate check for dustx/y normv*/

}

/******cond squares init*******/
void condsquares(void){
/*
  int i,j,k;
  double x,y;
  
//for(i=0; i<noofdusts; i++)
//  if(dtype[i]==D_COND)
//        {
//		noc++;
//	     }	
//condlut=(int *)malloc(noc*sizeof(int));
 csq=(condsq **)malloc(noofdusts*sizeof(condsq *));
 
 j=0;
 for(i=0; i<noofdusts; i++)
   if(dtype[i]==D_COND)
     {
       //		condlut[j]=i;
       //	j++;  
       //allocate memory for a given dust
       csq[i]=(condsq *)malloc(dpartlast[i]*sizeof(condsq));          
       //go through all virtual particles on a given dust
       //to calculate x,y and find where it is
       
       for(k=0; k<dpartlast[i]; k++)
	 {
	   x=dpart[i][k].x;
	   y=dpart[i][k].y;
	   //  printf("x and  y %E %E unitvec %E\n", x,y, unitvec[i][k].x);
	   //forward
	   csq[i][k].x1=x+unitvec[i][k].x*ddelta/2;
	   csq[i][k].y1=y+unitvec[i][k].y*ddelta/2;
	   //  printf("tutaj\n");
	   //right
	   csq[i][k].x2=x+orthvec[i][k].x*ddelta/2;
	   csq[i][k].y2=y+orthvec[i][k].y*ddelta/2;
	   //backward
	   csq[i][k].x3=x-unitvec[i][k].x*ddelta/2;
	   csq[i][k].y3=y-unitvec[i][k].y*ddelta/2;
	   //left
	   csq[i][k].x4=x-orthvec[i][k].x*ddelta/2;
	   csq[i][k].y4=y-orthvec[i][k].y*ddelta/2;
	   
	   csq[i][k].in1=csq[i][k].in2=csq[i][k].in3=csq[i][k].in4=0;	
	   //printf("x %E y %E, x1 %E y1 %E, x2 %E, y2 %E, x3 %E, y3 %E, x4 %E, y4 %E\n", x, y, csq[i][k].x1, csq[i][k].y1,csq[i][k].x2, csq[i][k].y2, csq[i][k].x3, csq[i][k].y3,csq[i][k].x4, csq[i][k].y4);
	   if(initpartcheck(csq[i][k].x1, csq[i][k].y1, ddelta/10000) % 2 == 1)
             csq[i][k].in1=1;	
	   // printf("done in1\n");
	   if(initpartcheck(csq[j][k].x2, csq[j][k].y2, ddelta/10000) % 2 == 1)
	     csq[i][k].in2=1;
	   // printf("done in2\n");
	   if(initpartcheck(csq[i][k].x3, csq[i][k].y3, ddelta/10000) % 2 == 1)
             csq[i][k].in3=1;
// printf("done in3\n");
	  if(initpartcheck(csq[i][k].x4, csq[i][k].y4, ddelta/10000) % 2 == 1)
             csq[i][k].in4=1;
		//-ddelta/4
//		printf("sq k %d, %d %d %d %d\n", k, csq[i][k].in1,csq[i][k].in2,csq[i][k].in3,csq[i][k].in4);

//getchar();
	 }
     }
 
 //now we do not need anymore unitvec and orthvec :)
 free(unitvec);
 free(orthvec);
 free(unitvecseg);
 free(orthvecseg);
 printf("Finished condsquares() \n");
}


		

void findnewpotentials(double interval, int collect, FILE *fpoint1)		
{
int i,j;
double totalcurrent; //double because ions may chave different charge -> calcul real ccurrent
for(i=0; i<noofdusts; i++)
{	
totalcurrent=0.0;
#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);	   
	MPI_Reduce(curr_av,rcurr_av,S,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
#endif
  if(rank==0)
    {
      for(j=0; j<S; j++)
	{
#ifdef MPI
	  totalcurrent+=charge[j]*rcurr_av[i][j];
#else
	  totalcurrent+=charge[j]*curr_av[i][j];
#endif
	}
      fprintf(fpoint1, "%E\t %E\n", dphifl[i]*normpot, (totalcurrent)/(collect*dt*normtime));
      printf("dust no %d ele on %E\t %ld\n",i, dphifl[i]*normpot, curr_av[i][0]);
      printf("dust no %d ion %E\t %ld\n", i,dphifl[i]*normpot, curr_av[i][1]);
      printf("dust no %d total current %E\t %E\n", i,dphifl[i]*normpot,(1.0*totalcurrent)/(collect*dt*normtime));   
	}   
    }
  if(totalcurrent < 0.0)
  dphifl[i]-=interval;
  else
  dphifl[i]+=interval;
  
  //clear curr_av
  for(j=0;j<S;j++)
  {
    curr_av[i][j]=0; 
	rcurr_av[i][j]=0;
   }
   */
}


//RETURNS THE TOTAL VOLUME OF ALL DUST IN THE UNITS OF GRID CELLS
double finddustvolume(int arc, char *arv[])
{
  double **dax, **day;
  double volume, totaldustvolume;
  double radius;
  double temp1, temp2;
  int dshapelocally;
  int ncount;
  char temp[80];
  int *nca;
  int i,j;
  FILE *dshape;
  /******open and read file**********/
  /******we can have many dust particles*******/
  /******create array for many dust particles********/

  dax=(double **)malloc((arc-1)*sizeof(double *));
  day=(double **)malloc((arc-1)*sizeof(double *));
  nca=(int *)malloc((arc-1)*sizeof(int));
  
  totaldustvolume=0.0;

  for(i=1; i<arc; i++)
    {
      dshape=my_file_open(arv[i], "r");
      
      while(fscanf(dshape, "%lf %lf %lf %d  %lf %lf %lf %lf %lf %lf", &temp1, &temp2, &temp2, &ncount, &temp1,&temp1, &temp1,&temp1,&temp1,&temp1) < 10)
	      {fscanf(dshape, "%s", temp);}
	  fscanf(dshape, "%s", temp); /*dust movement*/
	  fscanf(dshape, "%lf", &temp1);	  
      fscanf(dshape, "%s", temp); /*dust type*/
	  fscanf(dshape, "%lf", &temp1);	//printf("temp1 %s %E\n", temp, temp1);	  
      fscanf(dshape, "%s", temp); /*dust shape*/
	  fscanf(dshape, "%d", &dshapelocally); 
	
  
      printf("file %s open in dustvolume() dshapetype: %d \n", arv[i], dshapelocally);
      
	  if(dshapelocally==0) //spherical dust
	{
	fscanf(dshape, "%s", temp); //Radius given in grid points...
    fscanf(dshape, "%lf", &radius);
		volume=4*M_PI*radius*radius*radius/3.0;
		}	  
	if(dshapelocally==1) //polyhedral dust
		{
	  fscanf(dshape, "%s", temp); /*First contour*/
      fscanf(dshape, "%s", temp); /*Number of corners*/   
      fscanf(dshape, "%d", &nca[i-1]);     
      fscanf(dshape, "%s", temp); 
	  /********allocate memory for dust arrays**********/
	  
      for(j=0; j<nca[i-1]; j++)
		{
	  dax[i-1]=(double *)malloc((nca[i-1]+2)*sizeof(double));
	  day[i-1]=(double *)malloc((nca[i-1]+2)*sizeof(double));
		}
      //      printf("%f %f %d\n", dax[i-1][0], day[i-1][0], nca[i-1]);
 
      /*******read coordinates of the OUTER contour******/
/*
      for(j=0; j<nca[i-1]; j++)
		{	
	  fscanf(dshape, "%lf %lf", &dax[i-1][j], &day[i-1][j]);
	  //  printf("coord %f %f\n", dax[i-1][j], day[i-1][j]);
		}

      dax[i-1][nca[i-1]]=dax[i-1][0];
      day[i-1][nca[i-1]]=day[i-1][0];
      dax[i-1][nca[i-1]+1]=dax[i-1][1];
      day[i-1][nca[i-1]+1]=day[i-1][1];
	 
      *******move coordinates to x=x+1 and y=y+2*
         
	for(j=0; j<nca[i-1]+2; j++)
		{
	  dax[i-1][j]+=1;
	  day[i-1][j]+=2;
		}
      ******* we can now calculate the area*********
      area=0.0;

      for(j=1; j<=nca[i-1]; j++)
		{
	  area+= dax[i-1][j]*(day[i-1][j+1]-day[i-1][j-1]);
		}
	area=area/2.0;      
*/	  
	  free(dax);
      free(day);
      free(nca);
	  }
	  
      /*************close the dust input file **********/
      fclose(dshape);

	  totaldustvolume+=fabs(volume);
      printf("READY area %f and total area is %f -> file %s closed in dustarea()\n", fabs(volume), totaldustvolume, arv[i]);
      //  getchar();
    }
  return totaldustvolume;
}


/***calculate total area of dust grains in the sim. box.***/
double dustarea(int arc, char *arv[])
{
  double **dax, **day;
  double area, totaldustarea;
  double temp1, temp2;
  int ncount;
  char temp[80];
  int *nca;
  int i,j;
  FILE *dshape;
  /******open and read file**********/
  /******we can have many dust particles*******/
  /******create array for many dust particles********/

  dax=(double **)malloc((arc-1)*sizeof(double *));
  day=(double **)malloc((arc-1)*sizeof(double *));
  nca=(int *)malloc((arc-1)*sizeof(int));
  
  totaldustarea=0.0;

  for(i=1; i<arc; i++)
    {
      dshape=fopen(arv[i], "r");
      if(dshape== NULL)
	printf("Sorry, can not open the input file %s", arv[i]);
      while(fscanf(dshape, "%lf %lf %d  %lf %lf %lf %lf %lf", &temp1, &temp2, &ncount, &temp1,&temp1,&temp1,&temp1,&temp1) < 8)
	      {fscanf(dshape, "%s", temp);}
      fscanf(dshape, "%s", temp); /*dust type*/
	  fscanf(dshape, "%lf", &temp1);	//printf("temp1 %s %E\n", temp, temp1);
	    fscanf(dshape, "%s", temp); /*First contour*/
      fscanf(dshape, "%s", temp);    
      fscanf(dshape, "%d", &nca[i-1]);     
      fscanf(dshape, "%s", temp); 
  
      printf("file %s open in dustarea() \n", arv[i]);
      /********allocate memory for dust arrays**********/
      for(j=0; j<nca[i-1]; j++)
	{
	  dax[i-1]=(double *)malloc((nca[i-1]+2)*sizeof(double));
	  day[i-1]=(double *)malloc((nca[i-1]+2)*sizeof(double));
	}
      //      printf("%f %f %d\n", dax[i-1][0], day[i-1][0], nca[i-1]);
 
      /*******read coordinates of the OUTER contour******/

      for(j=0; j<nca[i-1]; j++)
	{	
	  fscanf(dshape, "%lf %lf", &dax[i-1][j], &day[i-1][j]);
	  //  printf("coord %f %f\n", dax[i-1][j], day[i-1][j]);
	}

      dax[i-1][nca[i-1]]=dax[i-1][0];
      day[i-1][nca[i-1]]=day[i-1][0];
      dax[i-1][nca[i-1]+1]=dax[i-1][1];
      day[i-1][nca[i-1]+1]=day[i-1][1];
	 
      /*******move coordinates to x=x+1 and y=y+2**/
         
      for(j=0; j<nca[i-1]+2; j++)
	{
	  dax[i-1][j]+=1;
	  day[i-1][j]+=2;
	}
      /******* we can now calculate the area**********/
      area=0.0;

      for(j=1; j<=nca[i-1]; j++)
	{
	  area+= dax[i-1][j]*(day[i-1][j+1]-day[i-1][j-1]);
	}
      area=area/2.0;
      /*************close the dust input file **********/
      fclose(dshape);
      totaldustarea+=fabs(area);
      printf("READY area %f and total area is %f -> file %s closed in dustarea()\n", fabs(area), totaldustarea, arv[i]);
      //  getchar();
    }
  free(dax);
  free(day);
  free(nca);

  return totaldustarea;
}


// Copyright 2000, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.


void average_current(void)
{
  int i,j;
  for(i=0; i<noofdusts; i++)		
    for(j=0; j<S; j++)
      {
	curr_av[i][j]+=current[i][j];
	current[i][j]=0;
      }
}


void drag_force_direct(double partvxnew, double partvynew, double partvznew, int particlespecie, int dno, double partxhit, double partyhit, double partzhit)
{
	
	double momentum;
	double dustcenter_x, dustcenter_y, dustcenter_z;
	double rx,ry,rz,rxu,ryu,rzu;
	double tang_x, tang_y, tang_z;
	double mymass, forcerzutnar;
	
	//partvynew=0;
	//partvznew=0;
	//partvxnew=-1;
	//partxhit=dustcxdx[dno]+dradiusdx[dno];
	//partyhit=dustcydx[dno];	
	//partzhit=dustczdx[dno];
	//mass is real * ratio
	//momentum in real units
	mymass=mass[particlespecie];
	momentum=mymass*sqrt(partvxnew*partvxnew+partvynew*partvynew+partvznew*partvznew)*normvel;
	dustcenter_x=dustcxdx[dno];
	dustcenter_y=dustcydx[dno];
	dustcenter_z=dustczdx[dno];
	
	rx=dustcenter_x-partxhit;
	ry=dustcenter_y-partyhit;
	rz=dustcenter_z-partzhit;
	//notrmal vector
	rxu=rx/dradiusdx[dno];
	ryu=ry/dradiusdx[dno];
	rzu=rz/dradiusdx[dno];
	//printf("\n tutu normal vec %E %E %E\n", rxu, ryu, rzu);
	
	//force acting on the centre
	//mass includes ratio, velocity is normalized -> this is in real values
	forcerzutnar=mymass*normvel*(rxu*partvxnew+ryu*partvynew+rzu*partvznew);
	//scalar becomes a vector
	drag_direct_x[dno]+=rxu*forcerzutnar;
	drag_direct_y[dno]+=ryu*forcerzutnar;
	drag_direct_z[dno]+=rzu*forcerzutnar;

	//printf("\n drag direct, %E %E %E\n", drag_direct_x[dno],  drag_direct_y[dno],  drag_direct_z[dno]);
	//rotational part	- normal to the surface
	tang_x=mymass*partvxnew*normvel-forcerzutnar*rxu;
	//mymass*partvxnew*normvel*(1-rxu);
	//*mymass*partvxnew*normvel;
	tang_y=mymass*partvynew*normvel-forcerzutnar*ryu;
	//*(1-ryu);
	//*mymass*partvynew*normvel;
	tang_z=mymass*partvznew*normvel-forcerzutnar*rzu;
	//*(1-rzu);
	//-rzu*mymass*partvznew*normvel;
	//tang_x=mymass*normvel*(ryu*partvznew-rzu*partvynew);
	//tang_y=-mymass*normvel*(rxu*partvznew-rzu*partvxnew);
	//tang_z=mymass*normvel*(rxu*partvynew-ryu*partvxnew);	
//	printf("\n tang  %E %E %E\n", tang_x, tang_y, tang_z);

	//around z axis
	//r_z_axis=sqrt((partxhit-dustcenter_x)*(partxhit-dustcenter_x)+(partyhit-dustcenter_y)*(partyhit-dustcenter_y));
	//tau around z axis
	
	//arm
	rx=-rx;
	ry=-ry;
	rz=-rz;
	
	//tau_z=rx*tang_y-ry*tang_x;
//	tau_y=-(rx*tang_z-rz*tang_x);
//	tau_x=ry*tang_z-rz*tang_y;
	
	if(rx>0)
	{
		//rotation around z axis z x1 point
		drot_z_x1[dno]+=rx*tang_y;
		drot_y_x1[dno]+=-rx*tang_z;  
	}
	else
	{
		drot_z_x2[dno]+=rx*tang_y;
        drot_y_x2[dno]+=-rx*tang_z;  
	}
	if(ry>0)
	{
		drot_z_y1[dno]+=-ry*tang_x;
		drot_x_y1[dno]+=ry*tang_z;
    }
	
	else		
	{
		drot_z_y2[dno]+=-ry*tang_x;
		drot_x_y2[dno]+=ry*tang_z;
	}
	if(rz>0)
	{
		drot_x_z1[dno]+=-rz*tang_y;
		drot_y_z1[dno]+=rz*tang_x;
	}	
	else		
	{
		drot_x_z2[dno]+=-rz*tang_y;
		drot_y_z2[dno]+=rz*tang_x;
	}
//	printf("\n drag direct rot, %E %E %E %E\n", drot_z_x1[dno], drot_z_x2[dno],drot_z_y1[dno], drot_z_y2[dno]);

	//I have a torque collected and redistributed on 6 extreme points drot_x,y,z 1,2[dno]
	//I have a force collected to:    drag_direct_x,y,z[dno]
	//since drag force is due to momentum change, I need to reduce it before printing ifdef MPI
	
}


/***electric component of the drag force***/
void drag_force_electric(void)
{
	int i,ii,j,k,l;
	double transx,transy, transz;
	double tempx, tempy, tempz;
	double sinbeta, cosbeta, r;
	double x,y,z,x1,y1,z1;
	double wjkl,wj1kl,wj1k1l,wjk1l,wjkl1,wj1kl1,wj1k1l1,wjk1l1; //weights
	double ftransX, ftransY, ftransZ, torque;
	double forcerzutnar;
	
	double fpar;
	double forceX, forceY, forceZ;
	double cosphi,sinphi;
	int jp1;
	
	double normforce=normmass*normx/(normtime*normtime);
	double normtorque=normforce*normx;
	
	double dustcenter_x, dustcenter_y, dustcenter_z;
	double rx,ry,rz,rxu,ryu,rzu;
	double tang_x, tang_y, tang_z;
	
	//  printf("in dust move\n noofdusts %d \n", noofdusts); getchar();
	
	for(i=0; i<noofdusts; i++)
    {
		
		//do i really need to nullify it here?
		ftransX=ftransY=ftransZ=torque=0.0;			
		elrot_z_x1[i]=elrot_z_x2[i]=elrot_y_x1[i]=elrot_y_x2[i]=elrot_z_y1[i]=elrot_z_y2[i]=elrot_y_z1[i]=elrot_y_z2[i]=elrot_x_y1[i]=elrot_x_y2[i]=elrot_x_z1[i]=elrot_x_z2[i]=0.0;
		drag_elect_x[i]=drag_elect_y[i]=drag_elect_z[i]=0.0;
		/*FIND FORCEs*/
		/*take all particles assigned to dust[i]-> weight force as usual*/				
		for(ii=0; ii<dpartlast[i]; ii++) 
		{  
			//no particles on rank 0 if MPI 			
			x=dpart[i][ii].x;
			y=dpart[i][ii].y;
			z=dpart[i][ii].z;
			
			//find forces from the closest grid points,and locate particle:)
			j=dpart[i][ii].x/dx;
			k=dpart[i][ii].y/dy;
			l=dpart[i][ii].z/dz;
			x1=dpart[i][ii].x-j*dx;
			y1=dpart[i][ii].y-k*dy;
			z1=dpart[i][ii].z-l*dz;
			/*find weights, E is multiplied already by dt/dxdy*/
			wjkl=(dx-x1)*(dy-y1)*(dz-z1);
			wj1kl=x1*(dy-y1)*(dz-z1);
			wj1k1l=x1*y1*(dz-z1);
			wjk1l=(dx-x1)*y1*(dz-z1);
			wjkl1=(dx-x1)*(dy-y1)*z1;
			wj1kl1=x1*(dy-y1)*z1;
			wj1k1l1=x1*y1*z1;
			wjk1l1=(dx-x1)*y1*z1;	  	  
	    	//force acting on part number ii
			forceX=dpart[i][ii].q*(Fs[ix(0,j,k,l)]*wjkl+Fs[ix(0,j+1,k,l)]*wj1kl+Fs[ix(0,j+1,k+1,l)]*wj1k1l+Fs[ix(0,j,k+1,l)]*wjk1l+Fs[ix(0,j,k,l+1)]*wjkl1+Fs[ix(0,j+1,k,l+1)]*wj1kl1+Fs[ix(0,j+1,k+1,l+1)]*wj1k1l1+Fs[ix(0,j,k+1,l+1)]*wjk1l1);
			forceY=dpart[i][ii].q*(Fs[ix(FsEy,j,k,l)]*wjkl+Fs[ix(FsEy,j+1,k,l)]*wj1kl+Fs[ix(FsEy,j+1,k+1,l)]*wj1k1l+Fs[ix(FsEy,j,k+1,l)]*wjk1l+Fs[ix(FsEy,j,k,l+1)]*wjkl1+Fs[ix(FsEy,j+1,k,l+1)]*wj1kl1+Fs[ix(FsEy,j+1,k+1,l+1)]*wj1k1l1+Fs[ix(FsEy,j,k+1,l+1)]*wjk1l1);
			forceZ=dpart[i][ii].q*(Fs[ix(FsEz,j,k,l)]*wjkl+Fs[ix(FsEz,j+1,k,l)]*wj1kl+Fs[ix(FsEz,j+1,k+1,l)]*wj1k1l+Fs[ix(FsEz,j,k+1,l)]*wjk1l+Fs[ix(FsEz,j,k,l+1)]*wjkl1+Fs[ix(FsEz,j+1,k,l+1)]*wj1kl1+Fs[ix(FsEz,j+1,k+1,l+1)]*wj1k1l1+Fs[ix(FsEz,j,k+1,l+1)]*wjk1l1);
			
			
			//mass is real * ratio
			//momentum in real units
			dustcenter_x=dustcxdx[i];
			dustcenter_y=dustcydx[i];
			dustcenter_z=dustczdx[i];
			
			rx=dustcenter_x-x;
			ry=dustcenter_y-y;
			rz=dustcenter_z-z;	
			//normal unit vectors in x,y,z
			rxu=rx/dradiusdx[i];
			ryu=ry/dradiusdx[i];
			rzu=rz/dradiusdx[i];
			
			//force acting on the centre
			forcerzutnar=(rxu*forceX+ryu*forceY+rzu*forceZ);

			drag_elect_x[i]+=rxu*forcerzutnar;
			drag_elect_y[i]+=ryu*forcerzutnar;
			drag_elect_z[i]+=rzu*forcerzutnar;
			
		
			//rotational part	- normal to the surface			
			//tang_x=(ryu*forceZ-rzu*forceY);
			//tang_y=-(rxu*forceZ-rzu*forceX);
			//tang_z=(rxu*forceY-ryu*forceX);	
			tang_x=forceX-rxu*forcerzutnar;
			tang_y=forceY-ryu*forcerzutnar;
			tang_z=forceZ-rzu*forcerzutnar;
			
			//around z axis
			//r_z_axis=sqrt((partxhit-dustcenter_x)*(partxhit-dustcenter_x)+(partyhit-dustcenter_y)*(partyhit-dustcenter_y));
			//tau around z axis
			rx=-rx;
			ry=-ry;
			rz=-rz;
			
			//tau_z=rx*tang_y-ry*tang_x;
			//tau_y=-(rx*tang_z-rz*tang_x);
			//tau_x=ry*tang_z-rz*tang_y;
			
			if(rx>0)
			{
				elrot_z_x1[i]+=rx*tang_y;
				elrot_y_x1[i]+=-rx*tang_z;  
			}
			else
			{
				elrot_z_x2[i]+=rx*tang_y;
				elrot_y_x2[i]+=-rx*tang_z;  
			}
			if(ry>0)
			{
				elrot_z_y1[i]+=-ry*tang_x;
				elrot_x_y1[i]+=ry*tang_z;
			}
			
			else		
			{
				elrot_z_y2[i]+=-ry*tang_x;
				elrot_x_y2[i]+=ry*tang_z;
			}
			if(rz>0)
			{
				elrot_x_z1[i]+=-rz*tang_y;
				elrot_y_z1[i]+=rz*tang_x;
			}	
			else		
			{
				elrot_x_z2[i]+=-rz*tang_y;
				elrot_y_z2[i]+=rz*tang_x;
			}
			
			//I have a torque collected and redistributed on 6 extreme points elrot_x,y,z 1,2[dno]
			//I have a force collected to:    drag_elect_x,y,z[dno]
			//since electric drag force is due to electric force acting on charge on each dust grain (differnt nodes)
			//I need to reduce it before printing ifdef MPI
		}
    }
}


void printdragforce(int timestep)
{
	FILE *fdrag;
	double nele,ndir;
	double meanradius;
	double dex[noofdusts],dey[noofdusts],dez[noofdusts],ddx[noofdusts],ddy[noofdusts],ddz[noofdusts];
	double exz1[noofdusts],exz2[noofdusts],exy1[noofdusts],exy2[noofdusts],eyx1[noofdusts],eyx2[noofdusts],eyz1[noofdusts],eyz2[noofdusts],ezx1[noofdusts],ezx2[noofdusts],ezy1[noofdusts],ezy2[noofdusts];
	double dxz1[noofdusts],dxz2[noofdusts],dxy1[noofdusts],dxy2[noofdusts],dyx1[noofdusts],dyx2[noofdusts],dyz1[noofdusts],dyz2[noofdusts],dzx1[noofdusts],dzx2[noofdusts],dzy1[noofdusts],dzy2[noofdusts];
	double extrax, extray, extraz;
	int dno;
	double rex[dno], rey[dno],rez[dno],rdx[dno], rdy[dno], rdz[dno];
#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);	
	for(dno=0; dno<noofdusts; dno++)
	{
		//printf("check1 before move rank %d dpart %E\n", rank, dpart[0][4].q);
		MPI_Reduce(&drag_elect_x[dno],&dex[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&drag_elect_y[dno],&dey[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&drag_elect_z[dno],&dez[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&drag_direct_x[dno],&ddx[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&drag_direct_y[dno],&ddy[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&drag_direct_z[dno],&ddz[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		//rotational electric		  
		MPI_Reduce(&elrot_x_z1[dno],&exz1[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&elrot_x_z2[dno],&exz2[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&elrot_x_y1[dno],&exy1[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&elrot_x_y2[dno],&exy2[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&elrot_y_x1[dno],&eyx1[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&elrot_y_x2[dno],&eyx2[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);  
		MPI_Reduce(&elrot_y_z1[dno],&eyz1[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&elrot_y_z2[dno],&eyz2[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);  
		MPI_Reduce(&elrot_z_x1[dno],&ezx1[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&elrot_z_x2[dno],&ezx2[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);  
		MPI_Reduce(&elrot_z_y1[dno],&ezy1[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&elrot_z_y2[dno],&ezy2[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); 
		//rotational direct  
		MPI_Reduce(&drot_x_z1[dno],&dxz1[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&drot_x_z2[dno],&dxz2[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&drot_x_y1[dno],&dxy1[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&drot_x_y2[dno],&dxy2[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&drot_y_x1[dno],&dyx1[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&drot_y_x2[dno],&dyx2[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);  
		MPI_Reduce(&drot_y_z1[dno],&dyz1[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&drot_y_z2[dno],&dyz2[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);  
		MPI_Reduce(&drot_z_x1[dno],&dzx1[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&drot_z_x2[dno],&dzx2[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);  
		MPI_Reduce(&drot_z_y1[dno],&dzy1[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&drot_z_y2[dno],&dzy2[dno],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); 
	} 
	MPI_Barrier(MPI_COMM_WORLD);	
	
#endif
	
#ifndef MPI
	for(dno=0; dno<noofdusts; dno++)
	{
		dex[dno]=drag_elect_x[dno];
		dey[dno]=drag_elect_y[dno];
		dez[dno]=drag_elect_z[dno];
		ddx[dno]=drag_direct_x[dno];
		ddy[dno]=drag_direct_y[dno];
		ddz[dno]=drag_direct_z[dno];
		exz1[dno]=elrot_x_z1[dno];
		exz2[dno]=elrot_x_z2[dno];
		exy1[dno]=elrot_x_y1[dno];
		exy2[dno]=elrot_x_y2[dno];
		eyx1[dno]=elrot_y_x1[dno];
		eyx2[dno]=elrot_y_x2[dno];  
		eyz1[dno]=elrot_y_z1[dno];
		eyz2[dno]=elrot_y_z2[dno];  
		ezx1[dno]=elrot_z_x1[dno];
		ezx2[dno]=elrot_z_x2[dno];  
		ezy1[dno]=elrot_z_y1[dno];
		ezy2[dno]=elrot_z_y2[dno];
		//direct
		dxz1[dno]=drot_x_z1[dno];
		dxz2[dno]=drot_x_z2[dno];
		dxy1[dno]=drot_x_y1[dno];
		dxy2[dno]=drot_x_y2[dno];
		dyx1[dno]=drot_y_x1[dno];
		dyx2[dno]=drot_y_x2[dno];  
		dyz1[dno]=drot_y_z1[dno];
		dyz2[dno]=drot_y_z2[dno];  
		dzx1[dno]=drot_z_x1[dno];
		dzx2[dno]=drot_z_x2[dno];  
		dzy1[dno]=drot_z_y1[dno];
		dzy2[dno]=drot_z_y2[dno];
	}			   			   
#endif
	if(rank==0)
	{ 
		for(dno=0; dno<noofdusts; dno++)
		{
			extrax=extray=extraz=0.0;
			
			//rotational
			rex[dno]=exz1[dno]+exz2[dno]+exy1[dno]+exy2[dno];
			rey[dno]=eyz1[dno]+eyz2[dno]+eyx1[dno]+eyx2[dno];
			rez[dno]=ezx1[dno]+ezx2[dno]+ezy1[dno]+ezy2[dno];
			
			//treat the rotational part 
			if(signof(ezx1[dno])!=signof(ezx2[dno]))
				extrax+=smaller_same_sign(ezx1[dno],-ezx2[dno]);
			if(signof(ezy1[dno])!=signof(ezy2[dno]))
				extray+=smaller_same_sign(ezy1[dno],-ezy2[dno]);
			if(signof(exz1[dno])!=signof(exz2[dno]))
				extraz+=smaller_same_sign(exz1[dno],-exz2[dno]);
			if(signof(exy1[dno])!=signof(exy2[dno]))
				extray+=smaller_same_sign(exy1[dno],-exy2[dno]);
			if(signof(eyx1[dno])!=signof(eyx2[dno]))
				extrax+=smaller_same_sign(eyx1[dno],-eyx2[dno]);
			if(signof(eyz1[dno])!=signof(eyz2[dno]))
				extraz+=smaller_same_sign(eyz1[dno],-eyz2[dno]);  
			
			//printf("signof %d\t%d\n", signof(-4), signof(4));
			//printf("smaller %E\n", smaller_same_sign(4,8));
			meanradius=M_PI*dradiusdx[dno]/4.0;   
			//factor 2 because we have two forces
			dex[dno]+=2*extrax/meanradius;
			dey[dno]+=2*extray/meanradius;
			dez[dno]+=2*extraz/meanradius;
			//printf("Drag force extra kick electric %E %E %E\n", extrax, extray, extraz);		   

			
			//rotational
			extrax=extray=extraz=0;		
			rdx[dno]=dxz1[dno]+dxz2[dno]+dxy1[dno]+dxy2[dno];
			rdy[dno]=dyz1[dno]+dyz2[dno]+dyx1[dno]+dyx2[dno];
			rdz[dno]=dzx1[dno]+dzx2[dno]+dzy1[dno]+dzy2[dno];
			printf("rot dir %E %E %E\n", rdx[dno], rdy[dno], rdz[dno]);
			//treat the rotational part 
			if(signof(dzx1[dno])!=signof(dzx2[dno]))
				extrax+=smaller_same_sign(dzx1[dno],-dzx2[dno]);
			if(signof(dzy1[dno])!=signof(dzy2[dno]))
				extray+=smaller_same_sign(dzy1[dno],-dzy2[dno]);
			if(signof(dxz1[dno])!=signof(dxz2[dno]))
				extraz+=smaller_same_sign(dxz1[dno],-dxz2[dno]);
			if(signof(dxy1[dno])!=signof(dxy2[dno]))
				extray+=smaller_same_sign(dxy1[dno],-dxy2[dno]);
			if(signof(dyx1[dno])!=signof(dyx2[dno]))
				extrax+=smaller_same_sign(dyx1[dno],-dyx2[dno]);
			if(signof(dyz1[dno])!=signof(dyz2[dno]))
				extraz+=smaller_same_sign(dyz1[dno],-dyz2[dno]);  
			
			meanradius=M_PI*dradiusdx[dno]/4.0;   
			ddx[dno]+=2*extrax/meanradius;
			ddy[dno]+=2*extray/meanradius;
			ddz[dno]+=2*extraz/meanradius;
		// 	printf("Drag force extra kick direct %E %E %E\n", extrax, extray, extraz);		   
			
			//normalize it all  
			//before it was only momentum change
			ndir=1.0/(dt*normtime);
			//before it was:
			nele=1.0*normEfield*ratio*normcharge;
			
		//	printf("time %E\n", timestep*dt*normtime);
		//	printf("Drag force direct part %E %E %E norm %E\n", ddx[dno],ddy[dno],ddz[dno], ndir);
		//	printf("Drag force electr part %E %E %E norm %E\n", dex[dno],dey[dno],dez[dno], nele); 
		//	printf("Drag force direct rot %E %E %E\n", rdx[dno],rdy[dno],rdz[dno]);
		//	printf("Drag force electr rot %E %E %E\n", rex[dno],rey[dno],rez[dno]);  
			
			
			fdrag=fopen("drag.txt", "a");
			fprintf(fdrag, "%E\t%d\t%E\t", timestep*dt*normtime, dno, meanradius*normx);
			fprintf(fdrag,"%E\t%E\t%E\t%E\t%E\t%E\t" , ddx[dno]*ndir, ddy[dno]*ndir, ddz[dno]*ndir, rdx[dno]*ndir*normx, rdy[dno]*ndir*normx, rdz[dno]*ndir*normx);
			fprintf(fdrag,"%E\t%E\t%E\t%E\t%E\t%E\n" , dex[dno]*nele, dey[dno]*nele, dez[dno]*nele, rex[dno]*nele*normx, rey[dno]*nele*normx, rez[dno]*nele*normx);
			fclose(fdrag);
			
		}	
	}
	//after printing clear the direct drag force collectors on all nodes
	
	for(dno=0; dno<noofdusts; dno++)
	{
	drot_z_x1[dno]=drot_z_x2[dno]=drot_y_x1[dno]=drot_y_x2[dno]=drot_z_y1[dno]=drot_z_y2[dno]=drot_y_z1[dno]=drot_y_z2[dno]=drot_x_y1[dno]=drot_x_y2[dno]=drot_x_z1[dno]=drot_x_z2[dno]=0.0;
	drag_direct_x[dno]=drag_direct_y[dno]=drag_direct_z[dno]=0.0;
	}
}

inline int signof(int a) 
{ 
return (a == 0) ? 0 : (a<0 ? -1 : 1); 
}


inline double smaller_same_sign(double a, double b) 
{ 
	if(a>0)
		return (a < b) ? a : b;
	else
		return (a >= b ) ? a : b;
}

