/* DUSTY PLASMA */
/* Various functions connected to the grid */
/* Author: Wojciech Jacek Miloch */
/* University of Oslo */
/* 2007 */

#include "const.h"
#include <math.h>
#include <stdlib.h>

/********Grid memory allocation ******************/
void memorygrid(void)
{
  int i;
  /*allocate memory for grid variables*/
  /*grid variables are stored in single arrays*/
  long int npoints=ngx*ngy*ngz; //all grid points
  phiMAX=phiavMAX=potconvMAX=PEtotalMAX=qdensMAX=npoints;
  PEMAX=npoints*S; //two species in one array
  PEMAXhalf=npoints;
  
  phi=dvecmem(0,phiMAX-1);  //potential
  phi_nodust=dvecmem(0,phiMAX-1);  //potential
  phiav=dvecmem(0,phiavMAX-1); //potential avareged in diagnostics
  PEtotal=dvecmem(0,PEtotalMAX-1); //total potential energy
  PE=dvecmem(0,PEMAX-1); //potential energy for each species
  qdens=dvecmem(0,qdensMAX-1); //charge density used in diagnostics (average)
  for(i=0; i<CONVTEST; i++)
    potconv[i]=dvecmem(0,npoints-1); //potential convegene in diagnostics
		  
  /*allocate memorty for vectors with grid data*/
  FsMAX=npoints*3*NOF; /*3 stands for 3D vector components of electric field*/
  FsEy=npoints; /*Ey offset*/
  FsEz=2*npoints; /*Ez offset*/
  Fs=dvecmem(0,FsMAX-1); //forces...
  Fs_nodust=dvecmem(0,FsMAX-1);
   //printf("FsEy %d, FsEz %d, FsMAX %d", FsEy, FsEz, FsMAX);
  //getchar();
  
  //Actual charge density
  rhoMAX=npoints; //same as rhoMAXhalf
  rhoMAXhalf=npoints; //we have no collected charge any more
  rho=dvecmem(0,rhoMAX-1); //on each node	
  rrho=dvecmem(0,rhoMAX-1); //reduced

  //each particle density (diagnostics)
  pdensMAX=npoints*S; /*S stands for el and prot*/
  pdens_off=npoints;
  pdens = dvecmem(0,pdensMAX-1); 
  rpdens = dvecmem(0,pdensMAX-1);
 
  //diagnostics: aerage velocity vectors
  vxvec=dvecmem(0,pdensMAX-1);	
  vyvec=dvecmem(0,pdensMAX-1);
  vzvec=dvecmem(0,pdensMAX-1);
  rvxvec=dvecmem(0,pdensMAX-1);
  rvyvec=dvecmem(0,pdensMAX-1);
  rvzvec=dvecmem(0,pdensMAX-1);
			
  KEMAX=npoints*S; /*S stands for el and prot*/
  KE_off=npoints;
  KE=dvecmem(0,KEMAX-1);
  rKE=dvecmem(0,KEMAX-1);
 
 /*CLEAR THE VECTORS INITIALLY*/
 for(i=0; i<rhoMAX; i++)
   {
     rho[i]=0.0;
     rrho[i]=0.0;
   }
   
   
  //create linkedlist  
if(((ngx-1)/4)==1)
  llngx=((ngx-1)/2)+1;
else
  llngx=((ngx-1)/4)+1;
  
if(((ngy-1)/4)==1)  
  llngy=((ngy-1)/2)+1;
else
  llngy=((ngy-1)/4)+1;
  
if(((ngz-1)/4)==1)    
  llngz=((ngz-1)/2)+1;
else
  llngz=((ngz-1)/4)+1;
  
  lldx=Lx/llngx;
  lldy=Ly/llngy;
  lldz=Lz/llngz;

  llsize=llngx*llngy*llngy;	  
  llmesh=ivecmem(0, llsize*S-1); //for each specie 
	
	
//PP(MD part)
	d_globallist=ivecmem(0,npoints-1);
}

/*********Free memory *************/
void memorygridfree(void)
{
int i;	
long int npoints=ngx*ngy*ngz;
  free_dvecmem(phi,0,phiMAX-1);
  free_dvecmem(phiav,0,phiMAX-1);
  free_dvecmem(PE,0,PEtotalMAX-1);
  free_dvecmem(qdens,0,qdensMAX-1);
  for(i=0; i<CONVTEST; i++)
      free_dvecmem(potconv[i],0,npoints-1);
//  free_dvecmem(potconv,0,);
  free_dvecmem(Fs,0,FsMAX-1);
  free_dvecmem(rho, 0, rhoMAX-1);
  free_dvecmem(pdens, 0, pdensMAX-1);
  free_dvecmem(KE,0,KEMAX);
  free_dvecmem(rrho,0,rhoMAX-1);
  free_dvecmem(rpdens,0,pdensMAX-1);
  free_dvecmem(rKE,0,KEMAX-1);
  free_ivecmem(llmesh, 0, llsize-1);
	free_ivecmem(d_globallist, 0, npoints-1);
}

/***********setting field values on the grid to zero before next run*******/
void cleargrid()
{
  int i;
 
  /*clear force vector*/
  for(i=0; i<FsMAX; i++)
    Fs[i]=0.0;

  for(i=0; i<rhoMAX; i++)
    rho[i]=0.0;
}

/***********setting field values on the grid to zero before next run*******/
void cleargrid2(void)
{
  int i;
  /*clear force vector*/
  for(i=0; i<FsMAX; i++)
    Fs[i]=0.0;

#ifdef INSU
  for(i=0; i<rhoMAXhalf; i++)
    rho[i]=rho[rhoMAXhalf+i];
#else
  for(i=0; i<rhoMAX; i++)
    rho[i]=0.0;
#endif
}

/****First order, linear weighting, to calculate charge dens on grid points****/
void weighting1(void)
{
  int i,ii; //particles
  int j,k,l; //grid
  double q;
  double x,y,z;
  for(i=0;i<S;i++)
    {
      q=chargeandnorm[i]/(dV);
      for(ii=0;ii<npart[i];ii++)
	{    	     	  
	  j=spec[i].part[ii].x/dx; //take the integer
	  k=spec[i].part[ii].y/dy;
	  l=spec[i].part[ii].z/dz;  
	  x=spec[i].part[ii].x-j*dx;
	  y=spec[i].part[ii].y-k*dy;
	  z=spec[i].part[ii].z-l*dz;
	  rho[ix(0,j,k,l)]+=q*(dx-x)*(dy-y)*(dz-z);
	  rho[ix(0,j+1,k,l)]+=q*x*(dy-y)*(dz-z);
	  rho[ix(0,j+1,k+1,l)]+=q*x*y*(dz-z);
	  rho[ix(0,j,k+1,l)]+=q*(dx-x)*y*(dz-z);
	  rho[ix(0,j,k,l+1)]+=q*(dx-x)*(dy-y)*z;
	  rho[ix(0,j+1,k,l+1)]+=q*x*(dy-y)*z;
	  rho[ix(0,j+1,k+1,l+1)]+=q*x*y*z;
	  rho[ix(0,j,k+1,l+1)]+=q*(dx-x)*y*z; 
	}
double test;
int ik,ij,il,cnt;
test=0.0;
cnt=0;
for(ik=1; ik<ngx-1; ik++)
for(ij=1; ij<ngy-1; ij++)
for(il=1; il<ngz-1; il++)
{
	test+=rho[ix(0,ik,ij,il)];	
	cnt++;

}	

	//printf("specie: %d on points 1: %E 2: %E 3: %E all %E cnt %d\n", i, rho[ix(0,10,10,10)], rho[ix(0,10,5,50)], rho[ix(0,5,10,10)], test/cnt, cnt);	 
    }
}

/************* generate boundaries ****************/
void gen_boundaries(void)
{
  int i,j,k;

  //Clear all grid points (note: we have Dirichlet BC, the potential is solved only on the inside points) 
  for(i=0; i< ngx; i++)
    for(j=0; j< ngy; j++)
    for(k=0; k< ngz; k++)
	      {     
				phi[ix(0,i,j,k)]=0.0;
			}
  /*clear vectors: Fs and rho*/
  for(i=0; i< FsMAX; i++)
    Fs[i] = 0.0;
  for(i=0; i< rhoMAX; i++)
    rho[i] = 0.0;
  
   //Distinguish boundaries
//  for(j=0; j<ngy; j++)
   // for(k=0; k<ngz; k++)
     {
	//    phi[ix(0,0,j,k)]=0.0;
	 //   phi[ix(0,ngx-1,j,k)]=0.0;
	 }

//  for(i=0; i<ngx; i++)
 //   for(k=0; k<ngz; k++)
     {
	//    phi[ix(0,i,0,k)]=0.0;
	//    phi[ix(0,i,ngy-1,k)]=0.0;
	 }

 // for(i=0; i<ngx; i++)
  //  for(j=0; j<ngy; j++)
     {
	//    phi[ix(0,i,j,0)]=0.0;
	//    phi[ix(0,i,j,ngz-1)]=0.0;
	 }
}

/*************generate probe -> new version***********************/
void gen_dust3D(int arc, char *arv[])
{
  /*this is the new version of the dust particle generator*/
  /*revised june 2007*/

  int i,j=0;
  int ct;
  //int ncorners;
 // double xmin, xmax, ymin, ymax;
  char temp[80];
  int ncont;
  int numberofpoints;
  FILE *file_dshape;
  /******open and read file**********/
  /******we can have many dust particles*******/
  /******create array for many dust particles********/
  noofdusts=arc-1;
  memorydust1_3D(noofdusts);
  cond_present=0;

   for(i=0; i<noofdusts; i++)
    {
      dustaccx[i]=0.0;
      dustaccy[i]=0.0;
	  dustaccz[i]=0.0;
      dustvxc[i]=0.0;
      dustvyc[i]=0.0;
      dustvzc[i]=0.0;
	  duste[i]=0.0;
      dustomega[i]=0.0;
	  dradius[i]=0.0; //just in case make it zero
	  dmass[i]=0.0;
	  dmomI[i]=0.0;
	  dmass_centr_x[i]=0.0;
	  dmass_centr_y[i]=0.0;
	  dmass_centr_z[i]=0.0;
    }

  for(i=1; i<arc; i++)
    {
      file_dshape=my_file_open(arv[i], "r");
     
      while(fscanf(file_dshape, "%lf %lf %lf %d %lf %lf %lf %lf %lf %lf", &dustcx[i-1], &dustcy[i-1], &dustcz[i-1], &ncont, &dustrho[i-1], &dustvxc[i-1], &dustvyc[i-1], &dustvzc[i-1], &dustomega[i-1], &dustworkfunct[i-1]) < 10)	
	     { fscanf(file_dshape, "%s", temp);}
      //	    printf("%s\n ", temp);}
	  
//	  printf("DUST PARAMETERS: %lf %lf %lf %d %lf %lf %lf %lf %lf %lf\n", dustcx[i-1], dustcy[i-1], dustcz[i-1], ncont, dustrho[i-1], dustvxc[i-1], dustvyc[i-1], dustvzc[i-1], dustomega[i-1], dustworkfunct[i-1]);	
//      printf("READ %f %f %d RHO %f workfunct %E\n", dustcx[i-1], dustcy[i-1],ncont, dustrho[i-1],dustworkfunct[i-1]);
      fscanf(file_dshape, "%s", temp); /*dust move*/
      fscanf(file_dshape, "%d", &dmove[i-1]);
	  fscanf(file_dshape, "%s", temp); /*dust type*/
      fscanf(file_dshape, "%d", &dtype[i-1]);
      fscanf(file_dshape, "%s", temp); /*dust shape*/
	  fscanf(file_dshape, "%d", &dshape[i-1]);
	  printf("file %s open in gen_dust3D() \n", arv[i]);
	    
	if(dshape[i-1]==0) //spherical dust
	  {
    	fscanf(file_dshape, "%s", temp); //read radius
        fscanf(file_dshape, "%lf", &dradius[i-1]);
		//now generate dust;
		//allocate memory for points on the dust
		//dradius[i-1]=sqrt(8);
    	numberofpoints=(int)(dradius[i-1]*M_PI*4*10);
		printf("numberofpoints on the dust object %d\n", numberofpoints);
		  getchar();
		memorydust2_3D(i-1,numberofpoints);
		//create points on the dust
		
        points_on_sphere(i-1,numberofpoints);	
      
		//move the points together with the centre
		dustcxdx[i-1]=dustcx[i-1]*dx;
		dustcydx[i-1]=dustcy[i-1]*dx;
		dustczdx[i-1]=dustcz[i-1]*dx;
		dradiusdx[i-1]=dradius[i-1]*dx;
		printf("dustcentre %E %E %E\n", dustcx[i-1], dustcy[i-1], dustcz[i-1]);
		printf("DDDDUUUSTT CENTRE %E %E %E %E %E %d\n", dustcxdx[i-1],  dustcydx[i-1],  dustczdx[i-1], dx, Lx, ngx);
		double normfact=normx*normx*normvel*normtime;
		printf("noof hits per cross section %E test %E dt: %E\n", M_PI*dradiusdx[i-1]*dradiusdx[i-1]*normx*normx, dens[1], M_PI*dradiusdx[i-1]*dradiusdx[i-1]*dt*vdriftx[1]*dens[1]*normfact);
printf("test2: %E\n", M_PI*dradiusdx[i-1]*dradiusdx[i-1]*dt*vdriftx[1]*dens[1]*(normx*normx*normx));

//getchar();
//	  printf("here\n");	
		for(ct=0; ct<numberofpoints; ct++)
		{
        dustxdx[i-1][ct]+=dustcxdx[i-1];
		dustydy[i-1][ct]+=dustcydx[i-1];
		dustzdz[i-1][ct]+=dustczdx[i-1];

		//IMPORTANT I AM OPERATING ON DPART STRUCTURE!!!
		dpart[i-1][ct].x=dustxdx[i-1][ct];
		dpart[i-1][ct].y=dustydy[i-1][ct];
		dpart[i-1][ct].z=dustzdz[i-1][ct];
		dpart[i-1][ct].q=0.0;	
		dpartq[i-1][ct]=0.0;	
//		printf("points %E %E %E test %E\n",dpart[i-1][ct].x/dx,dpart[i-1][ct].y/dx,dpart[i-1][ct].z/dx, dustcx[i-1]);
		}
		
		dpartlast[i-1]=numberofpoints;
		dpartmax[i-1]=numberofpoints;
		
		printf("Rank %d: Created %d points on the dust no. %d\n", rank, numberofpoints, i-1);								
	  }
	  
	  if(dshape[i-1]==1) //polyhedral dust
		{
	//		fscanf(file_dshape, "%s", temp); //read coordinates
	//		fscanf(file_dshape, "%d", &ncorners[i-1]);
	//		fscanf(file_dshape, "%s", temp); 
      
      /********allocate memory for dust arrays for each dust**********/
      //     for(j=0; j<ncorners[i-1]; j++)
	//		dnumber[i-1]=i+1;
	//		memorydust2_3D(i-1, ncorners[i-1]);
      
      /*print the centre of the dust*/
      // printf("%f %f %d\n", dustcx[i-1], dustcy[i-1], ncorners[i-1]);      
      /*******read coordinates******/
	//		for(j=0; j<ncorners[i-1]; j++)
	//		{	
	//		fscanf(file_dshape, "%lf %lf", &dustx[i-1][j], &dusty[i-1][j], &dustz[i-1][j]);
	//		}
	   }  
	       
      /*normalize initial velocities of the dust*/
      dustvxc[i-1]/=normvel;
      dustvyc[i-1]/=normvel;
      dustvzc[i-1]/=normvel;
	  dustomega[i-1]*=normx/normvel;
      
      /**close the dust input file***/
      fclose(file_dshape);
      
	  /***I CLOSED THE DUST***/
	  
      /*init potential COND*/
    //  if(dtype[i-1]==D_COND)
//	{dphifl[i-1]=Vpr; cond_present=1;}
  //    else
//	{dphifl[i-1]=0.0;} //insulator
      
	printf("File %s is read in genprobe2() and closed\n", arv[i]);
	  
	}
  /*Now we can calculate centre of mass, momentum of inertia, 
    and other static parameters for all the particles*/

for(i=0; i<noofdusts; i++)
{
   if(dshape[i]==0) // spherical dust
     {
	 dmass[i]=4*M_PI*dradius[i]*dradius[i]*dradius[i]*dustrho[i];
     dmass_centr_x[i]=dustcx[i]*dx;
     dmass_centr_y[i]=dustcy[i]*dy;
	 dmass_centr_z[i]=dustcz[i]*dz;
	}
	
   if(dshape[i-1]==0)//polyhedral dust
    {
   //calculate_staticparameters(arc, arv); 
//         for(j=0; j<ncorners[i]; j++)	{
//	  dustxdx[i][j]+=(dustcx[i]-dmass_centr_x[i]/dx);
//	  dustydy[i][j]+=(dustcy[i]-dmass_centr_y[i]/dy);  
//	  dustzdz[i][j]+=(dustcy[i]-dmass_centr_y[i]/d); 
	  //printf("dustx %E dusty %E\n", dustx[i][j], dusty[i][j]);	}
    } 
  /*******find the cetre of dust and reorganize the particle coordinates**/
  /*******new origin is at the centre of mass of the particle*****/
  /*******the particle center is then moved to dustcx,dustcy******/
	
  printf("Dust: %d mass %E mass center: %E %E\n", i, dmass[i], dmass_centr_x[i]/dx, dmass_centr_y[i]/dy);
      
      /******* we can now find a and b coefficients for dust no i**********/
      /******* and find real coordinates***/
  //    findabv(i);
    }
  /*WORK ON THIS PART*/
//  ortnormvec();
//  virtpart();
//   condsquares();
  
  /*******with all the dust particles scanned we can now mark them***/
  /*******on the grid*******/
  /******allocate memory for the look up table******/
  // lut=(int **)malloc(ngx*sizeof(int *));
  for(i=0; i<ngx; i++)
    {
    //  lut[i]=(int *)malloc(4*sizeof(int));
  //    for(j=0; j<4;j++)
//	        lut[i][j]=0;
    } 
//  markgriddust();
  
  //and now create current arrays
  create_currentarrays();

}


/***************************************************************/
void markgriddust(void)
{
  int i,j,k;
  int itemp, jtemp;
  int marker;

  for(i=0; i< ngx; i++)
    { 
      for(j=0; j< ngy; j++)
	{;
	//  if(gp[i][j].marker!=BOUND)     
	//    gp[i][j].marker=NORMAL;  
	}
	for(j=0; j<4 ; j++)	
      lut[i][j]=0;
    }
  /*mark the corners if they are laying on the grid*/
  for(i=0; i<noofdusts; i++)
    {
      for(j=0; j<ncorners[i]; j++)
	if(((dustx[i][j]-(int)dustx[i][j])==0)&&((dusty[i][j]-(int)dusty[i][j])==0))
	  {
	    itemp=(int)dustx[i][j];
	    jtemp=(int)dusty[i][j];
	    /***corner is a grid point***/
	//    gp[itemp][jtemp].marker=dnumber[i];
	    lut[itemp][0]++;
	    if(jtemp<lut[itemp][1] || lut[itemp][1]==0)
	      lut[itemp][1]=jtemp;
	    if(jtemp>lut[itemp][2])
	      lut[itemp][2]=jtemp+1;

	    /***see if there are any vertical points connected to line***/
	    /***if any, mark them***/
	    if(j+1<ncorners[i])
	      {
		if(dustx[i][j]==dustx[i][j+1])
		  {
		    for(k=jtemp; k<(int)dusty[i][j+1]; k++)
		  ;//    gp[itemp][k].marker=dnumber[i];
		    for(k=(int)dusty[i][j+1]; k<jtemp; k++)
		  ;//    gp[itemp][k].marker=dnumber[i];
		  }
	      }
	    else
	      {
		if(dustx[i][j]==dustx[i][0])
		  {
		    for(k=jtemp; k<(int)dusty[i][0]; k++)
		  ;//    gp[itemp][k].marker=dnumber[i];
		    for(k=(int)dusty[i][0]; k<jtemp; k++)
		   ;//  gp[itemp][k].marker=dnumber[i];
		  }
	      }
	  }
    }

  //  for(i=0; i<ngx; i++)
  // printf("before lut %d, %d, %d, %d\n", i, lut[i][0], lut[i][1], lut[i][2] );

  for(i=0; i<ngx; i++)
    checkcolcrossing(i);  /**pre
	pare lut for each column**/
  

  /*******some checks to be deleted later on********/
  for(i=0; i<ngx; i++)
    if(lut[i][0] % 2 != 0 )
      printf("something is wrong i %d -> %d !!!\n", i, lut[i][0]);
  
  //  for(i=0; i<ngx; i++)
     //  printf("lut %d, %d, %d, %d\n", i, lut[i][0], lut[i][1], lut[i][2] );


  /*****print grid!!!**/
  
  
  
   // for(j=0; j<4; j++)
   //   printf("corners %d: %E %E\n", j, dustx[0][j], dusty[0][j]);  
  //  getchar();

  /*****now the lut is created and we can check each grid point***/
  for(i=0; i<ngx; i++)
    if(lut[i][0]!=0)
      {
	marker=0;
       	for(j=lut[i][1]-1; j<=lut[i][2]+1; j++)
	  {
	   // if(gp[i][j].marker==BOUND)
	     // printf("DUST GRAIN REACHED THE BOUNDARY!!!\n");
	    // if(gp[i][j].marker==PROBE)
	    // ; //   marker++;
	   
	    if((marker % 2)==1)
	    ;//  {gp[i][j].marker=lut[i][3];}// printf("marker %d, marking %d %d\n",marker, i,j);} 
	    /***check the next segment***/
	    marker+=checkpointcrossing(i,j);	    
	    //   printf("marker %d %d %d\n", i,j, marker);
	  }
      }

  /*    for(j=ngy-1; j>=0; j--)
    {
      for(i=0; i<ngx; i++)
	if(gp[i][j].marker==PROBE)
	  printf("+ ");
	else
	  printf("%d ", gp[i][j].marker);
      printf("\n");
    }
  */
  //  getchar();

  // for(j=ngy-1; j>=0; j--)
    {
      //   for(i=0; i<ngx; i++)
	//if(gp[i][j].marker==PROBE)
	  //printf("+ ");
	//else
      //	  printf("%d ", gp[i][j].marker);
      //   printf("\n");
    }
// getchar();
 
}
/*******************************************/
void checkcolcrossing(int i)
{
  /*********'first check to create lut table*********/
  int k,l;
  for(l=0; l<noofdusts; l++)
    for(k=0; k<ncorners[l]; k++)
      {	
	if(k<(ncorners[l]-1))
	  {
	    if((dustx[l][k]<i && dustx[l][k+1]>i)||(dustx[l][k]>i && dustx[l][k+1]<i))
	      {
		//printf("cross in %d corner %d of %d\n", i, k, ncorners[l]-1);
		//getchar();
		/***we cross***/
		lut[i][0]++;
		if(dusty[l][k+1]<dusty[l][k])
		  {
		    if((dusty[l][k+1]<lut[i][1]) || (lut[i][1]==0))
		      lut[i][1]=(int)dusty[l][k+1];
		    if(dusty[l][k]>lut[i][2])
		      lut[i][2]=(int)dusty[l][k]+1;
		  }
		else
		  {
		    if((dusty[l][k]<lut[i][1]) || (lut[i][1]==0))
		      lut[i][1]=(int)dusty[l][k];
		    if(dusty[l][k+1]>lut[i][2])
		      lut[i][2]=(int)dusty[l][k+1]+1;
		  }
	      }
	  }
	else
	  {
	    if((dustx[l][k]<i && dustx[l][0]>i) ||(dustx[l][k]>i && dustx[l][0]<i))
	      {
		//printf("cross in %d corner %d of %d\n", i, k, ncorners[l]-1);
		/***we cross***/
		lut[i][0]++;
		if(dusty[l][0]<dusty[l][k])
		  {
		    if((dusty[l][0]<lut[i][1]) || (lut[i][1]==0))
		      lut[i][1]=(int)dusty[l][0];
		    if(dusty[l][k]>lut[i][2])
		      lut[i][2]=(int)dusty[l][k]+1;
		  }
		else
		  {
		    if((dusty[l][k]<lut[i][1]) || (lut[i][1]==0))
		      lut[i][1]=(int)dusty[l][k];
		    if(dusty[l][0]>lut[i][2])
		      lut[i][2]=(int)dusty[l][0]+1;
		  }
	      }
	  }
      }
}

/*******************************************/
int checkpointcrossing(int i, int j)
{
  int k,l;
  int cross=0;
  double yhelp;	
  int kp1, km1, kp2, km2;

  for(l=0; l<noofdusts; l++)
    for(k=0; k<ncorners[l]; k++)
      {	
	/*set indicies*/
	kp1=k+1;
	km1=k-1;
	kp2=k+2;
	km2=k-2;
	
	if(k==0)
	  {
	    km1=ncorners[l]-1;
	    km2=ncorners[l]-2;
	  }
	if(k==ncorners[l]-1)
	  {
	    kp1=0;
	    kp2=1;
	  }

	/*check if crossed the segment*/
	if((dustx[l][k]<i && dustx[l][kp1]>i)||(dustx[l][k]>i && dustx[l][kp1]<i))
	  {	  
	    
	    yhelp=dusta[l][k]*i+dustb[l][k];
	    if(yhelp>=j && yhelp<(j+1)) 
	      { cross++; lut[i][3]=dnumber[l];}
	    if(yhelp==j*1.0)
	      {	  
		  //
		  //gp[i][j].marker=dnumber[l];
	      }
	  }	

	/*****other corners on grid*****/
	//if(k>=1 && k<ncorners[l]-1)
	if((dustx[l][k]==i) && (dustx[l][kp1]!=i))
	  if(dusty[l][k]<j+1 && dusty[l][k]>=j)
	    {
	      if(((dustx[l][kp1]>dustx[l][k])&&(dustx[l][km1]<dustx[l][k])) || ((dustx[l][kp1]<dustx[l][k])&&(dustx[l][km1]>dustx[l][k])))
		{cross++; lut[i][3]=dnumber[l];} //normal corner skip convex corners on grid points
	    }
	
	/******vertical cornrers****/
	
	if(ncorners[l]>=3)
	  {	    
	    if(dusty[l][k]>dusty[l][km1]) //we go up
	      if((dustx[l][k]==i) && (dustx[l][kp1]!=i) && (dustv[l][km1]==1) &&  (dusty[l][k]>=j) && (dusty[l][k]<j+1))
		if((((dustx[l][km2]<dustx[l][k]) && (dustx[l][kp1]>dustx[l][k])) || ((dustx[l][km2]>dustx[l][k]) && (dustx[l][kp1]<dustx[l][k]))))
		 { cross++;lut[i][3]=dnumber[l];}
	    
	    if(dusty[l][kp1]<dusty[l][k] && (dusty[l][k]>=j) && (dusty[l][k]<j+1)) //we go down
	      {
		if((dustx[l][k]==i) && (dustx[l][km1]!=i) && (dustv[l][k]==1) &&  (dusty[l][k]>=j) && (dusty[l][k]<j+1))
		  if((((dustx[l][kp2]<dustx[l][k]) && (dustx[l][km1]>dustx[l][k])) || ((dustx[l][kp2]>dustx[l][k]) && (dustx[l][km1]<dustx[l][k]))))
		   { cross++;lut[i][3]=dnumber[l];}
		
	      }
	  }
	
      }  
  // printf("for i %d, j %d, cross %d\n", i, j, cross % 2 );
  return cross;
}
/********************************************/
/***find a,b,v coefficients and calculate real vertex coord and b*/
void findabv(int i)
{
  int j,jp1;
  for(j=0; j<ncorners[i]; j++)
    {
      jp1=j+1;
      if(j==(ncorners[i]-1))
	jp1=0;
      
      if(dustx[i][j]==dustx[i][jp1])
	{
	  dustv[i][j]=1;
	  dusta[i][j]=dustb[i][j]=0;
	}
      else
	{
	  dustv[i][j]=0;
	  dusta[i][j]=(dusty[i][jp1]-dusty[i][j])/(dustx[i][jp1]-dustx[i][j]);
	  dustb[i][j]=dusty[i][j]-dusta[i][j]*dustx[i][j];
	}
      dustxdx[i][j]=dustx[i][j]*dx;
      dustydy[i][j]=dusty[i][j]*dy;
      dustbdy[i][j]=dustb[i][j]*dy;
    }	

}
 
 /**********memory allocation dust 1***************/
 void memorydust1_3D(int no)
 {
 dpart=(d_particle **)malloc(no*sizeof(d_particle *));
 dpartq=(double **)malloc(no*sizeof(double *));
#ifdef MPI
  rdpart=(d_particle **)malloc(no*sizeof(d_particle *)); 
rdpartq=(double **)malloc(no*sizeof(double *));

#endif
 
  dustx=(double **)malloc(no*sizeof(double *));
  dusty=(double **)malloc(no*sizeof(double *));
  dustz=(double **)malloc(no*sizeof(double *));
  dustxdx=(double **)malloc(no*sizeof(double *));
  dustydy=(double **)malloc(no*sizeof(double *));  
  dustzdz=(double **)malloc(no*sizeof(double *));
  dustq=(double **)malloc(no*sizeof(double *));
  
  dustxdxold=(double **)malloc(no*sizeof(double *));
  dustydyold=(double **)malloc(no*sizeof(double *));
  dustzdzold=(double **)malloc(no*sizeof(double *));
  dustvx=(double **)malloc(no*sizeof(double *));
  dustvy=(double **)malloc(no*sizeof(double *));
  dustvz=(double **)malloc(no*sizeof(double *));
 
  dusta=(double **)malloc(no*sizeof(double *));
  dustb=(double **)malloc(no*sizeof(double *));
  dustbdy=(double **)malloc(no*sizeof(double *));
 
  dustv=(int **)malloc(no*sizeof(int *));
  dustcx=(double *)malloc(no*sizeof(double));
  dustcy=(double *)malloc(no*sizeof(double));
  dustcz=(double *)malloc(no*sizeof(double));
  dustpcx=(double *)malloc(no*sizeof(double));
  dustpcy=(double *)malloc(no*sizeof(double));
  dustpcz=(double *)malloc(no*sizeof(double));


  ncorners=(int *)malloc(no*sizeof(int));
  dtype=(int *)malloc(no*sizeof(int));
  dmove=(int *)malloc(no*sizeof(int));
  dshape=(int *)malloc(no*sizeof(int));
  dradius=(double *)malloc(no*sizeof(double));
  dradiusdx=(double *)malloc(no*sizeof(double));
  dustcxdx=(double *)malloc(no*sizeof(double));
  dustcydx=(double *)malloc(no*sizeof(double));
  dustczdx=(double *)malloc(no*sizeof(double));
 
  dmass=(double *)malloc(no*sizeof(double));
  dmomI=(double *)malloc(no*sizeof(double));
  dmass_centr_x=(double *)malloc(no*sizeof(double));
  dmass_centr_y=(double *)malloc(no*sizeof(double));
  dmass_centr_z=(double *)malloc(no*sizeof(double));
    
  dnumber=(int *)malloc(no*sizeof(int));

   /*for dust movement and initial values*/
  dustaccx=(double *)malloc(no*sizeof(double));
  dustaccy=(double *)malloc(no*sizeof(double));
  dustaccz=(double *)malloc(no*sizeof(double));
   dustvxc=(double *)malloc(no*sizeof(double));
  dustvyc=(double *)malloc(no*sizeof(double));
  dustvzc=(double *)malloc(no*sizeof(double));
   duste=(double *)malloc(no*sizeof(double));
  dustomega=(double *)malloc(no*sizeof(double));
  dphifl=(double *)malloc(no*sizeof(double));
  dustrho=(double *)malloc(no*sizeof(double));
  dustworkfunct=(double *)malloc(no*sizeof(double));
 if(dustrho==NULL) printf("RHO NUOT INUIOO\n");
 if(dphifl==NULL) printf(" NUOT INUIOO\n");
  
  /*for particle crossing Jan version*/
 
   daa=(double **)malloc(no*sizeof(double *));
  dbb=(double **)malloc(no*sizeof(double *));
  dcc=(double **)malloc(no*sizeof(double *));
  daa1y=(double **)malloc(no*sizeof(double *));
  daa1x=(double **)malloc(no*sizeof(double *));
  dbb1y=(double **)malloc(no*sizeof(double *));
  dbb1x=(double **)malloc(no*sizeof(double *));
  dcc1y=(double **)malloc(no*sizeof(double *));
  dcc1x=(double **)malloc(no*sizeof(double *));

  //for virtual particles
  ccorner=(int **)malloc(no*sizeof(int *));
 
  vipcorner=(int **)malloc(no*sizeof(int *));


  //ADDITIONAL
 /*generate "pointing variable" to last particle on dust*/
  dpartlast=(long int *)malloc(no*sizeof(long int));
  dpartmax=(long int *)malloc(no*sizeof(long int));
 
 //P-P (MD part)
  d_localmax=(int *)malloc(no*sizeof(int));
  d_locallist=(int **)malloc(no*sizeof(int *));
	 int i;
	 for(i=0; i< no; i++)
		 d_locallist[i]=(int *)malloc(LIST_SIZE*sizeof(int));
	 
	 //DRAG FORCE
	 drag_direct_x=(double *)malloc(no*sizeof(double));
	 drag_direct_y=(double *)malloc(no*sizeof(double));
	 drag_direct_z=(double *)malloc(no*sizeof(double));
	 drag_elect_x=(double *)malloc(no*sizeof(double));
	 drag_elect_y=(double *)malloc(no*sizeof(double));
	 drag_elect_z=(double *)malloc(no*sizeof(double));
	 drot_z_x1=(double *)malloc(no*sizeof(double));
	 drot_z_x2=(double *)malloc(no*sizeof(double));
	 drot_y_x1=(double *)malloc(no*sizeof(double));
	 drot_y_x2=(double *)malloc(no*sizeof(double));
	 drot_z_y1=(double *)malloc(no*sizeof(double));
	 drot_z_y2=(double *)malloc(no*sizeof(double));
	 drot_y_z1=(double *)malloc(no*sizeof(double));
	 drot_y_z2=(double *)malloc(no*sizeof(double));
	 drot_x_y1=(double *)malloc(no*sizeof(double));
	 drot_x_y2=(double *)malloc(no*sizeof(double));
	 drot_x_z1=(double *)malloc(no*sizeof(double));
	 drot_x_z2=(double *)malloc(no*sizeof(double));
	 
	 //electric
	 elrot_z_x1=(double *)malloc(no*sizeof(double));
	 elrot_z_x2=(double *)malloc(no*sizeof(double));
	 elrot_y_x1=(double *)malloc(no*sizeof(double));
	 elrot_y_x2=(double *)malloc(no*sizeof(double));
	 elrot_z_y1=(double *)malloc(no*sizeof(double));
	 elrot_z_y2=(double *)malloc(no*sizeof(double));
	 elrot_y_z1=(double *)malloc(no*sizeof(double));
	 elrot_y_z2=(double *)malloc(no*sizeof(double));
	 elrot_x_y1=(double *)malloc(no*sizeof(double));
	 elrot_x_y2=(double *)malloc(no*sizeof(double));
	 elrot_x_z1=(double *)malloc(no*sizeof(double));
	 elrot_x_z2=(double *)malloc(no*sizeof(double));	 
	 
 
 }


/**********memory dust allocation 2**********/
 void memorydust2_3D(int j, int nc){
	dustx[j]=(double *)malloc(nc*sizeof(double));
	dusty[j]=(double *)malloc(nc*sizeof(double));
	dustxdx[j]=(double *)malloc(nc*sizeof(double));
	dustydy[j]=(double *)malloc(nc*sizeof(double));
	dustzdz[j]=(double *)malloc(nc*sizeof(double));
	dustq[j]=(double *)malloc(nc*sizeof(double));
	dustxdxold[j]=(double *)malloc(nc*sizeof(double));
	dustydyold[j]=(double *)malloc(nc*sizeof(double));
	dustzdzold[j]=(double *)malloc(nc*sizeof(double));
	dustvx[j]=(double *)malloc(nc*sizeof(double));
	dustvy[j]=(double *)malloc(nc*sizeof(double));
	dustvz[j]=(double *)malloc(nc*sizeof(double));
	
	//FOR DUST COORDINATES, NOT OPTIMIZED!!!
	dpart[j]=(d_particle *)malloc(nc*sizeof(d_particle));	
	dpartq[j]=(double *)malloc(nc*sizeof(double));	
#ifdef MPI
    rdpart[j]=(d_particle *)malloc(nc*sizeof(d_particle));
    rdpartq[j]=(double *)malloc(nc*sizeof(double));	
#endif
	
	dusta[j]=(double *)malloc(nc*sizeof(double));
	dustb[j]=(double *)malloc(nc*sizeof(double));
	dustbdy[j]=(double *)malloc(nc*sizeof(double));
	dustv[j]=(int *)malloc(nc*sizeof(int));
	/*for dust particle hitting: jan version*/
	daa[j]=(double *)malloc(nc*sizeof(double));
	dbb[j]=(double *)malloc(nc*sizeof(double));
	dcc[j]=(double *)malloc(nc*sizeof(double));
	daa1y[j]=(double *)malloc(nc*sizeof(double));
	daa1x[j]=(double *)malloc(nc*sizeof(double));
	dbb1y[j]=(double *)malloc(nc*sizeof(double));
	dbb1x[j]=(double *)malloc(nc*sizeof(double));
	dcc1y[j]=(double *)malloc(nc*sizeof(double));
	dcc1x[j]=(double *)malloc(nc*sizeof(double));
	//for virt particles
	vipcorner[j]=(int *)malloc(nc*sizeof(int));	
	ccorner[j]=(int *)malloc(nc*sizeof(int));
      }
/***********introduce new potential for the probe*******************/
void new_probe_potential(double probepotential)
{
  int i,j;

  for(i=0; i<ngx; i++)
    for(j=0; j<ngy; j++)
//      if(gp[i][j].marker==PROBE)
      	{
//	  gp[i][j].phi=probepotential;
      	}
}

/********Normalize ******************/
void normalize(void)
{
  int i;
  cellvolume=dx*dy*dz; //the real cell volume
  
  /*dimensions*/
  normx=debye[0];
  Lx=Lx/normx;        
  Ly=Ly/normx;
  Lz=Lz/normx;
  dx=dx/normx;
  dy=dy/normx;
  dz=dz/normx;
  // printf("dx dy normalized %f %f \n", dx ,dy); //getchar();
  debyetotal=debyetotal/normx;

  /*time*/
  normtime=1/omegap[0];
  dt=dt/normtime;
  tmax=tmax/normtime;
  
  /*additional factors*/
  dV=dx*dy*dz;
  dxdy=dx*dy;
  dxdydt=dx*dy/dt;
  dVdt=dx*dy*dz/dt;
  dtdx=dt/dx;
  dtdy=dt/dy;
  dtdz=dt/dz;
  dydt=dy/dt;
  dxdt=dx/dt;
  dzdt=dz/dt;
  
  //Velocities, charge2mass, charge dens
  normvel=vthx[0];
  normqm=fabs(qm[0]); //absolute value
  normcharge=fabs(charge[0]); //charge of the real particle
  normmass=mass[0]/ratio; //real mass
  normdens=dens[0]*ratio; //real density
  normqdens=dens[0]*normcharge; //we have dens[0]*ratio*abs[charge]/ratio
	normPP=ratio/(4*M_PI*normdens*normx*normx*normx);
	
  for(i=0; i<S; i++)
    {
      vthx[i]/=normvel;
      vthy[i]/=normvel;
	  vthz[i]/=normvel;
      vdriftx[i]/=normvel;
      qm[i]/=normqm;
      //this is to be used for plasma particles only it is divided by density and already by dV???
		//this is divided by the number of particles in the cell calculated from the theoretical density
      chargeandnorm[i]=(charge[i]/normcharge)*(1/(dens[0]*dV*normx*normx*normx));
//     chargeandnorm[i]=charge[i]/normcharge;	
	   //this is in fact number of elementary charges and  equals the ratio   
      normalcharge[i]=ratio*charge[i]/normcharge;
		
      //now not any more
      // normalcharge[i]=charge[i]/normcharge;
    
	  printf("norm old %E new%E ratio %E dens[0] %E\n", chargeandnorm[i], normalcharge[i], ratio, dens[i]);
     // getchar();
    }

  /*Potential*/
  normpot=fabs(tempx[0]/charge[0]);
  Vpr_begin=Vpr_begin/normpot;
  Vpr_end=Vpr_end/normpot;
  Vpr_step=Vpr_step/normpot;
  Vbound=Vbound/normpot;
  TOLERANCE=TOLERANCE/normpot;
  tolfloating=tolfloating/normpot;


  //CHECK IT ALL!!!
  /*photon flux*/
  if(photons)
      ph_flux*=normx*normx*normtime;

  /*other for diagnostics*/
  normEfield=normpot/debye[0];
  normPE=normpot*normqdens;
  if(rank==0)
    {      
      //      printf("1/normpot %E\n", 1.0/normpot);
      // printf("dxdy %E\n", dxdy);
      //printf("normcharge %E\t %E\n", chargeandnorm[0], chargeandnorm[1]);
  // getchar();
    }
}
/****create current arrays**********/
void create_currentarrays(void)
{
  int i,j;
  current=(long int **)malloc(noofdusts*sizeof(long int *));
  for(i=0; i<noofdusts;i++)
    {
   current[i]=(long int *)malloc(S*sizeof(long int));
   for(j=0; j<S;j++)
       current[i][j]=0;
   }
  curr_av=(long int **)malloc(noofdusts*sizeof(long int *));
  for(i=0; i<noofdusts;i++)
    {
      curr_av[i]=(long int *)malloc(S*sizeof(long int));
      for(j=0; j<S;j++)
	curr_av[i][j]=0;
    }
  rcurr_av=(long int **)malloc(noofdusts*sizeof(long int *));
  for(i=0; i<noofdusts;i++)
    {
      rcurr_av[i]=(long int *)malloc(S*sizeof(long int));
      for(j=0; j<S;j++)
	rcurr_av[i][j]=0;
    }
}


