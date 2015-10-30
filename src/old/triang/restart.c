#include<stdio.h>
#include"const.h"

void dump(long int t)
{
FILE *fdump;
char sdump[45];
int i;
long int help[2];
int help2[2];
double help3[2];
printf("\nI am dumping now at rank%d.\n", rank);

if(rank<10)
   sprintf(sdump, "./data/dump_rank0%d.bin", rank);
else
   sprintf(sdump, "./data/dump_rank%d.bin", rank);


 fdump=fopen(sdump,"w");
//need to dump 
//real time
 help[0]=t;
 fwrite(help, sizeof(long int),1,fdump);
 //particles
 fwrite(npart, sizeof(long int),S,fdump);

 for(i=0; i< S; i++)
{
  //printf("dumped number of particles %d is %ld\n", i, npart[i]);
 fwrite(spec[i].part, sizeof(particle),npart[i],fdump);
//remaining flux
 fwrite(extrapart[i], sizeof(double), 6,fdump); 
}
//dust grains
 help2[0]=noofdusts;
 fwrite(help2, sizeof(int), 1, fdump);
 fwrite(dshape, sizeof(int), noofdusts, fdump);
 //centres	
 fwrite(dmass_centr_x, sizeof(double), noofdusts,fdump);
 fwrite(dmass_centr_y, sizeof(double), noofdusts,fdump);
 fwrite(dmass_centr_z, sizeof(double), noofdusts,fdump);
 fwrite(dustvxc, sizeof(double), noofdusts,fdump);
 fwrite(dustvyc, sizeof(double), noofdusts,fdump);
 fwrite(dustvzc, sizeof(double), noofdusts,fdump);
 fwrite(dustomega, sizeof(double), noofdusts,fdump);
 fwrite(dmass, sizeof(double), noofdusts,fdump);
 fwrite(dmomI, sizeof(double), noofdusts,fdump);
 fwrite(ncorners, sizeof(int), noofdusts,fdump);
 fwrite(dpartlast, sizeof(long int), noofdusts,fdump);
 fwrite(dpartmax, sizeof(long int), noofdusts,fdump);
 //dump dust workfunction
 fwrite(dustworkfunct, sizeof(double), noofdusts, fdump);
 printf("dustworkfunct %E\n", dustworkfunct[0]);
 for(i=0; i<noofdusts; i++)
   {
	 if(dshape[i]==0)
		{
		help3[0]=dradiusdx[i];
		fwrite(help3, sizeof(double), 1,fdump);
		help3[0]=dustcxdx[i];
		fwrite(help3, sizeof(double), 1,fdump);
		help3[0]=dustcydx[i];
		fwrite(help3, sizeof(double), 1,fdump);
		help3[0]=dustczdx[i];
		fwrite(help3, sizeof(double), 1,fdump);
		}
	 if(dshape[i]==1)
	 {
      fwrite(dustx[i], sizeof(double), ncorners[i],fdump);
      fwrite(dusty[i], sizeof(double), ncorners[i],fdump);
      fwrite(dustz[i], sizeof(double), ncorners[i],fdump);
	  fwrite(vipcorner[i], sizeof(int), ncorners[i], fdump);
      fwrite(ccorner[i], sizeof(int), ncorners[i], fdump);
	 }
	   
     fwrite(dpart[i], sizeof(d_particle), dpartlast[i],fdump);
     fwrite(current[i], sizeof(long int), S, fdump); 
	 fwrite(curr_av[i], sizeof(long int), S, fdump); 
     int j;
     j=2;
     // printf("radius %E\ and center %E %E %E\n", dradiusdx[i], dustcxdx[i], dustcydx[i], dustczdx[i]);
	   
	   //    printf("corner %E\t%E dpartlast %E %E\n", dustx[i][j], dusty[i][j], dpart[i][dpartlast[i]-1].x, dpart[i][dpartlast[i]-1].y);
     
   }
 //photons
 if(photons)
   {
     for(i=0; i<noofdusts; i++)
       {
	 fwrite(dustxnormv[i], sizeof(double), ncorners[i], fdump);
	 fwrite(dustynormv[i], sizeof(double), ncorners[i], fdump); 
	// fwrite(dustznormv[i], sizeof(double), ncorners[i], fdump);
       }
   }
 
 fclose(fdump);
 
 printf("Finished dumping data at rank%d for t=%E s.\n",rank,t*normtime);
 if(rank==0)
   printf("To start the simulation again from this point, folow instructions in the manual.\n");
 
 fdump=fopen(sdump,"r");
 //need to dump 
 //real time
 //help[1]=t;
 fread(help, sizeof(long int),1,fdump);
 //particles
 fread(npart, sizeof(long int),S,fdump);
 printf("HALLO help e here noofdusts at t:%d is %d\n", help[0],  npart[1]);
 i=1;
 printf("b Lx %E part i %E %E %E %E\n", Lx,spec[i].part[100].x, spec[i].part[100].y, spec[i].part[100].vx, spec[i].part[100].vy);
 for(i=0; i< S; i++)
   {
     fread(spec[i].part, sizeof(particle),npart[i],fdump);
     //remaining flux
     fread(extrapart[i], sizeof(double), 6,fdump); 
   }
 i=1;
 printf("Lx %E part i %E %E %E %E\n", Lx,spec[i].part[100].x, spec[i].part[100].y, spec[i].part[100].vx, spec[i].part[100].vy);
 
 fclose(fdump);
}

long int prog_restart()
{
  FILE *fdump;
  char sdump[45];
  int i;
  long int help[2];
  int help2[2];
   double help3[2];
  //printf("\nI am reading now at rank%d.\n", rank);

  if(rank<10)
    sprintf(sdump, "./data/dump_rank0%d.bin", rank);
  else
    sprintf(sdump, "./data/dump_rank%d.bin", rank);
  
  
  if((fdump=fopen(sdump,"r"))==NULL) 
    printf("cannot open file %s \n",sdump);
  printf("File %s open\n",sdump);
  
  fread(help, sizeof(long int),1,fdump); //this is time
  
  //particles
  fread(npart, sizeof(long int),S,fdump);
  printf("read no of part 0 %d 1 %d",npart[0], npart[1]);  
  for(i=0; i< S; i++)
    {
      fread(spec[i].part, sizeof(particle),npart[i],fdump);
      //remaining flux
      printf("Lx %E part i %E %E %E %E\n", Lx,spec[i].part[100].x, spec[i].part[100].y, spec[i].part[100].vx, spec[i].part[100].vy);
      fread(extrapart[i], sizeof(double), 6,fdump); 
    }
  //dust grains
  
  fread(help2, sizeof(int), 1, fdump);
  noofdusts=help2[0];
 //assign memory
 memorydust1_3D(noofdusts);
 
 memoryduststatic(noofdusts);
  printf("HALLO help e here noofdusts at t:%d is %d npart %d %d %d\n", help[0], help2[0], npart[0],npart[1], noofdusts);
 
 fread(dshape, sizeof(int), noofdusts, fdump);	
 //centres
 fread(dmass_centr_x, sizeof(double), noofdusts,fdump);
 fread(dmass_centr_y, sizeof(double), noofdusts,fdump);
 fread(dmass_centr_z, sizeof(double), noofdusts,fdump);
 fread(dustvxc, sizeof(double), noofdusts,fdump);
 fread(dustvyc, sizeof(double), noofdusts,fdump);
 fread(dustvzc, sizeof(double), noofdusts,fdump);
 fread(dustomega, sizeof(double), noofdusts,fdump);
 fread(dmass, sizeof(double), noofdusts,fdump);
 fread(dmomI, sizeof(double), noofdusts,fdump);
 fread(ncorners, sizeof(int), noofdusts,fdump);
 fread(dpartlast, sizeof(long int), noofdusts,fdump);
 fread(dpartmax, sizeof(long int), noofdusts,fdump);
 //dust type and workfunction 
 fread(dustworkfunct, sizeof(double), noofdusts, fdump);
 printf("dustworkfunction %E dpart %d ions %d, ncorners %d, xc %E\n", dustworkfunct[0], dpartlast[0], dpartlast[0], ncorners[0], dmass_centr_x[0]);
//  printf("corner %E\t%E dpartlast %E %E\n", dustx[i][j], dusty[i][j], dpart[i][dpartlast[i]-1].x, dpart[i][dpartlast[i]-1].y);
//memory for particles
/*I initialize dust particles here*/
  dpart=(d_particle **)malloc(noofdusts*sizeof(d_particle *));
  for(i=0; i<noofdusts; i++)
    { 
      dpart[i]=(d_particle *)malloc(dpartmax[i]*sizeof(d_particle));
  
    }
  if (dpart==NULL)
    { 
      printf("Error: Out of Memory rank:%d\n", rank);
      exit(1);
    }
	/*those I do not need**/
// d_polygon(arc, arv);  /***check it here***/ //to here ok
// d_centreofmass_and_momI(); 


  create_currentarrays();
  int j;
  for(i=0; i<noofdusts; i++)
    {
    
		if(dshape[i]==0)
		{
			memorydust2_3D(i,dpartlast[i]);
		fread(help3, sizeof(double), 1,fdump);
			printf("check here %E\n", help3[0]);
		dradiusdx[i]=help3[0];			
		fread(help3, sizeof(double), 1,fdump);
		dustcxdx[i]=help3[0];
		fread(help3, sizeof(double), 1,fdump);
		dustcydx[i]=help3[0];
		fread(help3, sizeof(double), 1,fdump);
		dustczdx[i]=help3[0];
		}
		if(dshape[i]==1)	
		{
			memorydust2_3D(i,ncorners[i]);
			printf("ncorners %d\n", ncorners[i]);
      fread(dustx[i], sizeof(double), ncorners[i],fdump);
      fread(dusty[i], sizeof(double), ncorners[i],fdump);
      fread(dustz[i], sizeof(double), ncorners[i],fdump);
	  fread(vipcorner[i], sizeof(int), ncorners[i], fdump);
      fread(ccorner[i], sizeof(int), ncorners[i], fdump);  
		}
	  fread(dpart[i], sizeof(d_particle), dpartlast[i],fdump);
      fread(current[i], sizeof(long int), S, fdump);
	  fread(curr_av[i], sizeof(long int), S, fdump);
      //for(j=0; j<ncorners[i]; j++)
      j=2;
		printf("radius %E\ and center %E %E %E\n", dradiusdx[i], dustcxdx[i], dustcydx[i], dustczdx[i]);
		
	 printf("corner %E\t%E dpartlast %E %E\n", dustx[i][j], dusty[i][j], dpart[i][dpartlast[i]-1].x, dpart[i][dpartlast[i]-1].y);
    //  findabv(i); 
    }
 
  if(photons)
    {
      printf("photons!!\n"); getchar();
      dustxnormv=(double **)malloc(noofdusts*sizeof(double *));
      dustynormv=(double **)malloc(noofdusts*sizeof(double *));
	 // dustznormv=(double **)malloc(noofdusts*sizeof(double *));

	for(i=0; i<noofdusts; i++)   
	{
	  dustxnormv[i]=(double *)malloc(ncorners[i]*sizeof(double));
	  dustynormv[i]=(double *)malloc(ncorners[i]*sizeof(double));	
	  //dustznormv[i]=(double *)malloc(ncorners[i]*sizeof(double));
	  fread(dustxnormv[i], sizeof(double), ncorners[i], fdump);
	  fread(dustynormv[i], sizeof(double), ncorners[i], fdump); 
	 // fread(dustznormv[i], sizeof(double), ncorners[i], fdump); 
	}            
    } 
  printf("closing the file at rank %d, %s\n", rank, sdump);
  fclose(fdump);
  printf("closed the file at rank %d\n", rank);
#ifdef MPI 
  printf("waiting rank %d\n\n", rank);
 MPI_Barrier(MPI_COMM_WORLD);

#endif

/*******with all the dust particles scanned we can now mark them***/
/******allocate memory for the look up table******/
  lut=(int **)malloc(ngx*sizeof(int *));
  for(i=0; i<ngx; i++)
    {
      lut[i]=(int *)malloc(3*sizeof(int));
      lut[i][0]=lut[i][1]=lut[i][2]=0;
    } 
  //markgriddust();
  
  shift_while_restarting(1, 1.0, 1.0 ,0);
  printf("Finished reading data at rank%d for t=%E s.\n",rank,help[0]*normtime);
  if(rank==0)
    printf("Restarted successfully.\n");
  // getchar();

  return help[0]+1;	//returning time
}


void shift_while_restarting(int dno, double x, double y, double z)
{
  double xdx, ydx, zdx;
  int i,j,ii;

  xdx=x*dx;
  ydx=y*dx;
  zdx=z*dx;

  printf("XDX %E %E %E\n", xdx, ydx, zdx);
  dmass_centr_x[dno]+=xdx;
  dmass_centr_y[dno]+=ydx;
  dmass_centr_z[dno]+=zdx;

  if(dshape[dno]==0)
    {
  dustcxdx[dno]+=xdx;
  dustcydx[dno]+=ydx;
  dustczdx[dno]+=zdx;
  //  dradiusdx[dno]*=1;
    }
  for(i=0; i<dpartlast[dno]; i++)
    {
      dpart[dno][i].x+=xdx;
      dpart[dno][i].y+=ydx;
      dpart[dno][i].z+=zdx;
    }

  if(rank>0)
    {
      for(i=0; i<S; i++)
	lostpart[i]=0;
      //delete particles within the dust
      for(i=0; i<S; i++)
	{
	  printf("checking particles rank %d npart %d\n", rank, npart[i]);
	  for(ii=0; ii<npart[i]; ii++)
	    {
	      //      printf("checking particles rank %E %E %E, dust cx %E %E %E radius %E\n",  spec[i].part[ii].x, spec[i].part[ii].y, spec[i].part[ii].z, dustcxdx[dno], dustcydx[dno],dustczdx[dno], dradiusdx[dno]);
	      if(initpartcheck_restart(dno, spec[i].part[ii].x, spec[i].part[ii].y, spec[i].part[ii].z, 0.0001*dx) % 2 != 0)
		{
		  lostpart[i]++;  
		  lostlist[i][lostpart[i]-1]=ii;	
		 		  printf("particle inside \n");
		}
	    }
	  
	  // printf("end\n");
	  /*shift particles*/
	  int allpart;
	  allpart=npart[i];
	  
	  
	  printf("rank %d specie %d allpart %d, lostpart %ld\n", rank, i, allpart, lostpart[i]);
	  // fprintf(history,"TIME %E, specie %d, allpart %d, lostpart %d\n",t*dt*normtime, i,  allpart, lostpart[i]);	  
	  
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
	  printf("rank %d specie %d allpart %d, lostpart %ld\n", rank, i, allpart, lostpart[i]);
	}
    }
} 

