/* DiP3D */
/* Main function */
/* Author: Wojciech Jacek Miloch */
/* University of Oslo, Norway */
/* 2009 */

#include "const.h"

int main(int argc, char *argv[])
{
/*****************INITIALIZATIONS********************/
/*checking simulation time*/

clockstart=clock();
(void) time(&timestart);


// t is the time step 
long int t,t_init=0,timeend;

int printing=1;
int iterationmax, iteration;
int ii; //dust numbers
double condcheck;
double cup,cdown,cstep;

if(rank==0)  
   history=my_file_open("./data/history.dat", "w");


  /*Initialize MPI*/
  numtasks=1;
  rank=0;
#ifdef MPI 
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks); //no of processesors
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); //numbers of processes
#endif

 if(rank==0) 
    printf("DiP3D, Dust in Plasma 3D simulation\n Author: Wojciech J. Miloch\n");

#ifdef EPRO
  if(rank==0)
    printf("Program is in the mode for a biased probe\n\n");

#ifdef CHAR
  if(rank==0)
    {
      printf("Program will perform the probe characteristic only\n");
      printf("!!!No other diagnostics will be plotted!!!\n\n\n");
    }
#endif 
#endif
  if(rank==0)
    printf("Program is in the mode for an insulator/conductor type dust particles\n\n");
	  
#ifdef BEAM
  if(rank==0)
    printf("!!!There is an additional ion beam!!!\n !!!Check if you use the proper input file!!!\n");
#endif  

/*************************INPUT ********************/
 if(rank==0)
    convert(); //input.c - convert input file
 
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD); //wait until the file is converted
#endif
  readdata(argc, argv); //input.c - read input param,setvariables
  printf("Rank %d, Initialization: I read the input data, now proceed with the flux calculation\n", rank);
  calculate_flux(); //flux.c - calcul. flux via boundaries on nodes 1,2,3..
 
/*************************GENERATE ********************/
  memorygrid(); //grid.c - allocate memory for the grid
  normalize(); //grid.c - normalize and find norm factors

//  if(photons) //NO PHOTONS YET , simple mode
 //   photonflux();
  gen_boundaries(); //grid.c - generate boundaries

 init_primeroot(rank/(numtasks*1.0));
#ifndef RESTART  
//IMPORTANT!!!!!
//not used at the moment -> we are in the test phase
  gen_dust3D(argc, argv);
  printf("rank %d finished generating dust \n", rank);

#else
  t_init=prog_restart(); //change time for starting and initialize!
  printf("finisherd restarting at rank %d", rank);
#endif

#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  init_newpart(); //flux.c -initialise param. for particle injection
  diagn_open(); //diagn.c - open files for diagnostics 
  pot_probes_init();
	
  if(rank==0)
     mglin_init(ngx,ngy,ngz); //initialize field solver
  	
  //for each probe potential do the following
#ifdef EPRO
  for(Vpr=Vpr_begin; Vpr<=Vpr_end; Vpr+=Vpr_step)
    {
	 cond_present=0; //we have probe, no iterative procedure for finding the floating potential
	  for(ii=0; ii<noofdusts; ii++)
	   dphifl[ii]=Vpr; //put potential on dusts/probes
#endif

//CHECK IF THERE ARE CONDUCTING DUSTS
    iterationmax=1;
 
    if(cond_present)
      {
	printf("CONDUCTING OBJECTS ARE PRESENT BUT THE PROGRAM IS NOT READY FOR THIS CHALLANGE!\n");
	exit(1);
	iterationmax=0;
	condcheck=fabs(Vpr_end-Vpr_begin);
	condcheck=condcheck/2;
	printf("initial %f\n", condcheck);
	while(condcheck>tolfloating) 
	  {
	    condcheck=condcheck/2;
	    iterationmax++;
	    if(rank==0)
	      printf("iter %d, Vpr (%E,%E)  condcheck %f, tolfloating %f\n", iterationmax, Vpr_end*normpot, Vpr_begin*normpot, condcheck*normpot, tolfloating*normpot);
	  }
	if(rank==0)
	  printf("I will do %d iterations to find the floating potential on the object\n", iterationmax);
	/*I calculate new potential, used only on rank 0*/
	cup=Vpr_end;
	cdown=Vpr_begin;
	cstep=(cup-cdown)/2;
	for(ii=0; ii<noofdusts; ii++)
	  dphifl[ii]=cdown+cstep;
	//  Vpr=cdown+cstep;
      }
    
    for(iteration=0; iteration<iterationmax; iteration++)
      {
	printf("I will make %d whole cycles\n", iterationmax);
	printing=0; 
	if(iteration==iterationmax-1)
	 {
	   printing=1;
	   printf("Diagnostics will be printed in this cycle\n");
	 }
	      
	if((numtasks>1 && rank!=0)||(numtasks==1 && rank==0))
#ifndef RESTART
//3D starts from here
	

	  gen_bgnd(); // generate.c - generate background on nodes > 0 if there are, else on rank0
//it is done soso Box Muller twice...
#else
	;	
	//GENERATE BACKGROUND FOR EPRO AND COND, FURTHER ITERATONS AND T_INIT=0
	//IF varialbe iterationmax is greater than 1
#endif
	cleargrid(); //grid.c -cleargrid,new poten. on probe
    weighting1(); //grid.c - linear weightning
		  weightingdust1(1);
		  
		  //		  	printf("parameters irho ngrid %E\n", irho[6][1][1][1]);	
#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
	printf("rank %d reducing density\n", rank);
	MPI_Reduce(rho,rrho,rhoMAX,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
	if(rank==0)
	  {
#ifdef MPI
	    mglin(rrho, NCYCLES); //OK
#else
	    mglin(rho,NCYCLES); //OK
#endif
	    electric_field(); // gauss.c  		
	  }
	if(rank==0)
	  timeend=(int)(tmax/dt);
	//MPI BROADCAST ELECTRIC FIELD DATA HERE
#ifdef MPI
	MPI_Bcast(Fs,FsMAX,MPI_DOUBLE,0,MPI_COMM_WORLD);  
	MPI_Bcast(&timeend,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
#endif
#ifndef RESTART  
      if((numtasks>1 && rank!=0)||(numtasks==1 && rank==0))
	   accel(-0.5); //accel2.c - find vel for t=-0.5dt" 
//	FOR EPRO AND COND, FURTHER ITERATONS AND T_INIT=0 
#endif
	 
		  
	/***********  MAIN LOOP *****************/
   
      for(t=t_init; t<=timeend; t++) /*MAIN LOOP*/
	{
	 if(rank==0)
	   {
	     printf("Rank %d Vpr %6.3E\tTIME %6.3E \t (iteration %ld). LEFT %ld iterations.\n\n\n", rank, Vpr*normpot, t*dt*normtime, t, timeend-t);
	     pot_probes(t);
		 printf("Rank %d, after pot_probes\n", rank);
	   }
		
	 /************DIAGNOSTICS****************/
#ifndef CHAR
	 if(printing==1)
         printgrid(t); //diagn.c - print grid quantities
				 
#endif

// printf("here\n");
	  if(rank==0)
	    {	  	
		printdustchargetime(dust_time, t, weight);
		printf("Rank %d after printdustchargetime\n", rank);
			  //    printdustshapetime(t);
	  
	    //  printdth(dhist,t);
     //printf("IN MAIN.c printdth\n");
	    }

	  printf("entering dmove rank %d\n", rank);
	 /**********MOVE -> TRAJECTORIES************/	     
	  // d_move(t);  //DUST MUST BE MOVED ON ALL ??? hmm... different nodes make it more difficult
	//	printf("dii %d %E\n", 41, dpart[0][41].q);
	//	getchar();

	 //printf("in markgriddust \n");
	// markgriddust();	   
	 printf("entering move rank %d\n", rank);	
	 //printf("\n\n\nMAIN: in diagnostics 2\n");	
		
		//PP part and MPI
	
		int dno;
#ifdef MPI
        printf("rank %d arrived in accel\n", rank);	
		MPI_Barrier(MPI_COMM_WORLD);
	//	printf("check1 before move rank %d dpart %E\n", rank, dpart[0][4].q);
		//dpart[0][4].q=3.0;
		for(dno=0; dno<noofdusts; dno++)
			for(ii=0; ii<dpartlast[dno]; ii++)
			{
				dpartq[dno][ii]=0.0;
				//printf("check1 before move rank %d dpart %E\n", rank, dpart[0][4].q);
				MPI_Reduce(&dpart[dno][ii].q,&dpartq[dno][ii],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			}
		
		MPI_Barrier(MPI_COMM_WORLD);
		for(dno=0; dno<noofdusts; dno++)
			for(ii=0; ii<dpartlast[dno]; ii++)
				MPI_Bcast(&dpartq[dno][ii],1,MPI_DOUBLE,0,MPI_COMM_WORLD); 
	//	printf("check before move rank %d dpartq %E\n", rank, dpartq[0][4]);
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		
#ifndef MPI	
		//assign the charge to the dust
		for(dno=0; dno<noofdusts; dno++)
			for(ii=0; ii<dpartlast[dno]; ii++)
				dpartq[dno][ii]=dpart[dno][ii].q;
#endif	
		
	 if((numtasks>1 && rank!=0)||(numtasks==1 && rank==0))
	   {
	     move(t); //accel3.c accel and move part.  
	     printf("if photons then start at t %d\n",(int)(0.0/(dt*normtime*omegap[1]/(2*M_PI))));
	      if(photons)
		{
		  // if(t*dt*normtime*omegap[1]/(2*M_PI)> 0.0)
		  if(t>2500 && t<5000)
		    {
		      if((numtasks>1 && rank!=0)||(numtasks==1 && rank==0))
			//only on nodes with particles plasma
			photoelectriceffect();
		      //getchar();
		    }
		  //printf("MAIN: finished photoelectric\n");
		}
	   }
		printf("before electric\n");		
		drag_force_electric();
		printf("after electric\n");
		printdragforce(t);
		printf("after printing drag\n");
		
//	 printf("MAIN: main check\n");	
#ifdef MPI 
	 MPI_Barrier(MPI_COMM_WORLD);
	 MPI_Reduce(KE,rKE,KEMAX, MPI_DOUBLE,MPI_SUM,1, MPI_COMM_WORLD);
#endif
     /***************DIAGNOSTICS 2 *********************/
#ifndef CHAR

	if((printing==1) && ((t % average) == 0))
	  printKEall(t);
#endif
//printf("Current el %d ion %d\n", current[0][0], current[0][1]);
	  //average current from last 1000 steps
	  //SUM ALL
		print_current(t);
		
		if(iteration<iterationmax-1)
	    {if(t>=(timeend-1000))
	      average_current(); 
	    }
#ifndef EPRO		
	  else
	    {
	      average_current();
	      if(((t % 1000) ==0) && (t>0))
		//if there are no conductors this is an empty statement
		;//	findnewpotentials(cstep,1000,curr);
	      //adjust floating potential
	    }
#endif		    	   
	  printf("MAIN: new particles & weighting\n");
	 /*************NEW PARTICLES & WEIGHTING******************/  
	//  init_primeroot(0);
	  if((numtasks>1 && rank!=0)||(numtasks==1 && rank==0))	  
	    newparticles(t); //generate.c - inject new part. 3D	  	  
	  
	  cleargrid(); //grid.c
	  weighting1(); //grid.c - linear weightning
	  weightingdust1(1);
	  
	  if(t==7000)
	    {
	      dump(t);
	      
		 //exit(1);
	    }
	  if(t==10)
	    ;//	    exit(1);
	  
	 /****************FIELDS**************/
#ifdef MPI
	 MPI_Barrier(MPI_COMM_WORLD);
	 printf("reducing density\n");
	 MPI_Reduce(rho,rrho,rhoMAX,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
        if(rank==0)
	  {
	    //   gauss_seidel(ngx,ngy,TOLERANCE); //gauss.c -solve poisson eq.
#ifdef MPI
	    mglin(rrho,NCYCLES);
#else
	    mglin(rho,NCYCLES); //take also rrho	     	      
#endif
	    printf("MAIN: potential solved on rank %d\n", rank);
	    electric_field(); //gauss.c 
	  }
	//MPI BROADCAST DATA
#ifdef MPI
	MPI_Bcast(Fs,FsMAX,MPI_DOUBLE,0,MPI_COMM_WORLD);  
	MPI_Barrier(MPI_COMM_WORLD);      
#endif	
	
	//	if(t==1000)
	//{	/*
	// 	 for(j=ngy-1; j>=0; j--)
        // {
	//for(i=0; i<ngx; i++)
	//  printf("%3.1E ", gp[i][j].phi);
	//printf("\n");
	//}*/
	    // i=0;
	   //	printf("b Lx %E part i %E %E %E %E\n", Lx,spec[i].part[100].x, spec[i].part[100].y, spec[i].part[100].vx, spec[i].part[100].vy);
	   // i=1;
       //printf("b Lx %E part i %E %E %E %E\n", Lx,spec[i].part[100].x, spec[i].part[100].y, spec[i].part[100].vx, spec[i].part[100].vy);

    //exit(1);
//}
//getchar();	

	}     
	/**********************END MAIN CYCLE*********************/

//FIND NEW POTENTIALS ON DUST GRAINS and print current!!!
cstep=cstep/2.0;
//findnewpotentials(cstep,1000,curr);
}

#ifdef EPRO
    }
#endif


/***********end main cycle for conductors*********/

/************************FINALIZE***********************/
if(rank==0)
  mglin_destroy();
  
 diagn_close(); //close diagnostics
 memorygridfree(); //free memory

int ti;
ti=0;
printf("rank %d FLUX for %d:\n c0 %E \n c1 %E c2 %E c3 %E c4 %E c5 %E\n", rank, ti, 1.0*c0[ti], 1.0*c1[ti], 1.0*c2[ti], 1.0*c3[ti], 1.0*c4[ti], 1.0*c5[ti]);
ti=1;
printf("rank %d FLUX for %d:\n c0 %E \n c1 %E c2 %E c3 %E c4 %E c5 %E\n",  rank, ti, 1.0*c0[ti], 1.0*c1[ti], 1.0*c2[ti], 1.0*c3[ti], 1.0*c4[ti], 1.0*c5[ti]);
#ifdef MPI
  MPI_Finalize();
#endif
  clockend=clock();
  (void) time(&timeending);
  timeelapsed=((double)(clockend-clockstart));
  if(rank==0)
    {
      fprintf(history, "CPU time: %E\t / %ld = %E\n", timeelapsed, CLOCKS_PER_SEC, timeelapsed/CLOCKS_PER_SEC);
      fprintf(history, "Total time %d seconds\n", (int)(timeending-timestart));
      fclose(history);
    }
  return 0;
}


