/*DiP3D */
/*Author: Wojciech Jacek Miloch */
/*University of Oslo, Norway */
/*2009*/
/*last revised 05.03.2009*/
#include <math.h>
#include "const.h"

/*function converting the intput.txt file into program readable file*/  
void convert(void)
{
  FILE *open, *write;
  char l;
  open=my_file_open("input.txt", "r");
  write=my_file_open("./src/newfile.dat", "w");  
  do 
    {
      l=getc(open);
      if(((l>='0') && (l <='9')) || l=='\n' || l=='E' || l=='.' || l=='-' || l==' ')
	putc(l,write);
    }
  while(l!=EOF);  
  fclose(open);
  fclose(write);
}

/*read data and print the initial values in print.txt file*/
void readdata(int arc, char *arv[]) 
{
  FILE *open,*write;
  int i;
  double Npel, Npio, Cs; /*plasma number, speed sound*/
  
#ifdef BEAM
double Npbe;
#endif  
  double partpercell, dustvolume; 
  double vxmax,vxmax1,vymax,vymax1,vzmax,vzmax1;

  
  open=my_file_open("./src/newfile.dat", "r");     
  if(rank==0)
    write=my_file_open("./data/print.txt", "w");
	
  fscanf(open,"%lf", &Lx);
  if(rank==0)
    fprintf(write,"System dimension, Lx:\t %f m\n", Lx);
  fscanf(open,"%lf", &Ly);
  if(rank==0)
    fprintf(write,"System dimension, Ly:\t %f m\n", Ly);
  fscanf(open,"%lf", &Lz);
  if(rank==0)
	fprintf(write,"System dimension, Lz:\t %f m\n", Lz);	
  if((Lx>Lx_MAX) || (Ly>Ly_MAX) || (Lz>Lz_MAX))
    {
      if(rank==0)
		{
		printf("ERROR, chamber dimensions bigger than maximum values: Lx_MAX: %d, Ly_MAX: %d, Lz_MAX: %d. Simulation aborted\n", Lx_MAX, Ly_MAX, Lz_MAX);
		fclose(open); fclose(write);
		}  
      exit(-1);
    }
  fscanf(open, "%d", &ngx);
  fscanf(open, "%d", &ngy);
  fscanf(open, "%d", &ngz);
  fscanf(open, "%lf", &dt); //time step 
  fscanf(open, "%lf", &tmax);
  
  if(rank==0)
    fprintf(write,"Number of grid points:\n x direction: %d\n y direction: %d\n z direction: %d\n", ngx,ngy,ngz);

  if(ngx>ngx_MAX)
    {
      if(rank==0)
	{
	  printf("ERROR, too many grid cells in x direction ngx=%d>ngx_MAX=%d\n", ngx, ngx_MAX);
	  fclose(open); fclose(write);
	}  
      exit(-1);
    } 

  if(ngy>ngy_MAX)
    {
      if(rank==0)
	{
	  printf("ERROR, too many grid cells in y direction ngy=%d>ngy_MAX=%d\n", ngy, ngy_MAX);
	  fclose(open); fclose(write); 
	} 
      exit(-1);
    }  

  if(ngz>ngz_MAX)
    {
      if(rank==0)
	{
	  printf("ERROR, too many grid cells in y direction ngy=%d>ngy_MAX=%d\n", ngz, ngz_MAX);
	  fclose(open); fclose(write); 
	} 
      exit(-1);
    }  

  dx=Lx*1.0/(ngx-1); //subtract 1 because ngx is number of gridpoints start from 0
  dy=Ly*1.0/(ngy-1);
  dz=Lz*1.0/(ngz-1);

  if(rank==0)
    {
      fprintf(write,"Grid spacing dx:%f m\nGrid spacing dy:%f m\nGrid spacing dz:%f\n", dx, dy, dz);
      fprintf(write,"Time step dt:%E  s \n", dt);
      fprintf(write,"Time max:    %E  s \n", tmax);
    }

  /*reading density of realpart particles*/
  fscanf(open,"%lf", &dens[0]);
  /*reading number of simulation particles*/
  fscanf(open,"%ld", &npart[0]);
  /*reading electron parameters*/  
  fscanf(open,"%lf", &qm[0]); 
  fscanf(open,"%lf", &debye[0]);  
  fscanf(open,"%lf", &vdriftx[0]);
  /*reading ion parameters*/
  fscanf(open, "%lf", &charge[1]); //no of elem charges per real ion
  fscanf(open, "%lf", &qm[1]);
  fscanf(open,"%lf", &vdriftx[1]);  
  fscanf(open, "%lf", &ti2te);
  ti2te=1.0/ti2te;
#ifdef BEAM
 /*reading ion beam*/
  fscanf(open, "%lf", &vdriftx[2]);
  fscanf(open, "%ld", &npart[2]);
#endif
  /*reading other parameters*/
  // fscanf(open,"%d", &dustmove); //I do not use it any more, data in the dust input file
  fscanf(open,"%lf %lf %lf", &Vpr_begin, &Vpr_end, &Vpr_step);
  fscanf(open,"%lf", &Vbound);
  fscanf(open,"%lf", &TOLERANCE);
  fscanf(open,"%lf", &tolfloating);
  /*reading photons*/
  fscanf(open,"%lf", &ph_flux);
  fscanf(open,"%lf", &ph_angle);
  fscanf(open,"%lf", &ph_energy);
  if(ph_flux != 0) 
    photons=1; 
  else 
    photons=0;
  /*reading diagnosics*/
  fscanf(open,"%d", &diagint);
  //allocate memory	
  diagint_av=ivecmem(0, diagint-1);
  diagint_st=ivecmem(0, diagint-1);
  for(i=0; i< diagint; i++)
	{     
      fscanf(open,"%d %d", &diagint_st[i], &diagint_av[i]);
      //temporary check in input.c
      printf("%d\n",diagint_st[i]);
	}
	
  /*set some parameters*/
  allpart=0;

  //DUST VOLUME   
	dustvolume=finddustvolume(arc,arv)*dx*dy*dz;
//	printf("dust volume: %E\n", dustvolume);
	//*dx*dy*dz;
	
   ratio=(dens[0]*(Lx*Ly*Lz-dustvolume))/npart[0]; //ratio real to sim, no particles introduced on the objects
   /*initiate npart on different nodes*/
   if(numtasks>1)
    {
      npartinit[0]=(int)((npart[0]*1.0)/(numtasks-1));
      if((npart[0] % (numtasks-1))!= 0)
	if(rank==0) printf("Warning: no of electrons != n*(numtasks-1), you have now %ld electrons\n", npartinit[0]*(numtasks-1));
    }
  else
    npartinit[0]=npart[0];

  /*calculating & printing electron parameters*/
   charge[0]=-Q;  
   omegap[0]=sqrt((charge[0]*qm[0]*dens[0]*1.0)/EPS0); 
   vthx[0]=vthy[0]=vthz[0]=omegap[0]*debye[0]; //isothermal
   mass[0]=charge[0]/qm[0];  
   tempx[0]=mass[0]*vthx[0]*vthx[0]; 
   tempy[0]=mass[0]*vthy[0]*vthy[0];
   tempz[0]=mass[0]*vthz[0]*vthz[0];
   Npel=dens[0]*pow(debye[0],3)/ratio; //plasma param. for sim part, there is a ratio
  /*check debye and dx dy, and Np*/
   if(rank==0)
     {
	 //commented warnings are not necesarry
     //  if((debye[0]>0.1*dx)|| (debye[0]>0.1*dy) || (debye[0]>0.1*dz) )
	 //printf("WARNING, Debye length too big for given grid cell!!! Choose bigger dx or dy!!!%f \n", debye[0]/(0.1*dx));
     //  if((3*debye[0]<dx)|| (3*debye[0]<dy) || (3*debye[0]<dz))
	 //printf("WARNING, The grid cell is bigger than 3*Debyelenght!!! Choose smaller dx or dy!!!\n");
       if(dt*omegap[0]>0.1)
	 printf("WARNING, Integration time too large! Choose smaller dt!!!\n");
       if(Npel<10)
	 printf("WARNING, Plasma number: Nd>>1 is not fulfilled (Np<10), we have collisions which are not implemented in the code!\n");
     }
   if((numtasks>1) && (rank >0))
     {
       if(npart[0]>NPART_MAX)
	 {
	   if(rank==1)
	     {
	       printf("ERROR! Number of particles %ld  exceeds the maximum value: %d\n Aborting the simulation \n", npart[0], NPART_MAX);
	     }
	   exit(-1);
	 }  
     }
   else
     if(npart[0]>NPART_MAX)
       {
	 if(rank==0)
	   {
		printf("ERROR! Number of particles %ld  exceeds the maximum value: %d\n Aborting the simulation \n", npart[0], NPART_MAX);
	   }
	  exit(-1);
       }  
   
    /*calculating ion parameters*/
  /*if ions++, we need 0.5 ions for charge neutrality*/
#ifdef BEAM 
  npart[1]=(npart[0]/charge[1]-npart[2]);
#else
  npart[1]=npart[0]/charge[1]; 
#endif
  dens[1]=npart[1]*ratio/(Lx*Ly*Lz-dustvolume);  //real density

  if(numtasks>1)
    {
      npartinit[1]=(int)((npart[1]*1.0)/(numtasks-1));
      if((npart[1] % (numtasks-1))!= 0)
	if(rank==0) 
	  printf("Warning: no of ions != n*(numtasks-1), you have now %ld ions\n", npartinit[1]*(numtasks-1));
    }
  else
    npartinit[1]=npart[1]; //for the loop with different potenetials
 
  charge[1]=-charge[1]*charge[0]; //switch to real charge value
  mass[1]=charge[1]/qm[1];
  //omegap[1]=omegap[0]*sqrt(mass[0]/mass[1]);
  omegap[1]=sqrt((charge[1]*qm[1]*dens[1]*1.0)/EPS0); 
  debye[1]=debye[0]*sqrt(ti2te);
  vthx[1]=vthy[1]=vthz[1]=omegap[1]*debye[1]; //isothermal ions
  tempx[1]=mass[1]*vthx[1]*vthx[1]; 
  tempy[1]=mass[1]*vthy[1]*vthy[1];
  tempz[1]=mass[1]*vthz[1]*vthz[1];
  Npio=dens[1]*pow(debye[1],3)/ratio; //plasma parameter, sim part

 /************BEAM**************/
#ifdef BEAM
 /*calculating and printing ion beam*/
  /*BEAM BEAM*/
  dens[2]=npart[2]*1.0*ratio/(Lx*Ly*Lz-dustvolume);
  //necesarry for the loop with different potenetials
 if(numtasks>1)
    {
      npartinit[2]=(int)((npart[2]*1.0)/(numtasks-1));
      if((npart[2] % (numtasks-1))!= 0)
	if(rank==0) printf("Warning: no of beamions != n*(numtasks-1), you have now %d electrons\n", npartinit[2]*(numtasks-1));
    }
  else
    npartinit[2]=npart[2];

  charge[2]=charge[1]; //switch for real charge value
  mass[2]=mass[1];
  qm[2]=qm[1];
  //omegap[2]=omegap[0]*sqrt(mass[0]/mass[1]);
  //  debye[2]=debye[0]*sqrt(ti2te);
  omegap[2]=sqrt((charge[2]*qm[2]*dens[2]*1.0)/EPS0); 
  vthx[2]=vthy[2]=vthz[2]=0; //isothermal cold beam
  tempx[2]=0; 
  tempy[2]=0;
  tempz[2]=0;
  Npbe=dens[2]*pow(debye[2],3)/ratio; //plasma parameter for sim part
#endif
/**************BEAM END**********/

	
  /*change for simulation particles and add all particles*/
  for(i=0; i<S; i++)
    {    
      mass[i]*=ratio;
      charge[i]*=ratio;
      tempx[i]*=ratio;
      tempy[i]*=ratio;
	  tempz[i]*=ratio;
      dens[i]/=ratio;         	  
      allpart+=npart[i]; //on each task
    }  
  if(numtasks>1)
    allpart*=(numtasks-1);

  /*total debye length*/
  debyetotal=sqrt(debye[0]*debye[0]*debye[1]*debye[1]/(debye[0]*debye[0]+debye[1]*debye[1]));
  /*Calculate particles per cell, homogenous in space*/
  partpercell=allpart*dx*dy*dz/(Lx*Ly*Lz-dustvolume); //Number of cells used: volume = total - probe and divided by cellvolume
  /*calculate Cs*/
  Cs=sqrt((tempx[0]+tempx[1]*5/3)/mass[1]);
  //and assign drift velocity on each task
    for(i=0; i<S; i++)
    {    
      vdriftx[i]*=Cs; //on each task
    }  
  //FINISHED READING
  fclose(open);


  /*check if time step is short enough for given parameters*/
  //NOTE: important to assume proper Vpr_begin and Vpr_end .
  for(i=0; i<S; i++)
    {
      vxmax=sqrt(vdriftx[i]*vdriftx[i]+4*vthx[i]*vthx[i]+fabs(qm[i]*2*Vpr_begin));
      vxmax1=sqrt(vdriftx[i]*vdriftx[i]+4*vthx[i]*vthx[i]+fabs(qm[i]*2*Vpr_end));
      vxmax = (vxmax > vxmax1) ? vxmax : vxmax1;
      
      vymax=sqrt(4*vthy[i]*vthy[i]+fabs(qm[i]*2*Vpr_begin));
      vymax1=sqrt(4*vthy[i]*vthy[i]+fabs(qm[i]*2*Vpr_end));
      vymax = (vymax > vymax1) ? vymax : vymax1;
      
	    vzmax=sqrt(4*vthz[i]*vthz[i]+fabs(qm[i]*2*Vpr_begin));
      vzmax1=sqrt(4*vthz[i]*vthz[i]+fabs(qm[i]*2*Vpr_end));
      vzmax = (vzmax > vzmax1) ? vzmax : vzmax1;
	
	if(rank==0){  
	  printf("vdrift[%d]: %E and in Cs: %E\n", i,vdriftx[i],vdriftx[i]/Cs);
		printf("vxmax[%d]: dx/dt %E vs %E\n", i,dx/dt,vxmax);
		}
		
	if((dx/dt)<=(vxmax))
	{
	  if(rank==0)
	    {
	      printf("!!!dt too large, particles can jump over one cell in x dir\n!!!put dt < %E\n",dx/vxmax); 
	    }
	  exit(1);
	}	
      if((dy/dt)*1<=(vymax))
	{
	  if(rank==0)
	    {
	  printf("!!!dt too large, particles jump over one cell in y dir\n!!!put dt < %E\n", dy/vymax); 
	    }
	  exit(1);
	}
    
	  if((dz/dt)*1<=(vzmax))
	{
	  if(rank==0)
	    {
	  printf("!!!dt too large, particles jump over one cell in z dir\n!!!put dt < %E\n", dz/vzmax); 
	    }
	  exit(1);
	}
	
	
    }


  /*print the print.txt file*/  

  if(rank==0)
  {
    
  /*print electron parameters*/
       fprintf(write, "\n\nINPUT PARAMETERS for electrons (3D system) %d\n", 0);
       fprintf(write, "Charge               \t %E C\n", charge[0]/(ratio)); 
       fprintf(write, "Charge to mass ratio \t %E C/kg\n", qm[0]);  
       fprintf(write, "Mass                 \t %E kg\n", mass[0]/(ratio));
       fprintf(write, "Plasma frequency:    \t %E 1/s\n",omegap[0]);
       fprintf(write, "Debye length:        \t %E m\n", debye[0]);
       fprintf(write, "Thermal velocity x:  \t %E m/s\n", vthx[0]);
       fprintf(write, "Thermal velocity y:  \t %E m/s\n", vthy[0]);
       fprintf(write, "Thermal velocity z:  \t %E m/s\n", vthz[0]);
	   fprintf(write, "Temperature x:       \t %E eV\n", tempx[0]/(Q*ratio));
       fprintf(write, "Temperature y:       \t %E eV\n", tempy[0]/(Q*ratio));
	   fprintf(write, "Temperature z:       \t %E eV\n", tempz[0]/(Q*ratio));
	   fprintf(write, "Drift vel. in x dir: \t %E m/s (%E Cs)\n", vdriftx[0], vdriftx[0]/Cs);
       fprintf(write, "Density              \t %E 1/m3\n", dens[0]*ratio);
       fprintf(write, "Charge density       \t %E 1/m3\n", charge[0]*dens[0]);
       fprintf(write, "Plasma number, Np:   \t %E \n",Npel);
       fprintf(write, "No. of sim. particles\t %ld \n", npart[0]);
       fprintf(write, "Plasma/Simulation particles pr species:%E\n", ratio);
	
   /*print ion parameters*/
      fprintf(write, "\n\nINPUT PARAMETERS for ions (3D system) %d\n", 1);
      fprintf(write, "Charges              \t %E C \n", charge[1]/(ratio)); 
      fprintf(write, "Charge to mass ratio \t %E C/kg \n", qm[1]);  
      fprintf(write, "Mass                 \t %E kg\n", mass[1]/(ratio));
      fprintf(write, "Plasma frequency:    \t %E 1/s\n",omegap[1]);
      fprintf(write, "Debye length:        \t %E m\n", debye[1]);
      fprintf(write, "Thermal velocity x:  \t %E m/s\n", vthx[1]);
      fprintf(write, "Thermal velocity y:  \t %E m/s\n", vthy[1]);
      fprintf(write, "Thermal velocity z:  \t %E m/s\n", vthz[1]);
	   fprintf(write, "Temperature x:	\t %E eV\n", tempx[1]/(Q*ratio));
      fprintf(write, "Temperature y:       \t %E eV\n", tempy[1]/(Q*ratio));
      fprintf(write, "Temperature z:       \t %E eV\n", tempz[1]/(Q*ratio));
	   fprintf(write, "Drift vel. in x dir: \t %E m/s (%E Cs) \n", vdriftx[1], vdriftx[1]/Cs);
      fprintf(write, "Density              \t %E 1/m3\n", dens[1]*ratio);
      fprintf(write, "Charge density       \t %E 1/m3\n", charge[1]*dens[1]);
      fprintf(write, "Plasma number, Np:   \t %E \n",Npio);
      fprintf(write, "No. of sim. particles\t %ld \n", npart[1]);
      fprintf(write, "Plasma/Simulation particles pr species:\t%E \n", ratio);
 #ifdef BEAM   
      fprintf(write, "\n\nINPUT PARAMETERS for beam (3D system) %d\n", 2);
      fprintf(write, "Charges              \t %E C \n", charge[2]/(ratio)); 
      fprintf(write, "Charge to mass ratio \t %E C/kg \n", qm[2]);  
      fprintf(write, "Mass                 \t %E kg\n", mass[2]/(ratio));
      fprintf(write, "Plasma frequency:    \t %E 1/s\n",omegap[2]);
      fprintf(write, "Debye length:        \t %E m\n", debye[2]);
      fprintf(write, "Thermal velocity x:  \t %E m/s\n", vthx[2]);
      fprintf(write, "Thermal velocity y:  \t %E m/s\n", vthy[2]);
      fprintf(write, "Thermal velocity z:  \t %E m/s\n", vthz[2]);
	  fprintf(write, "Temperature x:       \t %E eV\n", tempx[2]/(Q*ratio));
      fprintf(write, "Temperature y:       \t %E eV\n", tempy[2]/(Q*ratio));
	  fprintf(write, "Temperature z:       \t %E eV\n", tempz[2]/(Q*ratio));
	  fprintf(write, "Drift vel. in x dir: \t %E m/s (%E Cs)\n", vdriftx[2], vdriftx[2]/Cs);
      fprintf(write, "Density              \t %E 1/m3\n", dens[2]*ratio);
      fprintf(write, "Charge density       \t %E 1/m3\n", charge[2]*dens[2]);
      fprintf(write, "Plasma number, Np:   \t %E \n",Npbe);
      fprintf(write, "No. of sim. particles\t %ld \n", npart[2]);
      fprintf(write, "Plasma/Simulation particles pr species:\t%E \n", ratio);
#endif
  
      fprintf(write,"\nParticles per cell %E\n", partpercell);
      partpercell=npartinit[0]*dx*dy*dz/(Lx*Ly*Lz-dustvolume); //volume = total - dust
      if(numtasks>1)
	    partpercell*=(numtasks-1);
      fprintf(write,"\nElectrons per cell %E\n", partpercell);
      fprintf(write,"Floating potential %E\n", tempx[0]*log(sqrt((tempx[0]*mass[1])/(tempx[1]*mass[0])))/charge[0]);

      printf("floating potential %E\n log %E\n", tempx[0]/charge[0], log(sqrt((tempx[0]*mass[1])/(tempx[1]*mass[0]))));
      fprintf(write,"Total debye length %E\n", debyetotal);
      fprintf(write, "Plasma number, Np  %E\n", debyetotal*debyetotal*debyetotal*dens[0]);
      fprintf(write,"Cs (adiabatic ions) %E [m/s]\n", Cs);
   
   
   
      fprintf(write,"PROBE potential: initial:%f V, final: %f V, step: %f V\n", Vpr_begin, Vpr_end, Vpr_step);
      fprintf(write, "\nBOUNDARY potential: %f V\n", Vbound);
      fprintf(write, "TOLERANCE for Poisson's equation solver: %E\n", TOLERANCE);
      fprintf(write, "Tolerance for the floating potential (conductor): %E\n", tolfloating);
      printf("Completed: readdata(),\n Input parameters are printed in file ''print.txt''.\n");
         
    fclose(write);
   }
   
  /*No particles on rank 0*/
  if((numtasks>1) && (rank==0))
    {
      for(i=0; i<S; i++)
	npartinit[i]=npart[i]=0;
    }

  /*initiate frequently used variables*/
  pi=M_PI;
  sqrt_two=sqrt(2);
  sqrt_twopi=sqrt(2*M_PI);
  sqrt_pi=sqrt(M_PI);  
	

  //this is for ??
		int ti;
		ti=0;	
  c0[ti]=c1[ti]=c2[ti]=c3[ti]=c4[ti]=c5[ti]=0;
  	ti=1;	
  c0[ti]=c1[ti]=c2[ti]=c3[ti]=c4[ti]=c5[ti]=0;
 } 
 
