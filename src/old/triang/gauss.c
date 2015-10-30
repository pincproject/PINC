/* DiP3D*/
/* gauss.c file */
/* Author Wojciech Jacek Miloch */
/* University of Oslo */
 /* 2009  */
 /*last revised 09.03.2009*/
#include "const.h"

/**** FIND ELECTRIC FIELD ************/
void electric_field(void)
{
  int i,j,k;
  for(i=1; i<ngx-1; i++) //For now we avoid boundaries
    for(j=1; j<ngy-1; j++)
	  for(k=1; k<ngz-1; k++)
      {
	//if(gp[i][j].marker != BOUND)
	  {
	    Fs[ix(0,i,j,k)]=(phi[ix(0,i-1,j,k)]-phi[ix(0,i+1,j,k)])/(2*dx*dVdt);
	    Fs[ix(FsEy,i,j,k)]=(phi[ix(0,i,j-1,k)]-phi[ix(0,i,j+1,k)])/(2*dy*dVdt);
		Fs[ix(FsEz,i,j,k)]=(phi[ix(0,i,j,k-1)]-phi[ix(0,i,j,k+1)])/(2*dz*dVdt);
	    /* multiplied by dt/dV to facilitate accel() and move() */
	  }
	//if(gp[i][j].marker == D_COND)
	  {    
	    //      printf("Efield %d %3.60E\n Potential %3.60E\n", ix(FsEy,i,j)-FsEy, Fs[ix(FsEy,i,j)], gp[i][j].phi);
	    //  if(t==0)
	    // Fs[ix(0,i,j)]=1.0/3;
	    // Fs[ix(0,i+1,j)]=1.0/3;  Fs[ix(0,i+2,j)]=2.0/3.0;
	    // printf("sqrt %1.60E\n", (Fs[ix(0,i,j)]));
	    // printf("sqrt %1.60E\n", (Fs[ix(0,i+1,j)]));
	    // printf("sqrt %1.60E\n", (Fs[ix(0,i+2,j)]));
	    //	getchar();	
	  }
      }
  
  //TREAT BOUNDARIES 
  //NEED TO CHANGE ITTTT!!!
  /*
  if(cond_present)
    {
      printf("electric_field(): I am sorry but this needs to be revised/changed\n");
      printf("now the field may be solved wrong\nBye\n"); 
      exit(1); 
      for(i=1; i<ngx-1; i++)
	for(j=1; j<ngy-1; j++)
	  if(gp[i][j].marker == PROBE)
	    {
	      if(gp[i-1][j].marker != PROBE) //left
		{
		  Fs[ix(0,i,j)]=2*Fs[ix(0,i-1,j)]-Fs[ix(0,i-2,j)]; 
		  if((gp[i][j-1].marker == PROBE) && (gp[i][j+1].marker == PROBE))
		    Fs[ix(FsEy,i,j)]=0;
		}
	      if(gp[i+1][j].marker != PROBE) //right
		{
		  Fs[ix(0,i,j)]=2*Fs[ix(0,i+1,j)]-Fs[ix(0,i+2,j)];   
		  if((gp[i][j-1].marker == PROBE) && (gp[i][j+1].marker == PROBE))
		    Fs[ix(FsEy,i,j)]=0;
		}
	      if(gp[i][j+1].marker != PROBE) //top
		{
		  Fs[ix(FsEy,i,j)]=2*Fs[ix(FsEy,i,j+1)]-Fs[ix(FsEy,i,j+2)];  
		  if((gp[i-1][j].marker == PROBE) && (gp[i+1][j].marker == PROBE))
		    Fs[ix(0,i,j)]=0;
		}
	      if(gp[i][j-1].marker != PROBE) //bottom
		{
		  Fs[ix(FsEy,i,j)]=2*Fs[ix(FsEy,i,j-1)]-Fs[ix(FsEy,i,j-2)];   
		  if((gp[i-1][j].marker == PROBE) && (gp[i+1][j].marker == PROBE))
		    Fs[ix(0,i,j)]=0;
		}
	    }  
    }
	*/
  //Boundaries linear extrapolation E field
  //vertical walls z - direction
  for(k=1; k<ngz-1; k++) 
  {
    for(j=1; j<ngy-1; j++)
      {
	  //here only in x direction
       Fs[ix(0,0,j,k)]=2*Fs[ix(0,1,j,k)]-Fs[ix(0,2,j,k)];
       Fs[ix(FsEy,0,j,k)]=Fs[ix(FsEy,ngx-1,j,k)]=Fs[ix(FsEz,ngx-1,j,k)]=Fs[ix(FsEz,0,j,k)]=0;
       Fs[ix(0,ngx-1,j,k)]=2*Fs[ix(0,ngx-2,j,k)]-Fs[ix(0,ngx-3,j,k)];
	   
      }
   for(i=1; i<ngx-1; i++)
     {
	 //here only in y direction
       Fs[ix(FsEy,i,0,k)]=2*Fs[ix(FsEy,i,1,k)]-Fs[ix(FsEy,i,2,k)];
       Fs[ix(FsEy,i,ngy-1,k)]=2*Fs[ix(FsEy,i,ngy-2,k)]-Fs[ix(FsEy,i,ngy-3,k)];
       Fs[ix(0,i,0,k)]=Fs[ix(0,i,ngy-1,k)]=Fs[ix(FsEz,i,0,k)]=Fs[ix(FsEz,i,ngy-1,k)]=0;
     }
  }	
  //horizontal walls
   for(i=1; i<ngx-1; i++)
     for(j=1; j<ngy-1; j++)
	 {
	  Fs[ix(FsEz,i,0,0)]=2*Fs[ix(FsEz,i,j,1)]-Fs[ix(FsEz,i,j,2)];
	  Fs[ix(FsEz,i,j,ngz-1)]=2*Fs[ix(FsEz,i,j,ngz-2)]-Fs[ix(FsEz,i,j,ngz-3)];
	  Fs[ix(0,i,j,0)]=Fs[ix(0,i,j,ngz-1)]=Fs[ix(FsEy,i,j,0)]=Fs[ix(FsEy,i,j,ngz-1)]=0;
	 }
	
  //corners
  Fs[ix(0,0,0,0)]=Fs[ix(0,0,ngy-1,0)]=Fs[ix(0,ngx-1,0,0)]=Fs[ix(0,ngx-1,ngy-1,0)]=Fs[ix(0,0,0,ngz-1)]=Fs[ix(0,0,ngy-1,ngz-1)]=Fs[ix(0,ngx-1,0,ngz-1)]=Fs[ix(0,ngx-1,ngy-1,ngz-1)]=0;
  Fs[ix(FsEy,0,0,0)]=Fs[ix(FsEy,0,ngy-1,0)]=Fs[ix(FsEy,ngx-1,0,0)]=Fs[ix(FsEy,ngx-1,ngy-1,0)]=Fs[ix(FsEy,0,0,ngz-1)]=Fs[ix(FsEy,0,ngy-1,ngz-1)]=Fs[ix(FsEy,ngx-1,0,ngz-1)]=Fs[ix(FsEy,ngx-1,ngy-1,ngz-1)]=0;
  Fs[ix(FsEz,0,0,0)]=Fs[ix(FsEz,0,ngy-1,0)]=Fs[ix(FsEz,ngx-1,0,0)]=Fs[ix(FsEz,ngx-1,ngy-1,0)]=Fs[ix(FsEz,0,0,ngz-1)]=Fs[ix(FsEz,0,ngy-1,ngz-1)]=Fs[ix(FsEz,ngx-1,0,ngz-1)]=Fs[ix(FsEz,ngx-1,ngy-1,ngz-1)]=0;
 
  //TEST
 // for(i=0; i<FsMAX; i++)
 //   Fs[i]=0.0; 
	//printf("Forces zero, FsMAX %d\n", FsMAX);
	//getchar();
}
