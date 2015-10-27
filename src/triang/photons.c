#include"const.h"
#include<stdlib.h>

void photonflux(void)
{
/*
  //angle must be in [0,360)
  //find how many b parameters need to be drawn each time-step for a givenangle
  double twoDarea,sinangle,cosangle;
  double bmin,bmax;
  double xmin,xmax;
  double angle, phflux;
  double length;
  double c;
  c=3E8;
  //299792458; //m/s
  srand(0);
  srand48(0);

  phflux=ph_flux; //normalized in grid.c
  //go from three dimensions to two
  double nf,nf2d, nph;
  nf=phflux/(c/normvel); //ok
  //  nf2d=nf*nf;
  nf2d=pow(nf, 0.6667);  //ok
  //two dimensional density is here per m2
  //now number of photons entering each time step the unitary box
  nph=nf2d*c*dt/normvel;
  printf("nf2d %E\n\n\n\n nph %E\n dt %E", nf2d/(normx*normx), nf/(normx*normx*normx), normvel);
  //now divide by ratio
  nph/=ratio;
  //  phflux/=normx*normtime;
  //  printf("normalized %E, %E dt %E\n", ph_flux, phflux, dt);
  angle=(ph_angle*2*M_PI)/360; //angle in radians
  sinangle=sin((ph_angle*2*M_PI)/360);
  cosangle=cos((ph_angle*2*M_PI)/360);
  ph_sinangle=sinangle;
  ph_cosangle=cosangle;
  ph_angle_rad=angle;
  // ph_a=sinangle;
  if((ph_angle<0) || (ph_angle>=360))
    {
      printf("Photons angle out of range! Finishing!\n");
      exit(1);
    } 

  if((ph_angle==90) || (ph_angle==270))
    {
      printf("Photon flux is along y axis, b parameter undefined\n");
      ph_bmin=0;
      ph_bmax=Lx;  
      length=Lx;
      ph_length=Lx;         
      twoDarea=length; //this is in real quantities
      ph_fluxprdt=(long int)(nph*twoDarea); //number of photons per dt through unitary box
      printf("in photonfluxprdt %ld a %E input ph_flux %E\n", ph_fluxprdt, ph_a, ph_flux);

 printf("END 90-270 in photonfluxprdt %ld a %E input ph_flux %E\n bmin %E length %E\n", ph_fluxprdt, ph_a, ph_flux/(normx*normx*normtime), ph_bmin, ph_length);
      getchar();
      if(ph_angle==90)
	{
	  ph_vert=1;
	    return;
	}      
      if(ph_angle==270)
	{
	  ph_vert=-1;
	return;
	}
      printf("Photons, Something is wrong in vertical movement\n");
      exit(1);
    }
  else
    {          
      //sections for generation of photons b parameter
      ph_a=sinangle/cosangle; 
      ph_vert=0;
      if(sinangle<0)
	{
	  bmin=0;
	  bmax=Ly+Lx*fabs(ph_a);
	}
      if(sinangle==0)
	{
	  printf("zero pierwsze Lx %E\n", Lx); getchar();
	  bmin=0;
	  bmax=Ly;
	}
      if(sinangle>0)
	{
	  bmin=-Lx*fabs(ph_a); 
	  bmax=Ly;
	}   
      
      //hitting
      //define where is x starts
      if((ph_angle<90) || (ph_angle>270))
	//we go towards positive x
	{
	  xmin=0;
	  xmax=Lx;
	}
      if((ph_angle>90) && (ph_angle<270))
	 //we go towards negative x
	{
	  xmin=Lx;
	  xmax=0;
	}
    }

  ph_bmin=bmin;
  ph_bmax=bmax;
  ph_xmin=xmin;
  ph_xmax=xmax;
  length=(bmax-bmin);
  printf("length %E sinangle %E\n", bmax-bmin, sinangle);
  if(sinangle==0)
    printf("zero\n");
  ph_length=length;

  twoDarea=fabs(length*cosangle); //this is in normalized quantities


  ph_fluxprdt=(long int)(nph*twoDarea); //number of photons per dt through unitary box
  
  //take account of MPI
  if(numtasks>1)
    {
      ph_fluxprdt=(long int)(nph*twoDarea)/(numtasks-1);
      if(rank==0)
	ph_fluxprdt=0;
    }
  
  

  printf("END in photonfluxprdt %ld a %E input ph_flux %E\n bmin %E length %E\n", ph_fluxprdt, ph_a, ph_flux, ph_bmin, ph_length);
  getchar();

   return;
   */
}

void photoelectriceffect(void)
{
/*
  //draw number ph_fluxprdt of particles from U(ph_bmin,ph_bmax)
  int counthv;
  counthv=0;

  double b;
  double phx1,phy1,phx2,phy2;

  FILE *assigned;
  int i,l;
  int j,k;
  double x,y,ax,ay;
  double wjk,wj1k,wj1k1,wjk1; //weights
  double qmratio;
  double partxold, partyold;
  int vert,kk;
  double q;
  int offset;
  int kp1;
  double partvxold, partvyold;
  
  int dno;
  double par,par_tmin;
  double licz, mian,odj;
  double t_tomove, cosphi, sinphi, tempx, tempy;
  double quad_a, quad_b, quad_c,quad_del, t_cross, t_cross1, t_cross2; 
  double xhitp, yhitp;
  int dnoseghit;
  FILE *assignedphel, *assignedphelvel;

  //  srand(0);
  //this is not normalized yet here!
  for(i=0; i<ph_fluxprdt; i++)
    {
      b=drand48()*ph_length+ph_bmin;

      if(ph_vert==0) //there is an angle
	{	  
	  phx1=ph_xmin;
	  phx2=ph_xmax;
	  phy1=b;
	  //	  phy1=-50;
	  phy2=b+ph_a*(phx2-phx1);
	  //  printf("b %E phy1 %E, phy2 %E\n",b, phy1,phy2);
	  //  getchar();
	}
      else //phvert==1
	{
	  phx1=phx2=b;	  
	  if(ph_vert==1) // going up
	    {	    
	      //      printf("going up\n");
	      phy1=0;
	      phy2=Ly;
	    }
	  if(ph_vert==-1) //travelling down
	    {
	      phy1=Ly;
	      phy2=0;
	      //      printf("going down\n");
	    }	  
	}
      partxold=phx1;
      partvxold=(phx2-phx1)/dt;      
      partyold=phy1;
      partvyold=(phy2-phy1)/dt; 
      
      //MOOVING DUST	
      //check all dusts and all corner pairs
      //initialization
      int dnohit;
      dnohit=-1;
      double timehitmin=-10*dt;
      timehitmin=-10*dt;
      dnoseghit=-1;
      //check now
      for(dno=0; dno<noofdusts; dno++)
	{
	  for(kk=0; kk<ncorners[dno]; kk++)
	    {
	      quad_a=(daa[dno][kk]-(partvxold*daa1y[dno][kk]+partvyold*daa1x[dno][kk]));
	      quad_b=(dbb[dno][kk]-(partxold*daa1y[dno][kk]+partyold*daa1x[dno][kk]+partvxold*dbb1y[dno][kk]+partvyold*dbb1x[dno][kk]));
	      quad_c=(dcc[dno][kk]-(partxold*dcc1y[dno][kk]-partyold*dcc1x[dno][kk]));
	      //	           	printf("k a:%E\t b%E\t %E\n", quad_a, quad_b, quad_c);
	      t_cross=-1; //false time
	      
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
		    t_cross=-1; //no cross
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
	      
	      if((t_cross <= dt) && (t_cross >=0)) 
		// the photon may to hit the dust 
		// now lets find the hitting parameter 
		{		    
		  kp1=kk+1;
		  if(kk==(ncorners[dno]-1))
		    kp1=0;
		  licz=mian=par=0;
		  //I am finding the parameter p here
		  if((daa1y[dno][kk]*t_cross+dbb1y[dno][kk])==0);
		  
		  if(fabs(-daa1x[dno][kk]*t_cross-dbb1x[dno][kk]) > fabs(daa1y[dno][kk]*t_cross+dbb1y[dno][kk]))
		    {
		      odj=dustxdxold[dno][kk]+dustvx[dno][kk]*t_cross;
		      licz=partxold+partvxold*t_cross-odj;
		      //mian=dustxdxold[dno][kp1]+dustvx[dno][kp1]*t_cross-odj;
		      //this is doing the same
		      mian=-daa1x[dno][kk]*t_cross-dbb1x[dno][kk];
		    }
		  else
		    {//only parameter for y
		      
		      odj=dustydyold[dno][kk]+dustvy[dno][kk]*t_cross;
		      licz=partyold+partvyold*t_cross-odj;		   
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
		  
		  //now I have the parameter calculated
		  par=licz/mian;
		  //now I need to make decision rule
		  if(par>=0 && par<=1.0)
		    //particle is close to the dust
		    {		   
		      //Ready to assign particle :)
		      {
			if((t_cross<timehitmin) || (timehitmin< -5*dt))
			  {	
			    dnohit=dno; 
			    dnoseghit=kk;
			    timehitmin=t_cross;
			    par_tmin=par;			       
			    //printf("hits %d, %E\n", l, t_cross);getchar();
			  }
		      }
		    }
		}
	    }
	}           
      //photoionization  now
      if(dnohit>=0)
	{ 
	  //	  printf("hit\n"); getchar();
	  //	      printf("assigning normal\n");
	  if(dnoseghit<0)
	    {
	      printf("SOMETHING WRONG in photons part (accel).c\n");
	      getchar();
	    } 
	  t_cross=timehitmin;
	  assigned=fopen("photoionzation.dat", "a");
	  xhitp=partxold+partvxold*t_cross;
	  yhitp=partyold+partvyold*t_cross;
	  
	  fprintf(assigned, "%E\t%E\n", xhitp*normx, yhitp*normx); //getchar();
	  fclose(assigned);
	  
	  //Create particle with the MC decision rule 
	  //TO BE INCLUDED
	  //NOW ALL PHOTONS that are HITTING DUST CREATE ELECTRONS


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
	  // particles normalized
	  //create one ion on the surface
	  int licznik;
	  if(dnoseghit==0)
	    {
	      licznik=0;
	    }
	  else
	    {
	      licznik=ccorner[dnohit][dnoseghit-1];
	    }
	  
	  dpart[dnohit][licznik+podloga].q+=normalcharge[1]*(1.0-weight_low)/ratio;
	  dpart[dnohit][licznik+sufit].q+=normalcharge[1]*weight_low/ratio;
	  // printf("normal charge %E\n", normalcharge[1]);
	  counthv++;
	  //particle current
	  current[dnohit][1]++; //CORRECT THE CURRENT ARRAY
	  double phe_energy, phevtot;
	  //CREATE NEW ELECTRON
	  //find energy of the photon...
	  if(ph_energy>100)
	    ;//use algorithm 
	  else
	    {
	      //easy rule 
	      //work function is here workfunct (eV)
	      phe_energy=ph_energy-dustworkfunct[dnohit];
	    }
	  //use real mass of the photon to calculate total energy [in J]
	  //and normalize it
	  phevtot=sqrt(2*phe_energy*Q/Mass_0)/normvel;
	  if(phevtot*dt>dx)
	    {
	      //	      printf("photoelectrons are too fast, veltot %E!!! \n", phevtot*normvel);
	      //getchar();
	    }
	  //find x,y,vx,vy
	  double phex,phey,phevx,phevy;
	  //find position in some distance along the ray out from the surface
	  double r, phermax;
	  phermax=phevtot*dt;
	  // r=1.0*drand48()*phermax;
	  r=drand48();
	  //  phex=xhitp+r*dustxnormv[dnohit][dnoseghit];
	  // phey=yhitp+r*dustynormv[dnohit][dnoseghit];
	  
	  //find direction of the velocity
	  double velangle,cosvelangle,sinvelangle, kat;
	  velangle=drand48()*M_PI-0.5*M_PI;
	  // velangle+=atan(dusta[dnohit][dnoseghit]);
	  cosvelangle=cos(velangle);
	  sinvelangle=sin(velangle);
	  //  velangle+=ph_angle_rad;
	 

	  //  printf("velangle %E\n", velangle);

	  //	  if(velangle>=2*M_PI)
	  //  velangle-=2*M_PI;
	  
	  //rotate normal vector to find the direction of the velocity
	  phevx=cosvelangle*dustxnormv[dnohit][dnoseghit]-sinvelangle*dustynormv[dnohit][dnoseghit];
	  phevy=sinvelangle*dustxnormv[dnohit][dnoseghit]+cosvelangle*dustynormv[dnohit][dnoseghit];


	  //FIND the forces acting on the photoelectron...
	  int je,ke;
	  double xe,ye;
	  double wjk,wj1k,wj1k1,wjk1;
	  double ax, ay;

	  //find forces from the closest grid points,and locate particle
	  je=xhitp/dx;
	  ke=yhitp/dy;
	  xe=xhitp-je*dx;
	  ye=yhitp-ke*dy;
	
	  //find weights, Ex,Ey is multipl.with dt/dxdy
	  wjk=(dx-x)*(dy-y);
	  wj1k=x*(dy-y);
	  wj1k1=x*y;
	  wjk1=(dx-x)*y;

	  //take only r fraction
	  ax=r*qm[0]*(Fs[ix(0,je,ke)]*wjk+Fs[ix(0,je+1,ke)]*wj1k+Fs[ix(0,je+1,ke+1)]*wj1k1+Fs[ix(0,je,ke+1)]*wjk1);
	  ay=r*qm[0]*(Fs[ix(FsEy,je,ke)]*wjk+Fs[ix(FsEy,je+1,ke)]*wj1k+Fs[ix(FsEy,je+1,ke+1)]*wj1k1+Fs[ix(FsEy,je,ke+1)]*wjk1);
	  //accelerate!
	  //	 	  printf("RELATIVE CHANGE x: %E, y: %E\n", ax/phevx, ay/phevy);
	  phevx+=ax;
	  phevy+=ay;		

	  //FIXING DIRECTION
	  //	  phevx=-phevtot;
	  // phevy=0;

    
	  //and move photoelectron now...
	  phex=xhitp+phevx*r*dt;
	  phey=yhitp+phevy*r*dt;
	  //phevx=-phevtot*cos(velangle);
	  // phevy=-phevtot*sin(velangle);
	  //  printf("velangle %E\n", velangle);
	  long int last;
	  last=npart[0];
	  spec[0].part[last].x=phex;
	  spec[0].part[last].y=phey;
	  //multiply the velocity with the value
	  spec[0].part[last].vx=phevtot*phevx;
	  spec[0].part[last].vy=phevtot*phevy;
	  npart[0]++;
	  //  printf("GENERATED %E ann %E pi %E rand %E\n", velangle, velangle*360/(2*M_PI), M_PI, drand48());
	  assignedphelvel=fopen("photoelectronsvel.dat", "a");
	  assignedphel=fopen("photoelectrons.dat", "a");      
	  fprintf(assignedphel,"%E\t%E\n",spec[0].part[last].x*normx,  spec[0].part[last].y*normx);
	  fprintf(assignedphelvel,"%E\t%E\n",spec[0].part[last].vx*normvel,spec[0].part[last].vy*normvel); 
	  fclose(assignedphel);
	  fclose(assignedphelvel);

	}     		
    }      
  printf("finished photons rank %d and %d\n", rank, counthv);
  */
}
