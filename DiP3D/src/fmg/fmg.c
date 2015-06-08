/* Based on Numerical Recipes in C. */
/* To solve eliptic PDE on 2D grid with the Full Multi Grid method*/
/* If not possible to use fmm, then use Gauss-Seidel method*/

#include <stdio.h>
#include<math.h>
#define NRANSI
#include "nrutil.h"
#include "../const.h"
//#include "funct.h"
#define NPRE 2
#define NPOST 3
//#define NGMAX 15
inline int ix(int off, int i, int j, int k);

// we have external
//int fmg_ng, fmg_nnx, fmg_nny, fmg_mingridx, fmg_mingridy; 
//double **ires[NGMAX+1],**irho[NGMAX+1],**irhs[NGMAX+1],**iu[NGMAX+1];

/*3D*/
void mglin_init(int nx, int ny, int nz)
{  
  int nn,ng,ng1,nn1,nnx,nny,nnz,ngrid,nn1x,nn1y,nn1z, k,j;
  
  ng=0;
  fmg_nnx=nx;
  fmg_nny=ny;
  fmg_nnz=nz;
    
  if(nx==ny)
    { 
      if(nx==nz)
	{
	  nn=nx;
	  fmg_mingridx=fmg_mingridy=fmg_mingridz=3;
	  while (nn >>= 1) ng++;
	  fmg_ng=ng;
	  printf("start ng %d nn %d comp %d\n", ng,nn, 1+(1L<<ng)); getchar();
	  if (nx != 1+(1L << ng)) 
	    {
	      printf("WARNING: cells number is not a power of 2 in mglin. Will use a non-optimal field solver!\n");
	    //find minimal grid
	      nn1=nx;
	      ng1=1;
	      while((nn1 % 2) == 1)
		{
		  ng1++; nn1=1+(nn1-1)/2;
		}
	      fmg_mingridx=fmg_mingridy=fmg_mingridz=nn1;
	      fmg_ng=ng1;
	    }
	}
    }
  else
    {
      printf("WARNING: cells number is not a power of 2 in mglin. Will use a non-optimal field solver!\n");
      printf("Tutaj algorytm chyba nie dziala i zaraz padnie!!\n");
      nn1=nx;
      ng1=1;
      while((nn1 % 2) == 1)
	{
	  ng1++; nn1=1+(nn1-1)/2;
	}
      ng=ng1;
      nn1=ny;
      ng1=1;
      while((nn1 % 2) == 1)
	{
	  ng1++; nn1=1+(nn1-1)/2;
	  }
      if(ng1 < ng)
	ng=ng1;  
      
      ng1=1;
      nn1=nz;
      while((nn1 % 2) == 1)
	{
	  ng1++; nn1=1+(nn1-1)/2;
	}
      if(ng1 < ng)
	ng=ng1;  


      nn1x=nx;
      nn1y=ny;
      nn1z=nz;
	
      for(k=0; k<ng-1; k++)		
	{
	  nn1x=1+(nn1x-1)/2;
	  nn1y=1+(nn1y-1)/2;
	  nn1z=1+(nn1z-1)/2;
	}
      
      fmg_mingridx=nn1x;
      fmg_mingridy=nn1y;
      fmg_mingridz=nn1z;

      fmg_nnx=nx;
      fmg_nny=ny;
      fmg_nnz=nz;

      nnx=nx;
      nny=ny;
      nnz=nz;

      fmg_ng=ng;
    }
  printf("mglin init: ng %d mingrid %d and %d and %d  nnx %d %d %d\n", fmg_ng, fmg_mingridx, fmg_mingridy, fmg_mingridz, fmg_nnx, fmg_nny, fmg_nnz);  
  if (fmg_ng > NGMAX) nrerror("increase NGMAX in mglin.");
  
  // allocate memory for grid arrays
  if(fmg_ng>1)
    {
      nnx=fmg_nnx/2+1;
      nny=fmg_nny/2+1;
      nnz=fmg_nnz/2+1;
      ngrid=fmg_ng-1;			
      irho[ngrid]=d3tensor(0,nnx-1,0,nny-1,0,nnz-1);
		printf("rho ngrid %d : %E \n", ngrid, irho[ngrid][1][1][1]);
      while (nnx > fmg_mingridx) {
	nnx=nnx/2+1;
	nny=nny/2+1;
	nnz=nnz/2+1;
	irho[--ngrid]=d3tensor(0,nnx-1,0,nny-1,0,nnz-1);
      }
      
      nnx=fmg_mingridx;
      nny=fmg_mingridy;
      nnz=fmg_mingridz;
      iu[1]=d3tensor(0,nnx-1,0,nny-1,0,nnz-1);
      irhs[1]=d3tensor(0,nnx-1,0,nny-1,0,nnz-1);
      printf("rhs hmm %E \n", irhs[1][0][0][0]);
		
      for (j=2;j<=fmg_ng;j++) {
	nnx=2*nnx-1;
	nny=2*nny-1;
	nnz=2*nnz-1;
	//	printf("%d %d %d for %d\n", nnx, nny, nnz, j);
	//	getchar();
	iu[j]=d3tensor(0,nnx-1,0,nny-1,0,nnz-1);
	irhs[j]=d3tensor(0,nnx-1,0,nny-1,0,nnz-1);
	ires[j]=d3tensor(0,nnx-1,0,nny-1,0,nnz-1);
	//	printf("rhs %d \n", irhs[j][0][0][0]);	
      }	
    }
  
}

/*3D*/
void mglin_destroy()
{
  int nnx,nny,nnz,j;
  if(fmg_ng>1)
    {
      for(nnx=fmg_nnx,nny=fmg_nny,nnz=fmg_nnz,j=fmg_ng;j>=2;j--,nnx=nnx/2+1,nny=nny/2+1,nnz=nnz/2+1) 
	{
	  free_d3tensor(ires[j],0,nnx-1,0,nny-1,0,nnz-1);
	  free_d3tensor(irhs[j],0,nnx-1,0,nny-1,0,nnz-1);
	  free_d3tensor(iu[j],0,nnx-1,0,nny-1,0,nnz-1);
	  if (j != fmg_ng) free_d3tensor(irho[j],0,nnx-1,0,nny-1,0,nnz-1);
	}
      free_d3tensor(irhs[1],0,fmg_mingridx-1,0,fmg_mingridy-1,0,fmg_mingridz-1);
      free_d3tensor(iu[1],0,fmg_mingridx-1,0,fmg_mingridy-1,0,fmg_mingridz-1);
      free_d3tensor(irho[1],0,fmg_mingridx-1,0,fmg_mingridy-1,0,fmg_mingridz-1);
    }
}


void mglin(double *u, int ncycle)
{
  void addint(double ***uf, double ***uc, double ***res, int nfx, int nfy, int nfz);
  void copy(double ***aout, double ***ain, int nx, int ny, int nz);
  void copy0(double ***aout, double *ain, int nx, int ny, int nz);
  void copyfinal(double *aout, double ***ain, int nx, int ny, int nz);
  void fill0(double ***u, int nx, int ny, int nz);
  void interp(double ***uf, double ***uc, int nfx, int nfy, int nfz);
  void relax(double ***u, double ***rhs, int nx, int ny, int nz);
  void resid(double ***res, double ***u, double ***rhs, int nx, int ny, int nz);
  void rstrct(double ***uc, double ***uf, int ncx, int ncy, int ncz);
  void rstrct0(double ***uc, double *uf, int ncx, int ncy, int ncz);
  void slvsml(double ***u, double ***rhs);
  void slvsml2(double ***u, double ***rhs, int nx, int ny, int nz);
  
  unsigned int j,jcycle,jj,jpost,jpre,ng=0,ngrid;
  unsigned int  nnx, nny, nnz, nfx,nfy,nfz;
  
  printf("Solving field equations in mglin\n");

  /*initialize fmg*/
  ng=fmg_ng;
  nnx=fmg_nnx;
  nny=fmg_nny;
  nnz=fmg_nnz;
  
  if(ng>1)
    {
      nnx=nnx/2+1;
      nny=nny/2+1;
      nnz=nnz/2+1;
      ngrid=ng-1;			
	  printf("parameters ngrid %d %d %d %d\n", ngrid, nnx, nny, nnz);
	//	printf("parameters irho ngrid %E\n", irho[6][1][1][1]);
      rstrct0(irho[ngrid],u,nnx,nny,nnz);
 
      /*TUTAJ*/
 
      while (nnx > fmg_mingridx) {
	nnx=nnx/2+1;
	nny=nny/2+1;
	nnz=nnz/2+1;
	--ngrid;
	rstrct(irho[ngrid],irho[ngrid+1],nnx,nny,nnz);
		
      }
      nnx=fmg_mingridx;
      nny=fmg_mingridy;
      nnz=fmg_mingridz;
		
      if((fmg_mingridx==3) && (fmg_mingridy==3) && (fmg_mingridz==3))
	{slvsml(iu[1],irho[1]);  printf("I am using the full multigrid method\n");}
      else
	{
	  slvsml2(iu[1],irho[1],fmg_mingridx,fmg_mingridy,fmg_mingridz); 
	  printf("I am using the semi - multigrid method (not the best choice of grid size)\n");
	}
      ngrid=ng;
    
      for (j=2;j<=ngrid;j++) {
	nnx=2*nnx-1;
	nny=2*nny-1;
	nnz=2*nnz-1;
	printf("now interp?\n");
       	interp(iu[j],iu[j-1],nnx,nny,nnz);
	//	printf("rhs %d \n", irhs[2][0][0][0]); getchar();
	printf("I got here j: %d and ng %d\n",j,ng);
	if(j!=ngrid)
	  {
	    //  printf("try to copy %E\n", irhs[j][0][0][0]);
	  copy(irhs[j], irho[j], nnx,nny,nnz);
	  }
	else
	  {
	    // printf("j=%d\n",j);
	    copy0(irhs[j], u,nnx,nny,nnz);	
	    
	  }
	for (jcycle=1;jcycle<=ncycle;jcycle++) {
	  nfx=nnx;
	  nfy=nny;
	  nfz=nnz;
	  for (jj=j;jj>=2;jj--) {
	    for (jpre=1;jpre<=NPRE;jpre++)
	      relax(iu[jj],irhs[jj],nfx,nfy,nfz);
	    
	    resid(ires[jj],iu[jj],irhs[jj],nfx,nfy,nfz);
	    nfx=nfx/2+1;
	    nfy=nfy/2+1;
	    nfz=nfz/2+1;
	    rstrct(irhs[jj-1],ires[jj],nfx,nfy,nfz);
	    fill0(iu[jj-1],nfx,nfy,nfz);			
	  }
	  if((fmg_mingridx==3) && (fmg_mingridy==3) && (fmg_mingridz==3))
	    slvsml(iu[1],irhs[1]);
	  else
	    slvsml2(iu[1],irhs[1],fmg_mingridx,fmg_mingridy,fmg_mingridz);
	  nfx=fmg_mingridx;
	  nfy=fmg_mingridy;
	  nfz=fmg_mingridz;
	  for (jj=2;jj<=j;jj++) {
	    nfx=2*nfx-1;
	    nfy=2*nfy-1;
	    nfz=2*nfz-1;
	    addint(iu[jj],iu[jj-1],ires[jj],nfx,nfy,nfz);
	    
	    for (jpost=1;jpost<=NPOST;jpost++)
	      relax(iu[jj],irhs[jj],nfx,nfy,nfz);
	  }
	}
      }
      copyfinal(phi,iu[ngrid],fmg_nnx,fmg_nny,fmg_nnz);
    }
  else
{
  int i,ipass, isw,j,jsw, sweep, k;
  int iter=0;
  double h2,hx,hy,hz,error,errorcheck,toler, factor;  
  jsw=1;
  toler=0.000001;
  printf("Wrong grid size to use the multigrid method\nI am solving the potential with the slow converging Gauss-Seidel method\n");
  
  hx=Lx/(fmg_nnx-1);
  hy=Ly/(fmg_nny-1);
  hz=Lz/(fmg_nnz-1);
  
  h2=hx*hy;

  factor=1.0/6;
  do
    {
      iter++;
      errorcheck=0.5*toler;
      for (ipass=1;ipass<=2;ipass++,jsw=3-jsw) {
	sweep=jsw;	
	for (i=1;i<fmg_nnx-1;i++)
	  for (j=1;j<fmg_nny-1;j+=1)
	    for(k=(2-(i+j+sweep)%2); k<fmg_nnz-1; k+=2)	      
	      {
		error=phi[ix(0,i,j,k)];
		phi[ix(0,i,j,k)]=factor*(phi[ix(0,i+1,j,k)]+phi[ix(0,i-1,j,k)]+phi[ix(0,i,j+1,k)]+phi[ix(0,i,j-1,k)]+phi[ix(0,i,j,k-1)]+phi[ix(0,i,j,k+1)]+h2*u[ix(0,i,j,k)]);   
		//we have plus sign here, because the rho is positive here
		error=phi[ix(0,i,j,k)]-error;
		if(errorcheck < fabs(error))
		  errorcheck=fabs(error);				
	      }
      }
    }
  while(errorcheck>toler);
     
  printf("END %d iterations\n", iter);
}
}
#undef NPRE
#undef NPOST
#undef NGMAX
#undef NRANSI

/*3D*/
void rstrct(double ***uc, double ***uf, int ncx, int ncy, int ncz)
{
  int ic,iif,jc,jf,kc,kf,nccx, nccy, nccz;
  /*set ncc to be the dimension of the larger grid*/
  nccx = 2*ncx-2; /*previously*/
  nccy = 2*ncy-2; /*previously*/
  nccz = 2*ncz-2;
  double wf=1.0/12;

  for(kf=2,kc=1;kc<ncz-1;kc++,kf+=2)
    for (jf=2,jc=1;jc<ncy-1;jc++,jf+=2){
      for (iif=2,ic=1;ic<ncx-1;ic++,iif+=2){
	uc[ic][jc][kc]=0.5*uf[iif][jf][kf]+wf*(uf[iif+1][jf][kf]+uf[iif-1][jf][kf]+uf[iif][jf+1][kf]+uf[iif][jf-1][kf]+uf[iif][jf][kf+1]+uf[iif][jf][kf-1]);
      }
    }

  /*boundary points*/
  /*planes*/
  /*do bottom and top plane*/
  for (jf=2,jc=1;jc<ncy-1;jc++,jf+=2){
    for (iif=2,ic=1;ic<ncx-1;ic++,iif+=2){
      uc[ic][jc][0]=uf[iif][jf][0];
      uc[ic][jc][ncz-1]=uf[iif][jf][nccz];
    }
  }
  /*do face and back (depth -> y)*/
  for (kf=2,kc=1;kc<ncz-1;kc++,kf+=2){
    for (iif=2,ic=1;ic<ncx-1;ic++,iif+=2){
      uc[ic][0][kc]=uf[iif][0][kf];
      uc[ic][ncz-1][kc]=uf[iif][nccy][kf];
    }
  }
  /*do left and right (depth -> y)*/
  for (kf=2,kc=1;kc<ncz-1;kc++,kf+=2){
    for (jf=2,jc=1;jc<ncy-1;jc++,jf+=2){
      uc[0][jc][kc]=uf[0][jf][kf];
      uc[ncx-1][jc][kc]=uf[nccx][jf][kf];
    }
  }

  /*do edges*/
  //bottom and top along x
  for (jc=0,ic=0;ic<ncx;ic++,jc+=2) {
    uc[ic][0][0]=uf[jc][0][0];
    uc[ic][ncy-1][0]=uf[jc][nccy][0];
    uc[ic][0][ncz-1]=uf[jc][0][nccz];
    uc[ic][ncy-1][ncz-1]=uf[jc][nccy][nccz];
  }
  
  //bottom and top along y
  for (jc=0,ic=0;ic<ncy;ic++,jc+=2) {
    uc[0][ic][0]=uf[0][jc][0];
    uc[ncx-1][ic][0]=uf[nccx][jc][0];
    uc[0][ic][ncz-1]=uf[0][jc][nccz];
    uc[ncx-1][ic][ncz-1]=uf[nccx][jc][nccz];    
  }
  
  //vertically 
  for (jc=0,ic=0;ic<ncz;ic++,jc+=2) {
    uc[0][0][ic]=uf[0][0][jc];
    uc[ncx-1][0][ic]=uf[nccx][0][jc];
    uc[0][ncy-1][ic]=uf[0][nccy][jc];
    uc[ncx-1][ncy-1][ic]=uf[nccx][nccy][jc];    
  }
}
/*************************************************/
/*3D*/
void rstrct0(double ***uc, double *uf, int ncx, int ncy, int ncz)
{
  int ic,iif,jc,jf,kc,kf,nccx, nccy, nccz;
  /*set ncc to be the dimension of the larger grid*/
  nccx = 2*ncx-2; /*previously*/
  nccy = 2*ncy-2; /*previously*/
  nccz = 2*ncz-2;
  double wf=1.0/12;

  for(kf=2,kc=1;kc<ncz-1;kc++,kf+=2)
    for (jf=2,jc=1;jc<ncy-1;jc++,jf+=2){
      for (iif=2,ic=1;ic<ncx-1;ic++,iif+=2){
		  //	printf("%E\t%E\n",uf[ix(0,iif,jf,kf)], uc[ic][jc][kc]);
	//	  printf("%d %d %d %E\n", ic, jc, kc, uf[ix(0,iif,jf,kf)]);	  
	uc[ic][jc][kc]=0.5*uf[ix(0,iif,jf,kf)]+wf*(uf[ix(0,iif+1,jf,kf)]+uf[ix(0,iif-1,jf,kf)]+uf[ix(0,iif,jf+1,kf)]+uf[ix(0,iif,jf-1,kf)]+uf[ix(0,iif,jf,kf+1)]+uf[ix(0,iif,jf,kf-1)]);
	//	printf("%E\t%E\n",uf[ix(0,iif,jf,kf)], uc[ic][jc][kc]);
      }
    }
      printf("buss errrror byl here\n");
  /*boundary points*/
  /*planes*/
  /*do bottom and top plane*/
  for (jf=2,jc=1;jc<ncy-1;jc++,jf+=2){
    for (iif=2,ic=1;ic<ncx-1;ic++,iif+=2){
      uc[ic][jc][0]=uf[ix(0,iif,jf,0)];
      uc[ic][jc][ncz-1]=uf[ix(0,iif,jf,nccz)];
    }
  }
  /*do face and back (depth -> y)*/
  for (kf=2,kc=1;kc<ncz-1;kc++,kf+=2){
    for (iif=2,ic=1;ic<ncx-1;ic++,iif+=2){
      uc[ic][0][kc]=uf[ix(0,iif,0,kf)];
      uc[ic][ncz-1][kc]=uf[ix(0,iif,nccy,kf)];
    }
  }
  /*do left and right (depth -> y)*/
  for (kf=2,kc=1;kc<ncz-1;kc++,kf+=2){
    for (jf=2,jc=1;jc<ncy-1;jc++,jf+=2){
      uc[0][jc][kc]=uf[ix(0,0,jf,kf)];
      uc[ncx-1][jc][kc]=uf[ix(0,nccx,jf,kf)];
    }
  }

  /*do edges*/
  //bottom and top along x
  for (jc=0,ic=0;ic<ncx;ic++,jc+=2) {
    uc[ic][0][0]=uf[ix(0,jc,0,0)];
    uc[ic][ncy-1][0]=uf[ix(0,jc,nccy,0)];
    uc[ic][0][ncz-1]=uf[ix(0,jc,0,nccz)];
    uc[ic][ncy-1][ncz-1]=uf[ix(0,jc,nccy,nccz)];
  }
  
  //bottom and top along y
  for (jc=0,ic=0;ic<ncy;ic++,jc+=2) {
    uc[0][ic][0]=uf[ix(0,0,jc,0)];
    uc[ncx-1][ic][0]=uf[ix(0,nccx,jc,0)];
    uc[0][ic][ncz-1]=uf[ix(0,0,jc,nccz)];
    uc[ncx-1][ic][ncz-1]=uf[ix(0,nccx,jc,nccz)];    
  }
  
  //vertically 
  for (jc=0,ic=0;ic<ncz;ic++,jc+=2) {
    uc[0][0][ic]=uf[ix(0,0,0,jc)];
    uc[ncx-1][0][ic]=uf[ix(0,nccx,0,jc)];
    uc[0][ncy-1][ic]=uf[ix(0,0,nccy,jc)];
    uc[ncx-1][ncy-1][ic]=uf[ix(0,nccx,nccy,jc)];    
  }

	int i,j,k;
	for(i=0;i<ncx;i++)
	  for(j=0;j<ncy;j++)
	    for(k=0; k<ncz; k++)
	      {
		//	printf("in loop\n");
		//printf("%E\n", uc[i][j][k]);
		//	aout[i][j][k]=ain[i][j][k];
		//	printf("after %d %d\n", ain[i][j][k], aout[i][j][k]);
	      }
  //copy(uc, uc, ncx, ncy, ncz);
}
/*************************************************/

/*3D ??? test indices here else ti should be OK*/
void interp(double ***uf, double ***uc, int nfx, int nfy, int nfz)
{
  int ic,iif,jc,jf,kf,kc,ncx,ncy,ncz,jk;
  ncx=nfx/2+1;
  ncy=nfy/2+1;
  ncz=nfz/2+1;

  /*do copies*/  
  for(kc=0, kf=0; kc<ncz; kc++, kf+=2) //kc kf
    for (jc=0,jf=0;jc<ncy;jc++,jf+=2) //jc jf
      for (ic=0;ic<ncx;ic++)  //ic
	{
	
  uf[2*ic][jf][kf]=uc[ic][jc][kc];
  //  printf("copies %d %d %d %d %d %d %E %E\n", 2*ic, jf, kf, ic, jc, kc, uf[2*ic][jf][kf],uc[ic][jc][kc]);	
	}
  // printf("in interp OK\n");
  /*interp every second slice*/
  for(kf=0; kf<nfz; kf+=2)
    {
      for(jf=0;jf<nfy;jf+=2)
	for(iif=1;iif<nfx-1;iif+=2)
	  uf[iif][jf][kf]=0.5*(uf[iif+1][jf][kf]+uf[iif-1][jf][kf]);
      for(jf=1;jf<nfy-1;jf+=2)
	for(iif=0;iif <nfx;iif++)
	  uf[iif][jf][kf]=0.5*(uf[iif][jf+1][kf]+uf[iif][jf-1][kf]);
    }
  /*interpolate betwen slices in depth*/
  for(kf=1; kf<nfz-1; kf+=2)
    {
      for(jf=0;jf<nfy;jf++)
	for(iif=0; iif<nfx; iif++)
	  uf[iif][jf][kf]=0.5*(uf[iif][jf][kf-1]+uf[iif][jf][kf+1]);
    }
	int i,j,k;
	for(i=0;i<nfx;i++)
	  for(j=0;j<nfy;j++)
	    for(k=0; k<nfz; k++)
	      {
		//	printf("in loop\n");
		//	printf("%E\n", uf[i][j][k]);
		//	aout[i][j][k]=ain[i][j][k];
		//	printf("after %d %d\n", ain[i][j][k], aout[i][j][k]);
	      }
	//	printf("END interp\n");
}
/**************************************************/
/*3D*/
void addint(double ***uf, double ***uc, double ***res, int nfx, int nfy, int nfz)
{
  void interp(double ***uf, double ***uc, int nfx, int nfy, int nfz);
  int i,j,k;
  
  interp(res,uc,nfx,nfy,nfz);
  for(i=0;i<nfx;i++)
    for (j=0;j<nfy;j++)
      for(k=0;k<nfz;k++)
	uf[i][j][k] += res[i][j][k];
}
/**************************************************/
/*3D*/
void slvsml(double ***u, double ***rhs)
{	
  void fill0(double ***u, int nx, int ny, int nz);
  double h=Lx*0.5;
  fill0(u,3,3,3);
  u[1][1][1] = h*h*rhs[1][1][1]/6.0;
}

/**************************************************/
/*3D*/
void slvsml2(double ***u, double ***rhs, int nx, int ny, int nz)
{
  int i,ipass,isw,j,jsw=1,k, sweep;
  double hx,hy,h2;
  double toler, factor, error, errorcheck;
  void fill0(double ***u, int nx, int ny, int nz);
  int iter=0;
  fill0(u,nx,ny,nz);
  hx=Lx/(nx-1);
  hy=Ly/(ny-1);
  h2=hx*hy;
  toler=0.00000001;
  factor=1.0/6;
  do
    {
      iter++;
      errorcheck=0.5*toler;
      for (ipass=1;ipass<=2;ipass++,jsw=3-jsw) {
	sweep=jsw;	
	for (i=1;i<nx-1;i++)
	  for (j=1;j<ny-1;j+=1)
	    for(k=(2-(i+j+sweep)%2); k<nz-1; k+=2)	      
	      {
		error=u[i][j][k];
		u[i][j][k]=factor*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k]+u[i][j][k-1]+u[i][j][k+1]+h2*rhs[i][j][k]);
		error=u[i][j][k]-error;
		if(errorcheck < fabs(error))
		  errorcheck=fabs(error);				
	      }
      }
    }
  while(errorcheck>toler);
}
/**************************************************/
/*3D*/
void relax(double ***u, double ***rhs, int nx, int ny, int nz)
{
  int i,ipass,isw,j,jsw=1,k, sweep;
  double hx,hy,h2,factor;
  hx=Lx/(nx-1);
  hy=Ly/(ny-1);
  h2=hx*hy;
  factor=1.0/6;

  //go red first
  for (ipass=1;ipass<=2;ipass++,jsw=3-jsw) {
    sweep=jsw; 
    for (i=1;i<nx-1;i++)
      for (j=1;j<ny-1;j+=1)
	for(k=(2-(i+j+sweep)%2); k<nz-1; k+=2)
	  u[i][j][k]=factor*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k]+u[i][j][k-1]+u[i][j][k+1]+h2*rhs[i][j][k]);
  }
}

/**************************************************/
/*3D*/
void resid(double ***res, double ***u, double ***rhs, int nx, int ny, int nz)
{
  int i,j,k,jf,jc, iif, ic, kf, kc;
  double hx,hy,h2i;
  
  hx=Lx/(nx-1);
  hy=Ly/(ny-1);	
  h2i=1.0/(hx*hy);

  /*interior*/
  for (i=1;i<nx-1;i++)
    for (j=1;j<ny-1;j++)
      for(k=1;k<nz-1;k++)
	res[i][j][k] = h2i*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k]+u[i][j][k+1]+u[i][j][k-1]-6.0*u[i][j][k])+rhs[i][j][k];

  /*boundary points*/
  // for (i=0;i<nx;i++)
  //   res[i][0]=res[i][ny-1]=0.0;
  //  for (i=0;i<ny;i++)
  //   res[0][i]=res[nx-1][i]=0.0;

  /*boundary points*/
  /*planes*/
  /*do bottom and top plane*/
  for (jf=2,jc=1;jc<ny-1;jc++,jf+=2){
    for (iif=2,ic=1;ic<nx-1;ic++,iif+=2){
      res[ic][jc][0]=0.0;
      res[ic][jc][nz-1]=0.0;
    }
  }
  /*do face and back (depth -> y)*/
  for (kf=2,kc=1;kc<nz-1;kc++,kf+=2){
    for (iif=2,ic=1;ic<nx-1;ic++,iif+=2){
      res[ic][0][kc]=0.0;
      res[ic][nz-1][kc]=0.0;
    }
  }
  /*do left and right (depth -> y)*/
  for (kf=2,kc=1;kc<nz-1;kc++,kf+=2){
    for (jf=2,jc=1;jc<ny-1;jc++,jf+=2){
      res[0][jc][kc]=0.0;
      res[nx-1][jc][kc]=0.0;
    }
  }

  /*do edges*/
  //bottom and top along x
  for (jc=0,ic=0;ic<nx;ic++,jc+=2) {
    res[ic][0][0]=0.0;
    res[ic][ny-1][0]=0.0;
    res[ic][0][nz-1]=0.0;
    res[ic][ny-1][nz-1]=0.0;
  }
  
  //bottom and top along y
  for (jc=0,ic=0;ic<ny;ic++,jc+=2) {
    res[0][ic][0]=0.0;
    res[nx-1][ic][0]=0.0;
    res[0][ic][nz-1]=0.0;
    res[nx-1][ic][nz-1]=0.0;    
  }
  
  //vertically 
  for (jc=0,ic=0;ic<nz;ic++,jc+=2) {
    res[0][0][ic]=0.0;
    res[nx-1][0][ic]=0.0;
    res[0][ny-1][ic]=0.0;
    res[nx-1][ny-1][ic]=0.0;    
  }
}


/***************************************************/
/*3D*/
void copy(double ***aout, double ***ain, int nx, int ny, int nz)
{
  //  printf("in copy %E %E\n", ain[0][ny-1][nz-1], ain[0][ny-1][nz-1]);  
  //aout[1][ny-1][nz-1]=ain[1][ny-1][nz-1];  

  //  printf("in copy %E %E\n", aout[1][ny-1][nz-1], ain[1][ny-1][nz-1]);  
  int i,j,kk;
  for(i=0;i<nx;i++)
    for(j=0;j<ny;j++)
      for(kk=0; kk<nz; kk++)
	{
	  //	  printf("in loop ijk %d %d %d\n", i,j,kk);
	  //  printf("%E %E\n", ain[i][j][kk], aout[i][j][kk]);
	  aout[i][j][kk]=ain[i][j][kk];
	  //	printf("after %E %E\n", ain[i][j][k], aout[i][j][k]);
	}
  //  printf("finished\n");
}
/***************************************************/
/*3D*/
void copy0(double ***aout, double *ain, int nx, int ny, int nz)
{
  int i,j,k;
  for (i=0;i<nx;i++)
    for (j=0;j<ny;j++)
      for(k=0;k<nz;k++){
	//	printf("in copy %E %E\n", ain[ix(0,i,j,k)], aout[i][j][k]);
	    aout[i][j][k]=ain[ix(0,i,j,k)];}
  printf("copied 0\n");
  //minus sign is needed cause we have L*u=rhs, and our rhs should be -rho			
}
/***************************************************/
/*3D*/
void copyfinal(double *aout, double ***ain, int nx, int ny, int nz)
{
  int i,j,k;                    
  for (i=0;i<nx;i++)
    for (j=0;j<ny;j++)
      for(k=0;k<nz;k++)
	    aout[ix(0,i,j,k)]=ain[i][j][k];
}
/***************************************************/
/*3D*/
void fill0(double ***u, int nx, int ny, int nz)
{
  int i,j,k;
  for(i=0;i<nx;i++)
    for(j=0;j<ny;j++)
      for(k=0;k<nz;k++)
    	u[i][j][k]=0.0;
}
/***************************************************/


