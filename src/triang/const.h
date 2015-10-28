/* DiP3D program */
/* Definition of global variables and constants*/
/* Author: Wojciech Jacek Miloch */
/* University of Oslo */
/* July 2008 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "funct.h"

/*Parallel version library*/
#ifdef MPI
#include "mpi.h"
#endif

/*Main physical constants*/
#define EPS0 8.854E-12 /*epsilon zero,C^2/(N*m^2)*/
#define Mass_0 9.1095E-31 /*electron mass, kg*/
#define Mass_1 1.6726E-27 /*proton mass, kg*/
#define Q 1.602189E-19 /*elementary charge value, C*/
/*Definition of the upper limits*/
#define Lx_MAX 5 /*Maximum dimensions in meters*/
#define Ly_MAX 5
#define Lz_MAX 5

#define ngx_MAX 1024 /*Maximum number of grid points*/
#define ngy_MAX 1024
#define ngz_MAX 1024

#define NPART_MAX 10000000 /*Max no. of particles,each species per process*/

#define DIM 3 /*dimesnions*/
#define NOF 1 /*no of forces E=1, E+B=2*/

#ifdef BEAM
#define S 3
#else
#define S 2 /*number of species*/
#endif

/*Here are particles defined, index 0 electrons, 1 ions*/
#define electron 0
#define ion 1
#define beam 2
#define GONE -1

/*markers for the grid*/
#define NORMAL 0
#define BOUND 1
#define PROBE 2

#define D_INSU 3
#define D_COND 2
 
#define CONVTEST 5
#define PRSEG 4 //number of probesegments
#define CROSSFACTOR 0.00001 //CROSSFACTOR * dx -> accuracy for dust hitting

#define LIST_SIZE 1500

double TOLERANCE; //convergence limit for Poisson solver
double tolfloating; //toler. for float. pot.,condƒtor mode

/*particles*/
double dens[S]; /*number density*/
double debye[S]; /*debye length*/
double debyetotal; /*total debye length*/
double omegap[S]; /*plasmafrequency*/
double mass[S]; /*mass of the particles*/
double charge[S]; /*charge value*/
double normalcharge[S];
double tempx[S]; /*temperature in eV, x direction*/
double tempy[S]; /*temperature in eV, y direction*/
double tempz[S];
double vthx[S]; /*thermal velocity in x dir*/
double vthy[S]; /*thermal velocity in y dir*/
double vthz[S]; 

long int npart[S]; /*number of particles*/
long int npartinit[S];
double vdriftx[S]; /*drift of the specie in x dir*/
double qm[S]; /*charge to mass ratio for each species*/
double chargeandnorm[S]; /*charge times norm factor*/
double ti2te; /*temperature ratio*/
long int lostpart[S]; /*no. of lost particles,each species*/
long int **current;/*no.of part. lost on probe,each species*/
long int **curr_av; /*current averaged, each species*/
long int **rcurr_av;
double vmean[S];/*mean velocity*/
long int allpart; /*number of all particles*/
double ratio;

/*Structure for an individual particle*/
struct particle
{
  double x;
  double y;
  double z;
  double vx;
  double vy;
  double vz;
  double kenergy;
  int llnext;
  };
typedef struct particle particle;

struct d_particle
{
  double x;
  double y;
  double z;
  int spec;
  double q;
};
typedef struct d_particle d_particle;

struct species
{
  particle part[NPART_MAX];
};
typedef struct species species;

species spec[S]; /*particles are in seperate arrays*/
long int lostlist[S][NPART_MAX];/*list of lost part.,each specie*/
long int alllost; /*total number of lost particles*/

/*particles assigned to dust grains*/
d_particle **dpart, **rdpart;
double **dpartq, **rdpartq;
/*Structure for the grid*/
/*struct grrid
{ 
  int marker;
   double Ex;
   double Ey;
  double phi;
};
typedef struct grrid grrid;
*/

double *phi;
double *phiav;
double *PE;
double *PEtotal;
double *qdens;
double *potconv[CONVTEST];

double *Fs; /*array for forces*/
double *rho; /*array for charge density and collected charge on the dust*/ 
double *pdens; /*array for particle number density*/
double *KE; /*array for the kinetic energy*/
double *vxvec, *vyvec, *vzvec;

double *rrho, *rpdens, *rKE, *rvxvec, *rvyvec, *rvzvec; /*needed for data receiving with MPI*/

int FsMAX, rhoMAX, pdensMAX, KEMAX, phiMAX, phiavMAX, PEMAX, qdensMAX, potconvMAX, PEtotalMAX; /*array dimensions*/ 
int FsEy, rhoMAXhalf, pdens_off, KE_off, FsEz, PEMAXhalf;
/*chamber*/
double Lx,Ly,Lz; /*dimensions in x and y axis in meters*/
double Gx,Gy,Gz; /*no. of grid cells in x and y directions*/

/*time and grid step factors*/
double dt, tmax;
double dx,dy,dz; /*grid spacing*/
int ngx,ngy,ngz; /*number of grids*/
double dxdt, dydt, dzdt; /* dx/dt and dy/dt AS FOR NOW NOT USED*/ 
double dtdx, dtdy, dtdz; /* dt/dx and dt/dy */
double dxdy,dV;
double dxdz;
double dxdydt,dVdt;
/*normalization factors*/
double normtime;
double normvel;
double normx;
double normpot;
double normcharge;
double normmass;
double normqm;
double normEfield;
double normdens; 
double normqdens;
double normPE;
double cellvolume;
double normPP;

/*probe*/
int probe_version;
int probex, probey; /*probe dimensions in grid points*/
double Vpr_begin, Vpr_end,Vpr_step, Vpr; /*probe potential:begining,end,step,current*/ 
double probexmin, probexmax, probeymin, probeymax;
/*parameters for probe hitting*/
int *vertp;
double *ap,*bp,*x1p,*maxx, *minx,*maxy,*miny;
int probesegments;

/***Triangle structure for triangularisation***/
struct dtriangle
{  
  int pt1;
  int pt2;
  int pt3;  
  double tcx;
  double tcy;
  double area;
  double mass;
};
typedef struct dtriangle dtriangle;

/***probe new version Oct 06********/
double **dustx, **dusty, **dustz, **dusta, **dustb, **dustxdx, **dustydy, **dustzdz, **dustbdy;
double **dustq;

/*to determine the dust crossing if dust is moving*/
double ** dustxdxold, **dustydyold, **dustzdzold;
int **dustv; //vertical marker
double *dustcx, *dustcy, *dustcz, *dustpcx, *dustpcy, *dustpcz, *dustcxdx, *dustcydx, *dustczdx;
int *ncorners;
int **lut;
int *dtype;
int *dmove;
int *dshape;
double *dradius;
double *dradiusdx;
int *dnumber;
int noofdusts;

/***for triangularisation and movement***/
double *dmass, *dustrho, *dmass_centr_y, *dmass_centr_x, *dmass_centr_z, *dmomI, **dr2v2;
double *dphifl;
dtriangle **dtrian;
int *nooftriangles;
long int *dpartlast; /*poiner to the last particle assigned to dust*/
long int *dpartmax;
/*dust movement*/
double *dustvxc, *dustvyc, *dustvzc, *dustaccx, *dustaccy, *dustaccz, *duste, *dustomega;
double **dustvx, **dustvy, **dustvz;

/*Jan idea*/
double **daa,**dbb,**dcc,**daa1y,**daa1x,**dbb1y,**dbb1x,**dcc1y,**dcc1x; 

/*for photoionization*/
double *dustworkfunct;

/*part. flux through boundaries*/
double flux[S][6];
double fluxrest[S][6];
double extrapart[S][6];
double totalflux[S];
/*boundary potential*/
double Vbound;
/*for flux calculation with Newton-Raphson method*/
double llb[S], lzet[S], lv0[S]; //lower bound, partition funct, initial velocity
double rlb[S], rzet[S], rv0[S];
double rvdriftx[S];
/*Diagnostics, files*/
FILE * fp, *fp2; //OLD
FILE *curr, *poten, *potclr, *pot2Dav; //OLD
FILE *pot2Dclr,*pot2D,*efx,*efy, *efz,*frho,*idens,*edens;
FILE *epe,*eke,*ipe,*ike, *pe, *pe_time, *eke_time, *ike_time, *epe_time, *ipe_time; //NEW2D
FILE *dustcharge, *dust_time, *convergence;
FILE *dhist;
FILE *dustshape;
FILE *eavvel, *iavvel;
FILE *evxphs, *ivxphs;
int numberofprints;
double weight;
int average;

/*Prime root generator*/
#define BUCKETSIZE 1000
double primerootbucket[BUCKETSIZE];
double primerootno;
/*Some frequently used variables*/
double sqrt_two;
double sqrt_pi;
double sqrt_twopi;
double pi;

int particlesno;
//statistics
int takecut;
FILE *dens_err;


/*Parallel part*/

int rank, numtasks, mpicheck;

/*NEW DATATYPES*/

#ifdef MPI
MPI_Datatype gridtype, oldtypes[2];
MPI_Datatype gpcolumn;
int blockcounts[12];

MPI_Aint offsets[2], extent;
MPI_Status stat;
#endif


FILE * history;
clock_t clockstart, clockend;
time_t timestart, timeending;
double timeelapsed;


FILE *probes1, *probes2, *probes3, *dustshapet;
int timerprobes;
d_particle *tmp_dpart;


//VARIABLES AND PARAMETERS FOR THE FULL MULTIGRID POISSON SOLVER
#define NGMAX 15
#define NPRE 2
#define NPOST 3
#define NCYCLES 2
int fmg_ng, fmg_nnx, fmg_nny, fmg_nnz, fmg_mingridx, fmg_mingridy, fmg_mingridz; 
double ***ires[NGMAX+1],***irho[NGMAX+1],***irhs[NGMAX+1],***iu[NGMAX+1];

int diagint, * diagint_av, * diagint_st;

//conductor new version
struct condsq
{ 
double x1,x2,x3,x4;
double y1,y2,y3,y4;
int in1,in2,in3,in4;
};
typedef struct condsq condsq;

condsq  **csq;


struct vectorst
{ 
double x;
double y;
};
typedef struct vectorst vectorst;

vectorst ** unitvec, **orthvec,  ** unitvecseg, ** orthvecseg;
int cond_present,dustmove;

struct d_rho
{
  double x;
  double y;
  double rho;
  double rho_av;
};
typedef struct d_rho d_rho;
d_rho **drho;
double *rdrho, *tordrho;
int *drholast;
int rdrholast;
double ddelta, ddelta2;
int ** ccorner;
int ** vipcorner; //virtual particle corner


//PHOTONS
int photons;
long int ph_fluxprdt;
double ph_bmin,ph_bmax,ph_xmin,ph_xmax, ph_length, ph_a;
int ph_vert; //1 up, -1 down, 0 other
double ph_angle, ph_angle_rad, ph_flux, ph_energy, ph_cosangle, ph_sinangle;
double **dustxnormv, **dustynormv;


FILE *force_chk;


long int c0[2],c1[2],c2[2],c3[2],c4[2],c5[2];

/* P-P part, creating linked list */
int *llmesh;
double *phi_nodust;
double *Fs_nodust;
int llngx, llngy, llngz;
double lldx, lldy, lldz;
long int llsize;


/*P-P part, really implemented*/
int *d_globallist;
int *d_localmax;
int **d_locallist;


FILE *testowy;


double *drag_direct_x, *drag_direct_y, *drag_direct_z, *drag_elect_x, *drag_elect_y, *drag_elect_z;
double *drot_z_x1, *drot_z_x2, *drot_y_x1, *drot_y_x2, *drot_z_y1, *drot_z_y2, *drot_y_z1, *drot_y_z2, *drot_x_y1, *drot_x_y2, *drot_x_z1, *drot_x_z2;
double *elrot_z_x1, *elrot_z_x2, *elrot_y_x1, *elrot_y_x2, *elrot_z_y1, *elrot_z_y2, *elrot_y_z1, *elrot_y_z2, *elrot_x_y1, *elrot_x_y2, *elrot_x_z1, *elrot_x_z2;

