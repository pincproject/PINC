/*DUSTY PLASMA program*/
/*library of all functions*/
/*Author: Wojciech Jacek Miloch*/
/*University of Oslo, Norway*/
/*June 2006*/
/*input.c*/
void convert(void);
void readdata(int arc, char *arv[]);
/*flux.c*/
double erfcc(double x);
double cumfprim(double v0, int i);
double cumf(double v0, int i);
double zet1(int i);
void init_newpart(void);
void calculate_flux(void);
/*diagn.c*/
void diagn_open(void);
void diagn_close(void);
//2D diagnostics
void printgrid(int t);
void printdustshape(void);
void printdustshapetime(int t);
void printscale(FILE *fpointer);
void print_avpvel(FILE *fpointer, int k, int t, double weight);

void printqdensity(FILE *fpointer, int t, double weight);
void printdensity(FILE *fpointer, int k, int t, double weight);
void printpotential(FILE *fpointer, int t, double weight);
void printavpotential(FILE *fpointer, int t, double weight);
void printefield(FILE *fpointer, int help, int t, double wieght);
void printKEall(int t);
void printKE(FILE *fpointer, FILE *fpointer2, int specie, int t, double weight);
void printPE(FILE *fpointer, FILE *fpointer2, int specie, int t, double weight);
void printPEtotal(FILE *fpointer, FILE *fpointer2, int t, double weight);
void printdustcharge(FILE *fpointer, int t, double weight);
void printdustchargetime(FILE *fpointer, int t, double weight);
void print_current(int tid);

void printdth(FILE *fpointer, int t);
void printconvpot(FILE *fpointer, int t, int step);
//OLD DIAGN
void printpotcut(FILE *fpointer);
double printall(FILE *fpoint1, FILE *fpoint2, int collect);

void pot_probes_init(void);
void pot_probes(int t);
/*generate.c*/
void gen_bgnd(void);
void maxw_dist(int i, double vthx, double vthy, double driftx, double drifty);
void newparticles(int timestep);
void init_primeroot(double seed);
double primeroot(void);
int initpartcheck(double px, double py, double pz, double delta);
int initpartcheck_restart(int dno, double px, double py, double pz, double delta);
/*accel.c*/
void accel(float factor);
void move(int t);
void create_linkedlist(void);
/*grid.c*/
void memorygrid(void);
void memorygridfree(void);
void weighting1(void);
void cleargrid(void);
void cleargrid2(void);
void gen_boundaries(void);
void gen_probe(int version);
void gen_dust3D(int arc, char *arv[]);
void markgriddust(void);
void checkcolcrossing(int i);
int checkpointcrossing(int i, int j);
void findabv(int i);
void new_probe_potential(double probepotential);
void normalize(void);
void memorydust2_3D(int j, int nc);
void memorydust1_3D(int no);
void create_currentarrays(void);

 
/*gauss.c*/
void gauss_seidel(int nx, int ny, double tolerance);
void electric_field(void);
/*dustg.c*/
double finddustvolume(int arc, char *arv[]);
double dustarea(int arc, char *arv[]); 
void d_polygon(int arc, char *arv[]);
void calculate_staticparameters(int arc, char *arv[]);
void d_centreofmass_and_momI(void);
void d_move(int t);
void memoryduststatic(int no);
void memorydpart(int no, int dmax);
void weightingdust1(int ko);
void chargeoncond(int i);
void virtpart(void);
void ortnormvec(void);
void condsquares(void);
void average_current(void);
void findnewpotentials(double interval, int collect, FILE *fpoint1);

inline int signof(int a);
inline double smaller_same_sign(double a, double b);
void printdragforce(int timestep);
void drag_force_electric(void);
void drag_force_direct(double partvxnew, double partvynew, double partvznew, int particlespecie, int dno, double partxhit, double partyhit, double partzhit);


/*program for polygon triangulation*/

/*full multigrid method fgm.c*/
void mglin_destroy(void);
void mglin_init(int nx, int ny, int nz);
void mglin(double *u, int ncycle);


/*restart facility restart.c*/
void dump(long int t);
long int prog_restart(void);
void shift_while_restarting(int dno, double x, double y, double z);

void checkcond(void);


/*photons*/
void photonflux(void);
void photoelectriceffect(void);

/*shortcuts*/
double *dvecmem(long nl,long nh);
int *ivecmem(long nl,long nh);
void free_dvecmem(double *v, long nl, long nh);
void free_ivecmem(int *v, long nl, long nh);

//in fmg 
void nnrerror(char error_text[]);

FILE * my_file_open(const char * filename, const char * aarg);
inline int ix(int off, int i, int j, int k);

/*spherical*/
void points_on_sphere(int dustnumber, int numberofpoints);
