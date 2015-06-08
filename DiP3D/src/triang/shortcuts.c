/*DiP3D */
/*Functions to make the code more readible */

#include "const.h"
#include "./fmg/nrutil.h"
#include <stddef.h>
#include <stdlib.h>

/*a short procedure for the file opening*/
FILE * my_file_open(const char * filename, const char * aarg){

   FILE *fpointer;   
   if((fpointer=fopen(filename, aarg))==NULL)
	{
    	printf("I can not open the file %s\n", filename);
		printf("Sorry, but I am exiting program now\n"); 
	    exit(1);
	}
	return fpointer;
}


//BASED on NR in C.
#define NR_ENDD 1
#define FREE_ARGG char*

double *dvecmem(long nl,long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_ENDD)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvecmem()\n");
	return v-nl+NR_ENDD;
}

void free_dvecmem(double *v, long nl, long nh)
/* free a double vector allocated with dvecmem() */
{
	free((FREE_ARGG) (v+nl-NR_ENDD));
}

int *ivecmem(long nl,long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_ENDD)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivecmem()\n");
	return v-nl+NR_ENDD;
}

void free_ivecmem(int *v, long nl, long nh)
/* free a double vector allocated with dvecmem() */
{
	free((FREE_ARGG) (v+nl-NR_ENDD));
}


/****************finding an index ************************/
/*3D*/
inline int ix(int off, int i, int j, int k)
{
  return off+i*ngy*ngz+j*ngz+k;
}
