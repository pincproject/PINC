#include "halo.h"

double *getSlice(double *nextGhost, const double **valp, const long int *mul, const int *points){
	for(int j=0;j<*mul;j++) *(nextGhost++) = *((*valp)++);
	if(*mul!=1) for(int k=0;k<*points;k++) nextGhost = getSlice(nextGhost,valp,mul-1,points-1);
	else *valp += *points;
	for(int j=0;j<*mul;j++) *(nextGhost++) = *((*valp)++);
	return nextGhost;
}

inline void getHalo(double *halo, const double *val, const long  int *nGPointsProd, const int *nTGPoints,int nDims){
	getSlice(halo,&val,&nGPointsProd[nDims-1],&nTGPoints[nDims-1]);
}
