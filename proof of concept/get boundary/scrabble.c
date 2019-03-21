/*
 * METHOD A: Not general dimensionality
 */

for(int j=0;j<Ngx-1;j++){
	p = j;
	store(p);
} // p=4
for(int k=1;k<Ngy-2;i++){
	p++;		// p = 5, 10
	store(p);
	p+=Ngx-1;	// p = 9, 14
	store(p);
}
for(int j=0;j<Ngx-1;j++){
	p++;
	store(p);	// p = 15,16,17,18,19
}

/*
 * METHOD B: Using pre-computed indices
 */

nGI = ...
 int *ghostIndices = malloc(nGI*sizeof(int));
 int *nextGhost = ghostIndices;
 for(int p=0; p<nTot; p++){

	 int temp = p;

	 int i = temp % Nx;
	 temp /= Nx;
	 int j = temp % Ny;
	 temp /= Ny;
	 int k = temp % Nz;

	 if(i==0 || i==Nx-1 || j==0 || j==Ny-1 || k==0 || k==Nz-1){
		 *(nextGhost++) = p;
	 }

 }

 double *getHalo(...){

	 double *halo = gridQuantity->halo;
	 double *val = gridQuantity->val;
	 double *ghostIndices = gridQuantity->ghostIndices;
	 int nGI = grid->nGI; // Number of ghost indices

	 for(int p=0;p<nGI;p++){
		 halo[p] = val[ghostIndices[p]];
	 }

 }

 /*
  * METHOD C: Generalized B
  */

int *ghostIndices = malloc(...);

int *getGhostIndices(...){
	for(int p=0; p<nGPointsProd[nDims]; p++){

		int temp = p;
		for(int d=0; d<nDims; d++){

			// Find index d and prepare for next iteration
			int index = temp % nGPoints[d];
			temp /= nGPoints[d];

			if(index==0 || index==nGPoints[d]-1){
				*(nextGhost++) = index;
				break;	// avoids double entries of corners etc.
			}
		}
	 }
 }

 /*
  * METHOD D:
  */

p=-1;
for(int j=0;j<Ngx;j++) store(++p);	// After this: p=Ngx-1
for(int k=0;k<Ngy-1;k++){
	store(++p);
	store(p+=Ngx-1);
} // After this: p=(Ngy-1)*Ngx-1
for(int j=0;j<Ngx;j++) store(++p);

/*
 * METHOD E:
 */

p=0;
for(int j=0;j<nGPointsProd[1];j++) store(p+=nGPointsProd[0]);
for(int k=1;k<=nGPoints[1]-2;k++){
	store(p+=nGPointsProd[0]);
	store(p+=nGPoints[0]-1);
}
for(int j=0;j<nGPointsProd[1];j++) store(p+=nGPointsProd[0]);

l=0;

p=0;
int *mul = &grid->nGPointsProd[nDims-1];
int *points = &grid->nGPoints[nDims-1];
int *val = gridQuantity->val;

storeLayer(prod,nPoints,ghost,val);

function storeLayer(const int *prod,const int *nPoints,double *ghost, double *val){
	for(int j=0;j<prod[0];j++) *(nextGhost++) = *(val+=prod[1]);
	for(int k=0;k<prod[1];k++){
		storeLayer(prod+1,nPoints+1,nextGhost,val);
	}
	for(int j=0;j<prod[0];j++) *(nextGhost++) = val[ p+=prod[1] ];
}

getGhosts(double *nextGhost, const double *val, const int *mul, const int *points){
	for(int j=0;j<*mul;j++) *(nextGhost++) = *(++val);
	if(*mul!=1)
		for(int k=0;k<*points;k++){
			getGhosts(nextGhost,val,mul--,points--);
		}
	}
	val += *(mul) * (*(points)-1);	
	for(int j=0;j<*mul;j++) *(nextGhost++) = *(++val);
}
