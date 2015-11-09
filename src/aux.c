/**
 * @file		aux.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Small auxiliary functions
 * @date		29.10.15
 *
 * Small auxiliary functions.
 */

#include "pinc.h"

/******************************************************************************
 * ARRAY FUNCTIONS
 *****************************************************************************/

int *intArrMul(const int *a, const int *b, int nElements){

	int *result = malloc(nElements*sizeof(int));

	for(int i=0;i<nElements;i++){
		result[i]=a[i]*b[i];
	}

	return result;

}

int *intArrCumProd(const int *a, int nElements){

	int *result = malloc((nElements+1)*sizeof(int));
	result[0]=1;

	for(int i=1;i<nElements+1;i++){
		result[i]=result[i-1]*a[i-1];
	}

	return result;

}

int intArrProd(const int *a, int nElements){

	int result = 1;
	for(int i=0;i<nElements;i++){
		result *= a[i];
	}

	return result;
}

/******************************************************************************
 * TIMING FUNCTIONS
 *****************************************************************************/
