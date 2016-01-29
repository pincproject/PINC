/**
 * @file		aux.test.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 * @copyright	University of Oslo, Norway
 * @brief		Unit tests for aux.c
 * @date		16.12.15
 */

#include "pinc.h"
#include "test.h"

static int testAiProd(){
	int arr[] = {2,34,9,1,6,6};
	int prod = aiProd(arr,6);
	utAssert(prod==22032,"aiProd doesn't work");
	return 0;
}

static int testAEq(){
	int a[] = {2,3,4,5,6};
	int b[] = {2,3,4,5,6};
	long int al[] = {2,3,4,5,6};
	long int bl[] = {2,3,4,5,6};

	utAssert(aiEq(a,b,5),"aiEq is broken");
	utAssert(alEq(al,bl,5),"alEq is broken");
	utAssert(aEq((char*)a,(char*)b,5*sizeof(*a)),"aEq is broken 1");
	utAssert(aEq((char*)al,(char*)bl,5*sizeof(*al)),"aEq is broken 2");

	b[4] = 3;
	bl[4] = 3;

	utAssert(!aiEq(a,b,5),"aiEq is broken");
	utAssert(!alEq(al,bl,5),"alEq is broken");
	utAssert(!aEq((char*)a,(char*)b,5*sizeof(*a)),"aEq is broken 3");
	utAssert(!aEq((char*)al,(char*)bl,5*sizeof(*al)),"aEq is broken 4");


	return 0;
}

// All tests for aux.c is contained in this function
void testAux(){
	utRun(&testAiProd);
	utRun(&testAEq);
}
