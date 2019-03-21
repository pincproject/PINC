/**
 * @file		aux.test.c
 * @brief		Unit tests for aux.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
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

	b[4] = 3;
	bl[4] = 3;

	utAssert(!aiEq(a,b,5),"aiEq is broken");
	utAssert(!alEq(al,bl,5),"alEq is broken");

	return 0;
}

// All tests for aux.c is contained in this function
void testAux(){
	utRun(&testAiProd);
	utRun(&testAEq);
}
