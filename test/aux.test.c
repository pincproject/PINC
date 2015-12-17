/**
 * @file		aux.test.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 * @copyright	University of Oslo, Norway
 * @brief		Unit tests for aux.c
 * @date		16.12.15
 */

#include "pinc.h"
#include "test.h"

int testIntArrProd(){
	int arr[] = {2,34,9,1,6,6};
	int prod = intArrProd(arr,6);
	utAssert(prod==22032,"intArrProd doesn't work");
	return 0;
}

// All tests for aux.c is contained in this function
void testAux(){
	utRun(&testIntArrProd);
}
