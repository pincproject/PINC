#include "core.h"
#include "test.h"

static int testOSolFacingSurfaceNodes(dictionary *ini){
    
    Units *units=uAlloc(ini);
	uNormalize(ini, units);
    MpiInfo *mpiInfo = gAllocMpi(ini);
	Population *pop = pAlloc(ini,mpiInfo);
    Object *obj = oAlloc(ini,mpiInfo,units);              // for capMatrix - objects

    double expected[] = {10,11,12,13};

    utAssert(alEq(obj->exposedNodes, expected,4), "Solar facing nodes are incorrect");


}

void testObject(){
    utRun(&testOSolFacingSurfaceNodes);
}

