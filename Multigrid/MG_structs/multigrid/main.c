#include <stdio.h>
#include <stdlib.h>
#include "mgStructs.h"

int main(){


    struct Grid grid;
    struct Multigrid multigrid;

    int nDim = 2;
    int nValues = 1;
    int nNodes[2] = {4, 4};
    int nLevels = 2;

    gridInit(&grid, nNodes, nDim, nValues);

    multigridInit(&multigrid, &grid, nLevels);

    printf("Grid nodes in grid[0]: %d \n", multigrid.grids[0].nNodes[0]);
    printf("Grid nodes in grid[1]: %d \n", multigrid.grids[1].nNodes[0]);

    return 0;
}

