#!/usr/bin/env python3
import h5py
import numpy as np
import pyvista as pv

# 1. Read the HDF5 file
with h5py.File('../..data/object.grid.h5', 'r') as f:
    # dataset is stored in z, y, x order
    data = f['Object'][:]  
zdim, ydim, xdim = data.shape

# 2. Build a UniformGrid
#    VTK structured‐points (i.e. uniform grid) has NCells=(xdim,ydim,zdim)
#    and NPoints=(xdim+1, ydim+1, zdim+1)
grid = pv.UniformGrid()
grid.dimensions = (xdim+1, ydim+1, zdim+1)
grid.origin     = (0.0, 0.0, 0.0)   # match your coordinate origin
grid.spacing    = (1.0, 1.0, 1.0)   # cell size = 1 in each direction

# 3. Attach your data as CELL scalars
#    Flatten in C‐order so it matches how VTK expects the cells laid out
grid.cell_data["Object"] = data.flatten(order="C")

# 4. Write out to VTK
grid.save('object.vtk')
print("Wrote object.vtk — load this in ParaView to visualize your sphere.")