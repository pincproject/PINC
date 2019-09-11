"""
/**
* @file                GenGrid.py
* @author              Jan Deca <jandeca@gmail.com>
* @copyright           None
* @brief               File containing all tools to generate a grid with embedded objects for PINC
* @date                03.03.2016
*
* 
*/
"""

import sys as sys
import os.path
import numpy as np
import vtk as vtk
import h5py as h5py
from math import *

# Check the sanity of the input (to be extended as needed)
def checkSanity(infile,outfile):

    # Check whether we are not overwriting an input VTK file
    for i in range(0,len(infile)):
        if (infile[i]==outfile[0] or infile[i]==outfile[1]):
            print " The output file will overwrite one of the input vtk files. Aborting.\n"
            sys.exit()
        #else:
        #    

    # Check if output file already exists
    for i in range(0,len(outfile)-1):
        if outfile[i].endswith('.vtk') and os.path.isfile(outfile[i]):
            print " Output VTK file already exists. Aborting.\n"
            sys.exit()
        if outfile[i].endswith('.h5') and os.path.isfile(outfile[i]):
            print " Output H5 file already exists. Aborting.\n"
            sys.exit()

    # Continue happily
    print " Your input seems sane, continuing... (Honestly, didn't really check, to be implemented.)\n"

    return
            
# Read an unstructured grid VTK file
def readUnstructuredVTK(infile):
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(infile+".vtk")
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    
    polydata = reader.GetOutput()
    NumCells = polydata.GetNumberOfCells()
    np_pts = np.zeros((NumCells,3,3))
    
    for i in range(0,NumCells):
        pts = polydata.GetCell(i).GetPoints()
        np_pts[i]=np.array([pts.GetPoint(j) for j in range(pts.GetNumberOfPoints())])
        
    return np_pts

# Transform a triangle.
# Not fixing the size of the dummy will probably allow to merge with transformSeeds
def transformObject(obj_ori, transfo):
    dummy = np.empty((3,3))
    a = cos(transfo[3]*3.14159265359/180)
    b = sin(transfo[3]*3.14159265359/180)
    c = cos(transfo[4]*3.14159265359/180)
    d = sin(transfo[4]*3.14159265359/180)
    e = cos(transfo[5]*3.14159265359/180)
    f = sin(transfo[5]*3.14159265359/180)
    for t in range(len(obj_ori[:,0,0])):
        #Scale
        obj_ori[t,:,0] *= transfo[6]
        obj_ori[t,:,1] *= transfo[7]
        obj_ori[t,:,2] *= transfo[8]        
        #Rotate
        dummy[:,0] =  c*e*obj_ori[t,:,0] + (a*f+b*d*e)*obj_ori[t,:,1] + (b*f-a*d*e)*obj_ori[t,:,2]
        dummy[:,1] = -c*f*obj_ori[t,:,0] + (a*e-b*d*f)*obj_ori[t,:,1] + (b*e+a*d*f)*obj_ori[t,:,2]
        dummy[:,2] =  d*obj_ori[t,:,0] - b*c*obj_ori[t,:,1] + a*c*obj_ori[t,:,2]
        obj_ori[t,:,0] = dummy[:,0]
        obj_ori[t,:,1] = dummy[:,1]
        obj_ori[t,:,2] = dummy[:,2]
        #Translate
        obj_ori[t,:,0]+=transfo[0]
        obj_ori[t,:,1]+=transfo[1]
        obj_ori[t,:,2]+=transfo[2]

    return obj_ori

# Transform a seed.
# Probably merge this function with transformObject
def transformSeeds(content, transfo):
    a = cos(transfo[3]*3.14159265359/180)
    b = sin(transfo[3]*3.14159265359/180)
    c = cos(transfo[4]*3.14159265359/180)
    d = sin(transfo[4]*3.14159265359/180)
    e = cos(transfo[5]*3.14159265359/180)
    f = sin(transfo[5]*3.14159265359/180)

    for p in range(content[0]):
        seed = np.array(content[p+1], dtype=np.float)
        # Scale
        seed[0] *= transfo[6]
        seed[1] *= transfo[7]
        seed[2] *= transfo[8]
        # Rotate
        dummy0 =  c*e*seed[0] + (a*f+b*d*e)*seed[1] + (b*f-a*d*e)*seed[2]
        dummy1 = -c*f*seed[0] + (a*e-b*d*f)*seed[1] + (b*e+a*d*f)*seed[2]
        dummy2 =  d*seed[0] - b*c*seed[1] + a*c*seed[2]
        seed[0] = dummy0
        seed[1] = dummy1
        seed[2] = dummy2
        # Translate
        seed[0]+=transfo[0]
        seed[1]+=transfo[1]
        seed[2]+=transfo[2]
        # Place the seed back in the list
        content[p+1] = tuple(seed)
        
    return content

# Find the circumfering voxels for each object
def findCircVoxels(grid, gridpar, obj_trans, ident):
    print "  3.1 Circumference."
    xmin = gridpar[0]
    ymin = gridpar[2]
    zmin = gridpar[4]
    nnx = gridpar[6]
    nny = gridpar[7]
    nnz = gridpar[8]
    dx = (gridpar[1] - gridpar[0])/float(nnx)
    dy = (gridpar[3] - gridpar[2])/float(nny)
    dz = (gridpar[5] - gridpar[4])/float(nnz)
    
    for t in range(len(obj_trans)):
        # Find the bounding box of each triangle
        minx = int(round((obj_trans[t,:,0].min()-xmin)/dx))
        if (minx<0): minx=0
        if (minx>=nnx): minx=nnx-1
        maxx = int(round((obj_trans[t,:,0].max()-xmin)/dx))
        if (maxx<0): maxx=0
        if (maxx>=nnx): maxx=nnx-1
        miny = int(round((obj_trans[t,:,1].min()-ymin)/dy))
        if (miny<0): miny=0
        if (miny>=nny): miny=nny-1
        maxy = int(round((obj_trans[t,:,1].max()-ymin)/dy))
        if (maxy<0): maxy=0
        if (maxy>=nny): maxy=nny-1
        minz = int(round((obj_trans[t,:,2].min()-zmin)/dz))
        if (minz<0): minz=0
        if (minz>=nnz): minz=nnz-1
        maxz = int(round((obj_trans[t,:,2].max()-zmin)/dz))
        if (maxz<0): maxz=0
        if (maxz>=nnz): maxz=nnz-1

        # Tag node on the grid if it belongs to the circumference
        for i in range(minx,maxx+1):
            for j in range(miny,maxy+1):
                for k in range(minz,maxz+1):
                    P = np.array([i*dx+xmin,j*dy+ymin,k*dz+zmin])
                    dist, pp0 = pointTriangleDistance(obj_trans[t,:,:],P)
                    #if (dist<dx): # NOTE, this gives you trouble if not using a uniform grid !!!
                    #    grid[i,j,k] = ident
                    if (abs(P[0]-pp0[0])<dx and abs(P[1]-pp0[1])<dy and abs(P[2]-pp0[2])<dz): # This is better, possibly :-)
                        grid[i,j,k] = ident
                    if (i<=0 or i>=nnx-1 or j<=0 or j>=nny-1 or k<=0 or k>=nnz-1):
                        grid[i,j,k] = ident
                        
    return grid

# 3D Floodfill to find all voxels within the circumference 
def floodFill(grid, gridpar, content):
    print "  3.2 Internal volume."
    xmin = gridpar[0]
    ymin = gridpar[2]
    zmin = gridpar[4]
    nnx = gridpar[6]
    nny = gridpar[7]
    nnz = gridpar[8]
    dx = (gridpar[1] - gridpar[0])/float(nnx)
    dy = (gridpar[3] - gridpar[2])/float(nny)
    dz = (gridpar[5] - gridpar[4])/float(nnz)

    stacksize = nnx*nny*nnz # Stacksize for deep recursion. Probably too big, but safe.
    stack = np.zeros((stacksize,3))
    
    for i in range(content[0]):
        stack[i,0] = int(content[i+1][0]/dx-xmin/dx)
        stack[i,1] = int(content[i+1][1]/dy-ymin/dy)
        stack[i,2] = int(content[i+1][2]/dz-zmin/dz)

    nnn = content[0]-1
    while nnn>-1:
        stack,nnn = floodfill3Dstack(stack,grid,nnn, content[-1])
        
    return grid

# Write the requested output files
def writeOutput(grid, gridpar, outfile):
    for i in range(0,len(outfile)-1):
        if outfile[i].endswith('.vtk'):
            print "     ", outfile[i]
            writeLegacyVTK(grid, gridpar, outfile[i], outfile[-1])
        elif outfile[i].endswith('.h5'):
            print "     ", outfile[i]
            writeH5(grid, gridpar, outfile[i], outfile[-1])

# Write a legacy VTK (version 1.0)
def writeLegacyVTK(grid, gridpar, outfile, comment):
    nnx = gridpar[6]
    nny = gridpar[7]
    nnz = gridpar[8]

    # Convert grid to 1D
    grid1D = np.empty(nnx*nny*nnz)
    for k in range(0,nnz):
        for j in range(0,nny):
            for i in range(0,nnx):
                grid1D[k*nnx*nny+j*nnx+i] = grid[i][j][k]

    # Write the file in ASCII
    file = open(outfile,'w')
    file.write("# vtk DataFile Version 1.0\n")
    file.write(comment + "\n")
    file.write("ASCII\n")
    file.write("DATASET STRUCTURED_POINTS\n")
    file.write("DIMENSIONS " + str(nnx) + " " + str(gridpar[7]) + " " + str(gridpar[8]) + "\n")
    file.write("ORIGIN " + str(gridpar[0]) + " " + str(gridpar[2]) + " " + str(gridpar[4]) + "\n")
    file.write("SPACING " + str((gridpar[1] - gridpar[0])/float(nnx)) + " " + str((gridpar[3] - gridpar[2])/float(nny)) + " " + str((gridpar[5] - gridpar[4])/float(nnz)) + "\n")
    file.write("POINT_DATA " + str(nnx*nny*nnz) + "\n")
    file.write("SCALARS grid float\n")
    file.write("LOOKUP_TABLE default\n")
    for i in range(len(grid1D)):
        file.write(str(grid1D[i]) + "\n")
    file.close()

    return

# Write a H5
def writeH5(grid, gridpar, outfile,comment):
    grid4D = np.expand_dims(grid, axis=3)
    fileObj = h5py.File(outfile,'x')
    dset = fileObj.create_dataset("Object",(gridpar[6],gridpar[7],gridpar[8],1),dtype='i1', data=grid4D)
    fileObj.close()

    return
    
# Algorithm to find all voxels with the circumference of the object/geometry using a predefined stack to avoid deep recursion issues
def floodfill3Dstack(stack, grid,nnn,ident):
    if grid[int(stack[nnn,0]),int(stack[nnn,1]),int(stack[nnn,2])]==ident:
        nnn -=1
        return stack,nnn
    else:
        grid[int(stack[nnn,0]),int(stack[nnn,1]),int(stack[nnn,2])]=ident
        #down
        stack[nnn+1,0] = stack[nnn,0]-1
        stack[nnn+1,1] = stack[nnn,1]
        stack[nnn+1,2] = stack[nnn,2]
        #right
        stack[nnn+2,0] = stack[nnn,0]
        stack[nnn+2,1] = stack[nnn,1]+1
        stack[nnn+2,2] = stack[nnn,2]
        #left
        stack[nnn+3,0] = stack[nnn,0]
        stack[nnn+3,1] = stack[nnn,1]-1
        stack[nnn+3,2] = stack[nnn,2]
        #back
        stack[nnn+4,0] = stack[nnn,0]
        stack[nnn+4,1] = stack[nnn,1]
        stack[nnn+4,2] = stack[nnn,2]+1
        #front
        stack[nnn+5,0] = stack[nnn,0]
        stack[nnn+5,1] = stack[nnn,1]
        stack[nnn+5,2] = stack[nnn,2]-1
        #up
        stack[nnn,0] = stack[nnn,0]+1
        stack[nnn,1] = stack[nnn,1]
        stack[nnn,2] = stack[nnn,2]
        
        nnn +=5
        
        return stack,nnn

#Finds the bounding box for the grid, works for one object currently
def boundingBox(grid):

    out = np.zeros(grid.shape)
    x = np.any(grid, axis=(1, 2))
    y = np.any(grid, axis=(0, 2))
    z = np.any(grid, axis=(0, 1))

    xmin, xmax = np.where(x)[0][[0, -1]]
    ymin, ymax = np.where(y)[0][[0, -1]]
    zmin, zmax = np.where(z)[0][[0, -1]]

    out[xmin-1:xmax+1,ymin-1:ymax+1,zmin-1:zmax+1] = 2

    return out

# Algorithm to find the circumfering voxels of the object/geometry.
def pointTriangleDistance(TRI, P):
    # function [dist,PP0] = pointTriangleDistance(TRI,P)
    # calculate distance between a point and a triangle in 3D
    # SYNTAX
    #   dist = pointTriangleDistance(TRI,P)
    #   [dist,PP0] = pointTriangleDistance(TRI,P)
    #
    # DESCRIPTION
    #   Calculate the distance of a given point P from a triangle TRI.
    #   Point P is a row vector of the form 1x3. The triangle is a matrix
    #   formed by three rows of points TRI = [P1;P2;P3] each of size 1x3.
    #   dist = pointTriangleDistance(TRI,P) returns the distance of the point P
    #   to the triangle TRI.
    #   [dist,PP0] = pointTriangleDistance(TRI,P) additionally returns the
    #   closest point PP0 to P on the triangle TRI.
    #
    # Author: Gwolyn Fischer
    # Release date: 09/02/02
    #
    # The algorithm is based on
    # "David Eberly, 'Distance Between Point and Triangle in 3D',
    # Geometric Tools, LLC, (1999)"
    # http:\\www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
    #
    #        ^t
    #  \     |
    #   \reg2|
    #    \   |
    #     \  |
    #      \ |
    #       \|
    #        *P2
    #        |\
    #        | \
    #  reg3  |  \ reg1
    #        |   \
    #        |reg0\
    #        |     \
    #        |      \ P1
    # -------*-------*------->s
    #        |P0      \
    #  reg4  | reg5    \ reg6
    
    # rewrite triangle in normal form
    B = TRI[0,:]
    E0 = TRI[1,:] - B
    # E0 = E0/sqrt(sum(E0.^2)); %normalize vector
    E1 = TRI[2,:] - B
    # E1 = E1/sqrt(sum(E1.^2)); %normalize vector
    D = B - P
    a = np.dot(E0,E0)
    b = np.dot(E0,E1)
    c = np.dot(E1,E1)
    d = np.dot(E0,D)
    e = np.dot(E1,D)
    f = np.dot(D,D)
    
    #print "{0} {1} {2} ".format(B,E1,E0)
    det = a*c-b*b
    s = b*e-c*d
    t = b*d-a*e
    
    # Terible tree of conditionals to determine in which region of the diagram
    # shown above the projection of the point into the triangle-plane lies.
    if (s+t)<=det:
        if s<0.0:
            if t<0.0:
                # region4
                if d<0:
                    t = 0.0
                    if -d>=a:
                        s = 1.0
                        sqrdistance = a+2.0*d+f
                    else:
                        s = -d/a
                        sqrdistance = d*s+f
                else:
                    s = 0.0
                    if e>=0.0:
                        t = 0.0
                        sqrdistance = f
                    else:
                        if -e>=c:
                            t = 1.0
                            sqrdistance = c+2.0*e+f
                        else:
                            t = -e/c
                            sqrdistance = e*t+f
                            
                            # of region 4
            else:
                # region 3
                s = 0
                if e>=0:
                    t=0
                    sqrdistance = f
                else:
                    if -e>=c:
                        t = 1
                        sqrdistance = c+2.0*e+f
                    else:
                        t = -e/c
                        sqrdistance = e * t + f
                        # of region 3
        else:
            if t<0:
                # region 5
                t = 0
                if d>=0:
                    s = 0
                    sqrdistance = f
                else:
                    if -d>=a:
                        s = 1
                        sqrdistance = a+2.0*d+f;  # GF 20101013 fixed typo d*s ->2*d
                    else:
                        s = -d/a
                        sqrdistance = d*s+f
            else:
                # region 0
                invDet = 1.0/det
                s = s*invDet
                t = t*invDet
                sqrdistance = s*(a*s+b*t+2.0*d)+t*(b*s+c*t+2.0*e)+f
    else:
        if s<0.0:
            # region 2
            tmp0 = b+d
            tmp1 = c+e
            if tmp1>tmp0:  # minimum on edge s+t=1
                numer = tmp1-tmp0
                denom = a-2.0*b+c
                if numer>=denom:
                    s = 1.0
                    t = 0.0
                    sqrdistance = a+2.0*d+f;  # GF 20101014 fixed typo 2*b -> 2*d
                else:
                    s = numer/denom
                    t = 1-s
                    sqrdistance = s*(a*s+b*t+2*d)+t*(b*s+c*t+2*e)+f
                    
            else:  # minimum on edge s=0
                s = 0.0
                if tmp1<=0.0:
                    t = 1
                    sqrdistance = c+2.0*e+f
                else:
                    if e>=0.0:
                        t = 0.0
                        sqrdistance = f
                    else:
                        t = -e/c
                        sqrdistance = e*t+f
                        # of region 2
        else:
            if t < 0.0:
                # region6
                tmp0 = b+e
                tmp1 = a+d
                if tmp1>tmp0:
                    numer = tmp1-tmp0
                    denom = a-2.0*b+c
                    if numer>=denom:
                        t = 1.0
                        s = 0
                        sqrdistance = c+2.0*e+f
                    else:
                        t = numer/denom
                        s = 1-t
                        sqrdistance = s*(a*s+b*t+2.0*d)+t*(b*s+c*t+2.0*e)+f
                        
                else:
                    t = 0.0
                    if tmp1<=0.0:
                        s = 1
                        sqrdistance = a+2.0*d+f
                    else:
                        if d>=0.0:
                            s = 0.0
                            sqrdistance = f
                        else:
                            s = -d/a
                            sqrdistance = d*s+f
            else:
                # region 1
                numer = c+e-b-d
                if numer<=0:
                    s = 0.0
                    t = 1.0
                    sqrdistance = c+2.0*e+f
                else:
                    denom = a-2.0*b+c
                    if numer>=denom:
                        s = 1.0
                        t = 0.0
                        sqrdistance = a+2.0*d+f
                    else:
                        s = numer/denom
                        t = 1-s
                        sqrdistance = s*(a*s+b*t+2.0*d)+t*(b*s+c*t+2.0*e)+f
                        
    # account for numerical round-off error
    if sqrdistance < 0:
        sqrdistance = 0
        
    dist = sqrt(sqrdistance)
    
    PP0 = B+s*E0+t*E1
    return dist, PP0
