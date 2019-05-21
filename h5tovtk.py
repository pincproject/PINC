import h5py
import numpy as np
import sys

filePhi = h5py.File(sys.argv[1],'r')


for ts in range(int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])):
    phi = filePhi['/n=' + str(ts) + '.0']

    phi = phi[:,:,:,0]

#    print "Converting to 1D"
    H1D = np.empty(32*32*32)


    for k in range(0,32):
        for j in range(0,32):
            for i in range(0,32):
                H1D[k*32*32+j*32+i] = phi[i][j][k]

 #   print "Writing vtk"
    file = open(str(sys.argv[1])+'_'+ str(ts) + '.vtk','w')
    file.write("# vtk DataFile Version 1.0"+"\n")
    file.write("testrho"+"\n")
    file.write("ASCII"+"\n")
    file.write("DATASET STRUCTURED_POINTS"+"\n")
    file.write("DIMENSIONS 32 32 32"+"\n")
    file.write("ORIGIN 0 0 0"+"\n")
    file.write("SPACING 0.0104167 0.0104167 0.0104167"+"\n")
    file.write("POINT_DATA 32768"+"\n")
    file.write("SCALARS test float"+"\n")
    file.write("LOOKUP_TABLE default"+"\n")
    for i in range(len(H1D)):
        file.write(str(H1D[i])+"\n")


