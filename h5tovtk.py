import h5py
import numpy as np
import sys

filePhi = h5py.File(sys.argv[1],'r')

nxn = 32
nyn = 32
nzn = 32


for ts in range(int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])):
    phi = filePhi['/n=' + str(ts) + '.0']

    phi = phi[:,:,:,0]

#    print "Converting to 1D"
    H1D = np.empty(nxn*nyn*nzn)


    for k in range(0,nzn):
        for j in range(0,nyn):
            for i in range(0,nxn):
                H1D[k*nxn*nyn+j*nxn+i] = phi[i][j][k]

 #   print "Writing vtk"
    file = open(str(sys.argv[1])+'_'+ str(ts) + '.vtk','w')
    file.write("# vtk DataFile Version 1.0"+"\n")
    file.write("testrho"+"\n")
    file.write("ASCII"+"\n")
    file.write("DATASET STRUCTURED_POINTS"+"\n")
    file.write("DIMENSIONS " + str(nxn) + " " + str(nyn) +" " + str(nzn) +"\n")
    file.write("ORIGIN 0 0 0"+"\n")
    file.write("SPACING 0.0104167 0.0104167 0.0104167"+"\n")
    file.write("POINT_DATA "+ str(nxn*nyn*nzn)+"\n")
    file.write("SCALARS test float"+"\n")
    file.write("LOOKUP_TABLE default"+"\n")
    for i in range(len(H1D)):
        file.write(str(H1D[i])+"\n")


