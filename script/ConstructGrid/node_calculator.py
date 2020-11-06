import numpy as np
import csv

nodes = np.loadtxt('nodes.txt',dtype=int)
nodes = list(nodes)

#pos = []
pos2 = []
#sizeprod1 = [1,1,32,1024,32768]
sizeprod2 = [1,1,34,1156,39304]

# for p in nodes:
#     k = p/sizeprod1[3]
#     j = p/sizeprod1[2] - k*sizeprod1[3]/sizeprod1[2]
#     i = p/sizeprod1[1] - j*sizeprod1[2]/sizeprod1[1] - k*sizeprod1[3]/sizeprod1[1]
#     pos.append([i,j,k])
#     print(i,j,k)


# print('-------------')


for p in nodes:
    k = p/sizeprod2[3]
    j = p/sizeprod2[2] - k*sizeprod2[3]/sizeprod2[2]
    i = p/sizeprod2[1] - j*sizeprod2[2]/sizeprod2[1] - k*sizeprod2[3]/sizeprod2[1]
    pos2.append([i,j,k])
    print str((i,j,k))[1:-1]



# with open('new_nodes.csv', 'w') as csvfile:
#     csvwriter = csv.writer(csvfile)

#     for row in pos:
#         csvwriter.writerow(row)