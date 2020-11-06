import numpy as np
import matplotlib.pyplot as plt
file1 = open("../../data/PINC.out", "r")
file2 = open("../../data/PINC1.out", "r") # 64ppc


def getPotFromFile(infile):
	ObjPot = []
	for line in infile:
		words = line.split(" ")
		for i in range(2,len(words)):	
			if (words[i-2] == "object" and words[i-1] == "0" and words[i] == ":" ):
				ObjPot.append(words[i+1][:-1])
	return np.asarray(ObjPot, dtype = float)



arr1 = getPotFromFile(file1)
arr2 = getPotFromFile(file2) # 64 ppc


plt.plot(arr1, label="6ppc")
plt.plot(arr2, label="64ppc")

plt.legend()

plt.show()
