import numpy as np
import os

def main():
    files = []

    for fname in os.listdir("./potentials"):
        if fname.endswith("pot.txt"):
            files.append(fname)
    
    print(files)
    for fname in files:
        potential = np.loadtxt(f'./potentials/{fname}', delimiter='\n')
        potential = potential[999:]
        #print(len(potential))
        avgPot = np.mean(potential)
        print(f"Floating potential of {fname} is: {avgPot}")

if __name__ == "__main__":
    main()