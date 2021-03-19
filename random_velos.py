import sys
import os
import numpy as np
from pathlib import Path
from ase.io import vasp

np.random.seed(429)

def main(dirname, v, nOut = 1):
    # specify directory

    dirPath = Path(dirname)

    #specify velocity 

    v = float(v) #in vasp units

    #number of outputs per CONTCAR
    nOut = int(nOut)



    for fileName in os.listdir(dirPath):
        if ".out" in fileName:
            continue
        print(dirPath, fileName, Path(fileName))
            
        filePath = dirPath / Path(fileName)
        temp = vasp.read_vasp(filePath)
        headerLen = 9 #nlines of header block; index of first coordinate line

        Arindex = [atom.index for atom in temp if atom.symbol == 'Ar'][0]

        # accounts for header block, all positions + all velos + gap between, goes to last line of velo block
        veloLine = headerLen + len(temp) + Arindex + 1
        arLine = headerLen + Arindex

        for i in range(nOut):
            #phi is angle relative to z axis; theta is angle in xy plane
            #phi between 180 and 90 degrees, theta between 0 and 360
            phi, theta = np.random.uniform(-np.pi, -np.pi/2), np.random.uniform(0, 2*np.pi)

            vx = v * np.sin(phi) * np.cos(theta)
            vy = v * np.sin(phi) * np.sin(theta)
            vz = v * np.cos(phi)

            # random fractional coordinates for x and y positions of Ar
            x, y = np.random.random(size = 2)

            # write output files
            with open(filePath, 'r') as f:
                outputName = fileName + "-" + str(i) + ".out"
                with open(dirPath / Path(outputName), 'w') as fout:
                    for i, line in enumerate(f):
                        if i == veloLine:
                            fout.write(" ".join(np.array([vx, vy, vz], dtype = str)))
                        elif i == arLine:
                            z = line.split()[2] # same z position as input Ar
                            fout.write(" ".join(np.array([x, y, z, "T T T \n"], dtype = str)))
                        else:
                            fout.write(line)

if __name__ == "__main__":
    """
    Given a directory with some number of CONTCARs (after heating), set a random angle for the Ar ion 
    Assumes only one Ar ion
    Arg 1: directory with CONTCARs
    Arg 2: velocity in vasp units
    Arg 3: number of outputs per CONTCAR (optional, default 1)
    """
    args = sys.argv[1:]
    if len(args) > 3:
        print(args)
        raise Exception("Too many args")
    
    if len(args) < 2:
        print(args)
        raise Exception("Too few args, need dirname and v")
    main(*args)
