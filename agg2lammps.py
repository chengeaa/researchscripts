"""
Write 'aggcar' (positions +  velocities) to a lammps file, for visualizing heat dissipation
"""

#!/usr/bin/env python
import re
import os
import sys
import numpy as np


def main(filename):
    with open(filename, 'r') as infile:  # Use file to refer to the file object
        with open(filename +".out", 'w') as outfile:

            def tab():
                outfile.write(" "*4)
            def newline():
                outfile.write("\n")
            def write(text):
                outfile.write(text)
                newline()
            natoms = 0
            timestep = 0
            latticeparams = np.array([])
            for i, line in enumerate(infile):
                tokens = line.split()
                if len(tokens) == 1:
                    natoms = tokens[0]
                    write("ITEM: TIMESTEP")
                    write(str(timestep))
                    timestep += 1
                    write("ITEM: NUMBER OF ATOMS")
                    write(natoms)
                    write("ITEM: BOX BOUNDS pp pp pp")
                    for i in range(3):
                        write("-7.5 7.5")
                    write("ITEM: ATOMS element x y z vx vy vz")
                elif len(tokens) == 5 or len(tokens) == 4:
                    #skip lines with E and F data (5)
                    #and lines with no velocities (ie, last frame) (4)
                    continue
                elif len(tokens) == 7:
                    outfile.write(" ".join(tokens))
                    newline()
                else:
                    print("line length = ", len(tokens))
                    raise Exception("Unexpected line length at line " + str(i))




if __name__ == "__main__":
    """
    Arg 1: AGGCAR
    """
    args = sys.argv[1:]
    if len(args) > 1:
        print(args)
        raise Exception("Too many args")
    elif len(args) == 1:
        f = args[0]
    else:
        f = 'AGGCAR'
    main(f)
