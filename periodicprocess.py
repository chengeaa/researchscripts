"""
Not sure what this does honestly
"""

#!/usr/bin/env python
import re
import os
import sys
import numpy as np

import time

def main(filename, x, y):
    x, y = int(x), int(y)
    mults = np.array([x, y, 1])
    with open(filename, 'r') as infile:  # Use file to refer to the file object
        with open(filename +".xyz", 'w') as outfile:

            def tab():
                outfile.write(" "*4)
            def newline():
                tab()
                outfile.write("\n")
            def write(line):
                outfile.write(line)
                newline()
            def prettify(string):
                return "{: .16g}".format(float(string)).ljust(19, '0')

            
            # results = []
            donttouchlines = np.array([1, 2, 6]) #lines with like text that I don't wanna touch
            lveclines = np.array([3, 4, 5]) #lines with lattice vectors
            for i, line in enumerate(infile):
                # if i == 9:
                #     return
                tokens = line.split()
                i += 1 #vi line number <- 0 index conversion
                if i <= 2:
                    continue #ignore the first two comment lines
                if i in lveclines: #lattice vectors
                    result = [float(token) for token in tokens] * np.array([x, y, 1])
                    if i == 3:
                        xvec = np.array(result)
                    if i ==4:
                        yvec = np.array(result)
                    if i ==5:
                        zvec = np.array(result)

                elif i == 6:
                    atomlabels = np.array(tokens)

                elif i == 7: #atom counts
                    result = [int(token) * x * y for token in tokens] 

                    atomlabels = np.array([((label+" ") * int(count)) for label, count in zip(atomlabels, result)])
                    atomlabels = np.concatenate([string.split() for string in atomlabels])

                elif "Direct" in tokens:
                    write(str(len(atomlabels)))
                    write("Frame:" + tokens[-1])
                    atomcount = 0

                elif len(tokens) == 3 and i > 8 and "Direct" not in tokens:
                    def dir2cart(atom):
                        """
                        Atom should be a 3 element numpy array 
                        """
                        basis = np.array([xvec, yvec, zvec])
                        return np.matmul(basis.T, atom.T)
                    base = np.array([float(token) for token in tokens]) * 1/mults
                    
                    newatoms = [dir2cart(base)]
                    xshifted = []
                    for j in range(1, x):
                        for atom in newatoms[:]:
                            xshifted +=  [(np.array(atom) + np.array([j, 1, 1]) * xvec/x)]
                    newatoms += xshifted
                    yshifted = []
                    for j in range(1, y):
                        for atom in newatoms[:]:
                            yshifted +=  [(np.array(atom) + np.array([1, j, 1]) * yvec/y)]
                    newatoms += yshifted
                    newatoms = np.reshape(newatoms, (-1, 3))

                    for atom in newatoms:
                        write(atomlabels[atomcount] + " " * 4 + " ".join([str(coord) for coord in atom]))
                        atomcount += 1
                else:
                    raise Exception("unrecognized line format at line", i)
                
                # print("Line %d: %s" % (i, "hi"))

            # for i in results:
            #     outfile.write(i)
            #     newline()



if __name__ == "__main__":
    """
    Arg 1: XDATCAR
    Arg 2, 3: xmult, ymult
    No zmult allowed for now
    """
    args = sys.argv[1:]
    if len(args) > 3:
        print(args)
        raise Exception("Too many args")
    main(*args)
