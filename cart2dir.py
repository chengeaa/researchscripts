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
                tab()
                outfile.write("\n")
                
            latticeparams = np.array([])
            for i, line in enumerate(infile):
                tokens = line.split()
                try:
                    tokens = np.array(tokens, dtype = float)
                except:
                    pass

                if i in range(0, 2): #don't touch first two lines
                    outfile.write(line)
                if i in range(2, 5): #modify lattice params
                    tab()
                    params = np.array([float(i) for i in tokens])
                    outfile.write(" ".join([str(i) for i in tokens]))
                    latticeparams = np.append(latticeparams, params)
                    newline()
                if i == 5: #don't touch atom labels
                    latticeparams = np.reshape(latticeparams, (3,3))
                    mult = np.linalg.inv(latticeparams)
                    print(mult)
                    outfile.write(line)
                if i == 6: #multiply atom numbers
                    tab()                   
                    outfile.write(" ".join([str(int(i)) for i in tokens]))
                    newline()
                if i == 7: #set as direct 
                    outfile.write("Direct \n")
                if i == 8:
                    outfile.write(line)
                if i > 8:
                    inlen = len(tokens)
                    tokens = tokens[:3]
                    try:
                        tokens = np.array(tokens, dtype = float)
                        result = np.matmul(tokens, mult)
                        outfile.write(" ".join([str(i) for i in result]))
                        if inlen > 3:
                            outfile.write("   T"*3)
                        newline()
                    except:
                        outfile.write(line)




if __name__ == "__main__":
    """
    Arg 1: POSCAR
    """
    args = sys.argv[1:]
    if len(args) > 1:
        print(args)
        raise Exception("Too many args")
    main(*args)