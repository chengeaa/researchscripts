#!/usr/bin/env python
'''Performs some kind of parsing on POSCAR file!...?

Erik won't tell me what this does, so all I can say is that it parses the file passed in
and writes some modified version to an outfile.

Args:
    POSCAR - filepath for POSCAR file

Returns:
    None; writes results to "{POSCAR}.out" based on the POSCAR filepath provided as an argument.

'''

import re
import os
import sys
import numpy as np
import argparse


def main(filename):
    with open(filename, 'r') as infile:  # Use file to refer to the file object
        with open(filename + ".out", 'w') as outfile:

            def tab():
                outfile.write(" "*4)
            def newline():
                tab()
                outfile.write("\n")
                
            latticeparams = np.array([])
            for i, line in enumerate(infile):
                tokens = line.split()
                try:
                    tokens = np.array(tokens, dtype=float)
                except:
                    pass

                if i in range(0, 2):  # Don't touch first two lines
                    outfile.write(line)
                if i in range(2, 5):  # Modify lattice params
                    tab()
                    params = np.array([float(i) for i in tokens])
                    outfile.write(" ".join([str(i) for i in tokens]))
                    latticeparams = np.append(latticeparams, params)
                    newline()
                if i == 5:  # Don't touch atom labels
                    latticeparams = np.reshape(latticeparams, (3, 3))
                    mult = np.linalg.inv(latticeparams)
                    print(mult)
                    outfile.write(line)
                if i == 6:  # Multiply atom numbers
                    tab()                   
                    outfile.write(" ".join([str(int(i)) for i in tokens]))
                    newline()
                if i == 7:  # Set as direct 
                    outfile.write("Direct \n")
                if i == 8:
                    outfile.write(line)
                if i > 8:
                    inlen = len(tokens)
                    tokens = tokens[:3]
                    try:
                        tokens = np.array(tokens, dtype=float)
                        result = np.matmul(tokens, mult)
                        outfile.write(" ".join([str(i) for i in result]))
                        if inlen > 3:
                            outfile.write("   T"*3)
                        newline()
                    except:
                        outfile.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Performs some kind of parsing on POSCAR file(?)")
    parser.add_argument('POSCAR', help='filepath of POSCAR input file', type=str)
    args = parser.parse_args()
    main(*args)
