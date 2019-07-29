#!/usr/bin/env python
import re
import os
import sys
import numpy as np


def main(filename, source, x, y, z):
    with open(filename, 'r') as infile:  # Use file to refer to the file object
        with open(filename.replace(".xyz", ".out"), 'w') as outfile:
            with open(source, 'r') as origin: 
                
                def tab(n = 4):
                    outfile.write(" "*n)
                def newline():
                    outfile.write("\n")
                # def prettify(string):
                #     return ("%7f" % float(string)).ljust(18, '0')
                def prettify(string):
                    return "{: .16g}".format(float(string)).ljust(18, '0')
                
                mults = np.array([x,y,z], dtype = float)

                for i, line in enumerate(origin):
                    tokens = line.split()
                    try:
                        tokens = np.array(tokens, dtype = float)
                    except:
                        pass

                    if i in range(0, 2): #don't touch first two lines
                        outfile.write(line)
                    if i in range(2, 5): #modify lattice params
                        tab(6)
                        outfile.write(" ".join([prettify(i) for i in tokens * mults]))
                        newline()
                    if i == 5: #don't touch atom labels
                        outfile.write(line)
                    if i == 6: #multiply atom numbers
                        tab()                   
                        outfile.write(" ".join([str(int(i)) for i in tokens * np.prod(mults)]))
                        newline()
                    if i == 7: #don't touch dynamics type
                        outfile.write(line)
                outfile.write("Cartesian") #xyz file doesn't do fractional
                newline()
                count = 0
                for i, line in enumerate(infile):
                    tokens = line.split()
                    if i > 1: #don't need to care about first two lines now
                        new = "  " + " ".join([prettify(i) for i in tokens[1:]])
                        if count == 0:
                            new += "   F" * 3
                        else:
                            new += "   T" * 3
                        outfile.write(new)
                        newline()
                        count += 1
                    else:
                        continue
                newline()
                for i in range(count):
                    outfile.write("  " + (("%4e" % float("0")) + " ") * 3)
                    newline()




if __name__ == "__main__":
    """
    Arg 1: .xyz from VESTA and turns into POSCAR.
    Arg 2: Unit cell POSCAR
    Arg 3, 4, 5: xmult, ymult, zmult
    Specifically for unit cell->supercell conversions
    """
    args = sys.argv[1:]
    if len(args) > 5:
        print(args)
        raise Exception("Too many args")
    main(*args)
