"""
Translate a POSCAR some amount in x,y,z directions for easier visualization
"""

#!/usr/bin/env python
import re
import os
import sys
import numpy as np


def main(filename, x, y, z):
    with open(filename, 'r') as infile:  # Use file to refer to the file object
        with open(filename +".out", 'w') as outfile:

            def tab():
                outfile.write(" "*4)
            def newline():
                tab()
                outfile.write("\n")
                
            def prettify(string):
                return "{: .16g}".format(float(string)).ljust(18, '0')
            translations = np.array([x,y,z], dtype = float)

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
                    outfile.write(" ".join([str(i) for i in tokens]))
                    newline()
                if i == 5: #don't touch atom labels
                    outfile.write(line)
                if i == 6: #multiply atom numbers
                    tab()                   
                    outfile.write(" ".join([str(int(i)) for i in tokens]))
                    newline()
                if i == 7: #don't touch dynamics type
                    outfile.write(line)
                if i == 8:
                    outfile.write(line)
                if i > 8:
                    inlen = len(tokens)
                    tokens = tokens[:3]
                    try:
                        tokens = np.array(tokens, dtype = float)
                        
                        if inlen > 3:
                            outfile.write("  " + " ".join([str(i) for i in tokens + translations]))
                            outfile.write("   T"*3)
                        else:
                            outfile.write(" ".join([str(i) for i in tokens]))
                        newline()
                    except:
                        outfile.write(line)




if __name__ == "__main__":
    """
    Arg 1: POSCAR
    Arg 2, 3, 4: xtrans, ytrans, ztrans 
    """
    args = sys.argv[1:]
    if len(args) > 4:
        print(args)
        raise Exception("Too many args")
    main(*args)
