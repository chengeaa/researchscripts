#!/usr/bin/env python
import re
import os
import sys
import numpy as np


def main(filename, reference, metric):
    with open(filename, 'r') as infile:  # Use file to refer to the file object
        with open(filename.replace(".xyz", "2.xyz"), 'w') as outfile:
            xmin, xmax, ymin, ymax = float(xmin), float(xmax), float(ymin), float(ymax) 
            #utility functions
            def tab(n = 4):
                outfile.write(" "*n)
            def newline():
                outfile.write("\n")
            def prettify(string):
                return "{: .16g}".format(float(string)).ljust(18, '0')
            
            atoms = [] #god this is space inefficient but oh well
            count = 0
            max_z_seen = -np.infty

            for i, line in enumerate(infile):
                tokens = line.split()
                if i == 0:
                    continue #do nothing with the first line, 
                    #count needs to be updated anyway
                if i == 1:
                    temp += ["truncated" + line] #comment line stays
                if i > 1: 
                    elem, coords = tokens[0], tokens[1:]
                    x, y, z = [float(i) for i in coords]
                    atoms += [line]
                    count += 1
                    if z > max_z_seen:
                        max_z_seen = z

            for z_curr in np.arange(max_z_seen, 0, stepsize):
                temp = []
                for line in atoms:
                    elem, coords = line[0], tokens[1:]
                    x, y, z = [float(i) for i in coords] 
                    if z < z_curr:
                        temp += [line]

                #TODO: finish this
                metric(temp)

if __name__ == "__main__":
    """
    arg 1: ion damaged file (assume: .xyz)
    arg 2: crystalline reference (we'll see if I use this)
    arg 3: distance metric program (i suppose)?
    """
    args = sys.argv[1:]
    if len(args) > 3:
        print(args)
        raise Exception("Too many args")
    main(*args)

