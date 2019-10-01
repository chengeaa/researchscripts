#!/usr/bin/env python
import re
import os
import sys
import numpy as np


def main(filename, xmin, xmax, ymin, ymax):
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
            
            temp = []#god this is space inefficient but oh well
            count = 0
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
                    if xmin <= x <= xmax and ymin <= y <= ymax:
                        temp += [line]
                        count += 1

    with open(filename.replace(".xyz", "2.xyz"), 'w') as outfile:
        outfile.write(str(count))
        newline()
        for line in temp:
            outfile.write(line)

if __name__ == "__main__":
    """
    Arg 1: .xyz 
    Arg 2, 3, 4, 5: xmin, xmax, ymin, ymax
    For trimming xyz files within x y bounds.
    """
    args = sys.argv[1:]
    if len(args) > 5:
        print(args)
        raise Exception("Too many args")
    main(*args)

