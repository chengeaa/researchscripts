#!/usr/bin/env python3
'''
Convert from .gen to POSCAR
'''

#imports
import sys
from ase.io import gen, vasp, xyz, extxyz
from ase.visualize import view



def main(file, output):
    vasp.write_vasp(output, gen.read_gen(file), sort = True, vasp5 = True)

if __name__ == "__main__":
    """
    Takes in one argument:
    file: name of file (or path to file)
    output: name of output file (or path to output file)
    """

    args = sys.argv[1:]
    if len(args) > 2:
        print(args)
        raise Exception("No more than 1 arguments allowed")
    main(*args)
