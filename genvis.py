#!/usr/bin/env python3
'''
View .gen file
'''

#imports
import sys
from ase.io import gen, vasp, xyz, extxyz
from ase.visualize import view



def main(file):
    view(gen.read_gen(file))

if __name__ == "__main__":
    """
    Takes in one argument:
    file: name of file (or path to file)
    """

    args = sys.argv[1:]
    if len(args) > 1:
        print(args)
        raise Exception("No more than 1 arguments allowed")
    main(*args)
