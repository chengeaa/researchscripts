#!/usr/bin/env python
import re
import os
import sys
import numpy as np
from ase.io import xyz, vasp
from ase.cell import Cell



def main(filename, source, x=1, y=1, z=1):
    temp = list(xyz.read_xyz(filename))[0]
    
    if source == "NA":
        temp.cell = Cell([[x, 0, 0], [0, y, 0], [0, 0, z]])

    else:
        try:
            temp.cell = vasp.read_vasp(source).cell
        except:
            raise Exception("source is neither valid POSCAR nor NA")

    vasp.write_vasp("../POSCAR_ec", temp, sort = True, vasp5=True)

if __name__ == "__main__":
    """
    Arg 1: .xyz from VESTA and turns into POSCAR.
    Arg 2: Path for POSCAR source OR NA, let args 3-5 be cubic cell inputs
    Arg 3, 4, 5: xmult, ymult, zmult
    Specifically for unit cell->supercell conversions
    """
    args = sys.argv[1:]
    if len(args) > 5:
        print(args)
        raise Exception("Too many args")
    main(*args)
