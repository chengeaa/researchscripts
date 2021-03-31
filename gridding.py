# imports 
# base python
import os
import copy
from sys import getsizeof

# scipy
import numpy as np

#ase
from ase.io import gen, vasp, xyz, extxyz
from ase.visualize import view
from ase.build import make_supercell
from ase.visualize.plot import plot_atoms
from ase.build import add_adsorbate
from ase.geometry.analysis import Analysis


adsorbate = vasp.read_vasp("adsorbate here ")
slab = vasp.read_vasp("slab here ")


i = 1
for _x in np.arange(0, 1, .2):
    for _y in np.arange(0, 1, .2):
        s = slab.copy()
        x, y, z = slab.cell.cartesian_positions([_x, _y, 0])
        add_adsorbate(s, adsorbate, height = 2, position = (x, y))
        s = Atoms(sorted(s, key = lambda x: x.symbol))
        s.cell = slab.cell
        vasp.write_vasp("POSCAR" + str(i), s, 
                        label = 'label goes here', 
                        vasp5=True)
        
        i += 1
