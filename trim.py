#!/usr/bin/env python3
'''
Designed to work on the cluster, removing 'ejected species' after each bomb or quench step
'''
# imports 
# base python
import os
import copy
import re
import sys
import random

# scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt,mpld3
import seaborn as sns
from pathlib import Path
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d, interp2d
from sklearn.metrics.pairwise import rbf_kernel
from sklearn import preprocessing
import matplotlib.tri as tri

#ase
from ase.io import gen, vasp, xyz, extxyz, dftb
from ase.io.dftb import read_dftb_velocities, write_dftb_velocities
from ase.calculators.dftb import Dftb
from ase import Atoms, Atom
from ase.constraints import FixAtoms
from ase.visualize import view
from ase.build import make_supercell
from ase.visualize.plot import plot_atoms
from ase.build import add_adsorbate
import nglview
from ase.geometry.analysis import Analysis



def main(datadir = "temp/", #data files, structured as datadir/output$i-$j.gen and datadir/velos$i-$j
    outputdir = "../temp.new/", #files for output
    surftype = 'nh'):

    surfzmaxes = {'nrichhterm': 14.405, 'sirichfterm': 13.229} #xtl reference slab zmax
    surfabbr = {'nh':'nrichhterm', 'sf': 'sirichfterm'}
    surftype = surfabbr[surftype]
    xtlzmax = surfzmaxes[surftype] 
    zcutoff = xtlzmax + 4 # approximating H bonding range as 4 AA
    zmincutoff = 0.1 #somewhat arbitrary value to get rid of atoms that have gone into bulk

    numsims = 170 # number of simulations per batch

    ##############################
    ### Read in geometry files ###
    ##############################

    geometries = {}
    for i in os.listdir(datadir):
        if "output" in i:
            key = re.search(r"\d+-\d+", i)
            if key:
                key = key.group(0)
                geometries[key] =  gen.read_gen(datadir + i)

    ##########################
    ### Read in velocities ###
    ##########################
    velos = dict()
    for i in os.listdir(datadir):
        if "velos" in i:
            key = re.search(r"\d+-\d+", i)
            if key:
                key = key.group(0)
                velos[key] = pd.read_csv(datadir + i, header = None, dtype = float, sep = "\s+")


    ################
    ### trimming ###
    ################


    trimmedgeoms = dict()
    trimmedvelos = dict()

    removedspecies = dict()

    for key, geom in geometries.items(): 
        aboveZindices = set()
        belowZindices = set()
        otherindices = set()
        nearbyatoms = set()
        removedatoms = {'Si': 0, 'N': 0, 'H': 0, 'Ar': 0, 'F':0, 'C':0}
        for atom in geom:
            if atom.position[2] > zcutoff:
                aboveZindices.add(atom.index)
            elif atom.position[2] < zmincutoff:
                belowZindices.add(atom.index)
            else:
                otherindices.add(atom.index)
        # iteration through geom once guarantees uniqueness in aboveZindices and otherindices 
        
        for i in aboveZindices:
            _dists = geom.get_distances(i, list(otherindices))
            nearbyatoms.update(np.array(list(otherindices))[_dists < 2]) 
            # add indices where distance to i is less than 2
        atomstoremove = aboveZindices.union(nearbyatoms).union(belowZindices)
        
        for idx in atomstoremove:
            removedatoms[geom[idx].symbol] += 1 #tally removed atoms by species
        
        geomcopy = geom.copy()
        del geomcopy[[atom.index for atom in geomcopy if atom.index in atomstoremove]]
        

        removedspecies[key] = pd.Series(removedatoms)
        trimmedgeoms[key] = geomcopy
        trimmedvelos[key] = velos[key][[i not in atomstoremove for i in np.arange(len(velos[key]))]]

        
    # collect all removed species series into a df and write as csv
    pd.DataFrame(removedspecies).to_csv("removedspecies.csv")

    for key, geom in trimmedgeoms.items():
        gen.write_gen("%sinput%s.gen" % (outputdir, key), 
            geom)
    for key, v in trimmedvelos.items():
        v.to_csv("%svelos%s.in" % (outputdir, key), 
            sep = " ", index = False, header = False)

if __name__ == "__main__":
    """
    Takes in three arguments:
    datadir: where the input data is
    outputdir: where the output goes
    surftype: original slab type (used for zmax ref, some slabs are taller than others)
    """

    args = sys.argv[1:]
    if len(args) > 3:
        print(args)
        raise Exception("No more than 3 arguments allowed")
    main(*args)
