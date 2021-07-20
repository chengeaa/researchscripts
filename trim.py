#!/usr/bin/env python3
'''
Designed to work on the cluster, removing 'ejected species' after each bomb or quench step
'''
# imports 

import os
import re
import sys

# scipy
import numpy as np
import pandas as pd

#ase
from ase.io import gen, vasp, xyz, extxyz, dftb



def main(
    surftype, #type of surface as corresponding to keys in surfzmaxes
    datadir = "temp/", #data files, structured as datadir/output$i-$j.gen and datadir/velos$i-$j
    outputdir = "temp.new/",  #files for output
    hbondrange = 4, #offset from surface corresponding to Hbond range
    zmincutoff = 0.1, #somewhat arbitrary value to get rid of atoms that have gone into bulk
    output_geom_name = "output",  #prefix for output geometry files
    output_velos_name = "velos" #prefix for output velocity files
    ):

    surfzmaxes = {'nrichhterm': 14.405, 'sirichfterm': 13.229, 'modded':15.00723826, 'amorphous':18.622, 'tol_amorphous':32.979, 'pure_si': 22.96, 'sirichfterm_amorphous': 18.3} #xtl reference slab zmax
    surfabbr = {'nh':'nrichhterm', 'sf': 'sirichfterm', 'a': 'amorphous', 'tol': 'tol_amorphous', 'pure_si':'pure_si', 'sif_a':'sirichfterm_amorphous'}
    if surftype in surfabbr.keys():
        surftype = surfabbr[surftype]
        xtlzmax = surfzmaxes[surftype] 
    else:
        try:
            xtlzmax = float(surftype)
        except TypeError:
            print("surftype arg must either be a key within the known surface types, or a float representing the max surface height")
            raise
    zcutoff = xtlzmax + hbondrange # approximating H bonding range as 4 Å

    ##############################
    ### Read in geometry files ###
    ##############################


    geometries = {}
    for i in os.listdir(datadir):
        if output_geom_name in i:
            key = re.search(r"\d+", i)
            if key:
                key = key.group(0)
                geometries[key] =  gen.read_gen(datadir + i)

    ##########################
    ### Read in velocities ###
    ##########################
    velos = dict()
    for i in os.listdir(datadir):
        if output_velos_name in i:
            key = re.search(r"\d+", i)
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

    #write 
    for key, geom in trimmedgeoms.items():
        gen.write_gen("%sinput%s.gen" % (outputdir, key), 
            geom)
    for key, v in trimmedvelos.items():
        v.to_csv("%s%s%s.in" % (outputdir, output_velos_name, key), 
            sep = " ", index = False, header = False)

if __name__ == "__main__":
    """
    Takes in three arguments:
    surftype: original slab type (used for zmax ref, some slabs are taller than others)
        if this is a number, use this as the zmax ref
    datadir: where the input data is
    outputdir: where the output goes
    """

    args = sys.argv[1:]
    if len(args) > 3:
        print(args)
        raise Exception("No more than 3 arguments allowed")
    main(*args)
