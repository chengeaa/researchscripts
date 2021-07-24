#!/usr/bin/env python3
'''
Designed to work on the cluster, removing 'ejected species' after each bomb or quench step
Make sure the source files are named in the form `output%d.gen`
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
from ase.geometry.analysis import Analysis

#researchscripts
from researchscripts.structure import Graph


def main(
    surftype, #type of surface as corresponding to keys in surfzmaxes
    datadir = "temp/", #data files, structured as datadir/output$i-$j.gen and datadir/velos$i-$j
    outputdir = "temp.new/",  #files for output
    hbondrange = 3, #offset from surface corresponding to Hbond range
    zmincutoff = 0.1, #somewhat arbitrary value to get rid of atoms that have gone into bulk
    output_geom_name = "output",  #prefix for output geometry files
    output_velos_name = "velos" #prefix for output velocity files
    ):

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
        removedatoms = {'Si': 0, 'N': 0, 'H': 0, 'Ar': 0, 'F':0, 'C':0}

        # construct graph 
        adjmat = Analysis(geom).adjacency_matrix[0]
        numnodes = adjmat.shape[0]
        g = Graph(numnodes)
        for i in range(numnodes):
            for j in range(numnodes):
                if adjmat[i,j]:
                    g.addEdge(i,j)
        cc = g.connectedComponents()

        #identify slab, and max height of slab
        maingraph = np.array([i for i in cc if 0 in i][0])
        slab = geom[[atom.index for atom in geom if atom.index in maingraph]]
        gen.write_gen(outputdir + "slab{}.gen".format(key), slab)
        zcutoff = np.max([atom.position[2] for atom in slab]) + hbondrange
        
        # isolate fragments and identify which to remove
        fragGraphs = [i for i in cc if 0 not in i]
        fragZs = [[geom[i].position[2] for i in frag] for frag in fragGraphs]
        removeFrag = [np.all(np.array(i) > zcutoff) or np.all(np.array(i) < zmincutoff) 
                for i in fragZs]
        atomsToRemove = [i for g,r in zip(fragGraphs, removeFrag) if r for i in g]

        for idx in atomsToRemove:
            removedatoms[geom[idx].symbol] += 1 #tally removed atoms by species
        
        geomcopy = geom.copy()
        del geomcopy[[atom.index for atom in geomcopy if atom.index in atomsToRemove]]
        
        removedspecies[key] = pd.Series(removedatoms)
        trimmedgeoms[key] = geomcopy
        trimmedvelos[key] = velos[key][[i not in atomsToRemove for i in np.arange(len(velos[key]))]]

        
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
