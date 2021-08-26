#!/usr/bin/env python3
'''
Designed to work on the cluster, removing 'ejected species' after each equilibration step
Make sure the source files are named in the form `output%d.gen`
'''
# imports 
import os
import re
import sys

import numpy as np
import pandas as pd

#ase
from ase.io import gen, vasp, xyz, extxyz, dftb
from ase.build import add_adsorbate
from ase.geometry.analysis import Analysis

#researchscripts
from researchscripts.structure import Graph


def main(
    numsofar, #use the run number you're seeding for
    batch, #current batch number
    velo, # velocity of incident Ar in Å/ps
    datadir = "temp/", #data files, structured as datadir/output$i-$j.gen and datadir/velos$i-$j
    outputdir = "temp.new/", #files for output
    hbondrange = 3,
    zmincutoff = 0.1, #somewhat arbitrary value to get rid of atoms that have gone into bulk
    numperbatch = 17,
    numbatches = 10 
    ):

    numsofar = int(numsofar)
    batch = int(batch)
    velo = int(velo)

    ##############################
    ### Read in geometry files ###
    ##############################

    geometries = {}
    for i in os.listdir(datadir):
        if "output" in i:
            key = re.search(r"\d+", i)
            if key:
                key = key.group(0)
                geometries[key] =  gen.read_gen(datadir + i)

    ##########################
    ### Read in velocities ###
    ##########################
    velos = dict()
    for i in os.listdir(datadir):
        if "velos" in i:
            key = re.search(r"\d+", i)
            if key:
                key = key.group(0)
                velos[key] = pd.read_csv(datadir + i, header = None, dtype = float, sep = "\s+")

    # to account for seed behavior from first numssofar sets of runs
    # numssofar can also be interpreted as = current run seeding for
    np.random.seed(429)
    for b in range(batch + numsofar * numbatches):
        for i in range(numperbatch):
            x_rand, y_rand, z_rand = np.append(np.random.random(size = 2), 0)

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

        x_rand, y_rand, z_rand = geomcopy.cell.cartesian_positions(np.append(np.random.random(size = 2), 0))

        add_adsorbate(geomcopy, adsorbate = 'Ar', height = 7, position = (x_rand, y_rand))

        removedspecies[key] = pd.Series(removedatoms)
        trimmedgeoms[key] = geomcopy
        trimmedvelos[key] = velos[key][[i not in atomsToRemove for i in np.arange(len(velos[key]))]]
        trimmedvelos[key] = trimmedvelos[key].append(pd.Series([0, 0, -velo]), ignore_index=True)


        
    # collect all removed species series into a df and write as csv
    pd.DataFrame(removedspecies).to_csv("removedspecies.csv")

    #write 
    for key, geom in trimmedgeoms.items():
        gen.write_gen("%sinput%s.gen" % (outputdir, key), 
            geom)
    for key, v in trimmedvelos.items():
        v.to_csv("%svelos%s.in" % (outputdir, key), 
            sep = " ", index = False, header = False)

if __name__ == "__main__":
    """
    Takes in four arguments:
    numsofar: number of runs so far; alternatively, run number seeding for
    batch: batch number (this script is designed for 2-level batch/run structure)
    velo: velocity 
    datadir: where the input data is
    outputdir: where the output goes
    """

    args = sys.argv[1:]
    if len(args) > 5:
        print(args)
        raise Exception("No more than 5 arguments allowed")
    main(*args)
