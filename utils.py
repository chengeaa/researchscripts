import os
import numpy as np
import pandas as pd
import re
from ase.io import gen


def normalize(y,x):
    """
    Takes y, x of data and returns normalized y
    """
    return y/np.trapz(y,x)

def flatten(l):
    """
    flattens a 2d list
    """
    return [item for sublist in l for item in sublist] 

def KE(v_tot):
    """
    Returns KE of Ar+ in eV given total velocity
    """
    return 6.24E18 * 0.5 * 1.66E-27*39.95*(v_tot*1E5)**2

def v_from_KE(E):
    """
    Returns v(z) of Ar+ in eV given KE
    """
    return np.sqrt(E/(6.24E18 * 0.5 * 1.66E-27*39.95))/1E5

def readStructs(datadir, shallow = True, name = "output"):
    """
        Currently designed for output from single layer directory trees.
        Reads in final adsorption geometries and energy data, returns dataframe with geometry and energy data

        Input:
            datadir: string that points to directory containing the following:
                - convergence: each line i has convergence status of run i
                - energies: each line i has total energy and ads energy from run i
                - output{indices}.gen: final geometries for each index
                
            slabEnergy: energy of slab
            adsorbateEnergy: energy of the adsorbate in the system

        Returns:
            output: pd Dataframe with:
                - index: indices for runs that worked
                - geometry: final geometry of run
                - total energy: raw energy from file
                - adsorption energy: energy as adjusted by adsorbate_energy
    """
    geometries = {}
    if shallow:
        pattern = r"{}(\d+).gen".format(name)
    else:
        pattern = r"{}(\d+-\d+).gen".format(name)
    files = os.listdir(datadir) 

    if "energies" in files and "convergence" in files:
        convergence = pd.read_csv(datadir + "convergence", header = None)
        energies = pd.read_csv(datadir + "energies", header = None)
        output =  pd.concat([energies, convergence], axis = 1)
        output.columns = ["E", "E_ads", "conv"]

        for i in files:
            key = re.search(pattern, i)
            if key:
                if shallow:
                    key = int(key.group(1))
                else:
                    key = key.group(1)
                geometries[key] =  gen.read_gen(datadir + i)
        output['geom'] = pd.Series(geometries)

        output = output[output['conv'] == "Geometry converged"]
        output = output.drop("conv", axis = 1)

    else:
        for i in files:
            key = re.search(pattern, i)
            if key:
                if shallow:
                    key = int(key.group(1))
                else:
                    key = key.group(1)
                geometries[key] =  gen.read_gen(datadir + i)
        output = pd.DataFrame(pd.Series(geometries))
        output.columns = ['geom']
    return output

def readData(outdir, indir, useRemoval = True, useFrags = True, useBonds = True, wrapStructs = True):
    """
    Specifically made for data with BOTH input and output structures
    Read structures; append removal data, fragment data, and bond data optionally
    Wrap structures optionally

    Returns:
        df of the above with columns (struct, in), (struct, out), 
        optionally, (frags, <fragnames>), (bonds, <bondnames>), (removal, <removedelems>)
        where each <x> corresponds to a column label for a specific subgroup (ie, NH3, N-H, or NH3 resp.)
    """

    # read structures, initialize df

    data = pd.concat([readStructs(indir, shallow = False, name = 'input'),
                         readStructs(outdir, shallow = False)
                        ], axis = 1)
    data.columns = ['in', 'out']
    for i in data['in']:
        i.wrap()

    if useRemoval:
        bombdata, quenchdata, eqdata = postprocessResults(outdir)

        newdata = {}
        for idx in data.index:
            i, j = idx.split("-") # for i-j type output
            i = int(i)
            newdata[idx] = data.loc[idx].append(bombdata[i][j])
        data = pd.DataFrame(newdata).T
        ncols = data.shape[1]
        data.columns = pd.MultiIndex.from_arrays(
            [["struct"] * 2 + ["removal"] * (ncols - 2), data.columns],
            names = ['source', 'data']
        )
    else:
        data.columns = pd.MultiIndex.from_arrays(
            [["struct"] * 2, data.columns],
            names = ['source', 'data']
        )


    if useFrags:
        infrags = pd.read_csv(indir+"fragdata.csv", index_col=0)
        outfrags = pd.read_csv(outdir+"fragdata.csv", index_col=0)

        fragdiffs = outfrags.subtract(infrags, fill_value = 0)
        ncols = fragdiffs.shape[1]
        fragdiffs.columns = pd.MultiIndex.from_arrays(
            [['frags'] * ncols, fragdiffs.columns], names = ['source', 'data']
        )

        data = pd.concat([data, fragdiffs], axis = 1)


    if useBonds:
        outbonds = pd.read_csv(outdir + "bondcounts.csv", index_col=0)

        inbonds = pd.read_csv(indir + "bondcounts.csv", index_col=0)

        bondDiffs = outbonds.subtract(inbonds)

        ncols = bondDiffs.shape[1]
        bondDiffs.columns = pd.MultiIndex.from_arrays(
            [['bonds'] * ncols, bondDiffs.columns], names = ['source', 'data']
        )

        data = pd.concat([data, bondDiffs], axis = 1)
    return data

def calculateVASPtemp(velos, nfixed, POSCAR):
    """
    Take in a 3N matrix of velocities, subtract of nfixed from DOF
    POSCAR (Atoms object) used for getting masses from atoms
    Return T as calculated by VASP
    """
    masses = [atom.mass for atom in POSCAR]
    KE = 0
    dof = len(velos) * 3 - nfixed   #natoms * 3 - 6 (fixed atoms) * 3
    for i, v in enumerate(velos):
        vx, vy, vz = v
        m = masses[i]
        v = np.sqrt(vx**2 + vy**2  + vz**2)
        KE += 1/2 * m * v**2
    KE *= 103.642697 # amu * Å^2/fs^2 to eV conversion
    return 2 * KE/kb_evK/dof

#############
# constants #
#############
amu2g = 1.66054e-24 #multiply by this
a2cm = 1e-8 # multiply by this
kb_evK = 8.617333262145E-5 #ev/K
