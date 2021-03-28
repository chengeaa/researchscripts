import numpy as np
import pandas as pd
import re
from ase.io import gen


def normalize(y,x):
    """
    Takes y, x of data and returns normalized y
    """
    return y/np.trapz(y,x)

def KE(v_tot):
    "Returns KE of Ar+ in eV given total velocity"
    return 6.24E18 * 0.5 * 1.66E-27*39.95*(v_tot*1E5)**2

def v_from_KE(E):
    "Returns v(z) of Ar+ in eV given KE"
    return np.sqrt(E/(6.24E18 * 0.5 * 1.66E-27*39.95))/1E5

def readStructs(datadir, shallow = True, name = "output"):
    """
        Currently designed for output from single layer directory trees
        Reads in final adsorption geometries and energy data,
            returns dataframe with geometry and energy data

        Input:
            datadir: string that points to directory containing the following:
                - convergence: each line i has convergence status of run i
                - energies: each line i has total energy and ads energy from run i
                - output{indices}.gen: final geometries for each index
            slabEnergy: energy of slab
            adsorbateEnergy: energy of the adsorbate in the system
        Output:
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
 
#############
# constants #
#############
amu2g = 1.66054e-24 #multiply by this
a2cm = 1e-8 # multiply by this
