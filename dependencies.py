###########
# imports # 
###########
# base python
import os
import copy
import re
from sys import getsizeof
import random

# scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt,mpld3
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import seaborn as sns
from pathlib import Path
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d, interp2d
from sklearn.metrics.pairwise import rbf_kernel
from sklearn import preprocessing
import matplotlib.tri as tri

plt.rcParams["figure.figsize"] = (20,7)

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

#dscribe
from dscribe.descriptors import SOAP
from dscribe.descriptors import MBTR
from dscribe.kernels import REMatchKernel
from dscribe.kernels import AverageKernel

from sklearn import preprocessing


#quippy 
from ase.build import bulk
from ase.optimize import LBFGS
from ase.visualize import view
from quippy.potential import Potential


#misc
import similaritymeasures

#############
# functions #
#############

def show_atoms_grid(data, rotation = '-0x,0y,0z', save= False, filename = 'grid_configs'):
    '''
    Where data is list of Atoms objects
    '''
    dim = int(np.ceil(np.sqrt(len(data))))
    fig, axarr = plt.subplots(dim, dim, figsize=(25, 25))
    for i, config in enumerate(data):
        plot_atoms(config, axarr[i%dim,i//dim], rotation = rotation)
    if save:
        fig.savefig(filename + ".png")

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

def coordLabeller(atoms, image, fullCoordinations = {"Si": 4, "N":3, "H":1, "F": 1},
                  minAngles = {"Si": 90, "N": 109.5, "H": 360, "F": 360}, # 360 for atoms with max 1 bond
                  maxBonds_per_element = {"Si": 6, "N":4, "H":1, "F":1},
                  angle_tolerance = 0, #tolerance for valid bond angles
                  bond_tolerance = 0.5, #tolerance for valid bond angles
                  minz = 0 #minimum height above which to compute
                 ):
    """
    Takes a structure, returns two dictionaries, the keys of which are identical:
        the index of the atom for which the statistic is calculated
    relativeCoordinations: -1 if atom i is undercoordinated, 1 if overcoordinated, 0 else
    bonds: list of bonds for atom i
    """
#     atomAnalysis = Analysis(atoms)
#     atomBonds = atomAnalysis.all_bonds[0]
#     neighbors = build_neighbor_list(atoms, bothways = True, self_interaction = False)
    nl = NewPrimitiveNeighborList(
        cutoffs = np.array(natural_cutoffs(atoms)) * (bond_tolerance + 1),
        bothways = True,
        self_interaction = False)
    nl.build(pbc = [True, True, False] , cell = atoms.cell, positions = atoms.positions)
    coordinations = {}
    relativeCoordinations = {}
    newBonds = {}
    for atom in atoms:
        idx = atom.index
        if atom.symbol == 'Ar' or atom.position[2] < minz:
            coordinations[idx] = 0
            relativeCoordinations[idx] = 0
            newBonds[idx] = []
            continue
        minAngle = minAngles[atom.symbol] * (1 - angle_tolerance) # set minimum required angle to keep
        bonds = nl.get_neighbors(idx)[0]
        maxBonds = maxBonds_per_element[atom.symbol]

        bonds = sorted(bonds,
#                        key = lambda b: atomAnalysis.get_bond_value(image, [idx, b]),
                       key = lambda b: atoms.get_distance(idx, b, mic = True),
                       reverse = False) #sort bonds list in ascending order of bond length
        keptBonds = set()
        skipBonds = set()
        if len(bonds) == 0:
            print("No bonds detected for atom %d" % idx)
        else:
            keptBonds.add(bonds[0])
        if len(bonds) > 1:
            for i, bond1 in enumerate(bonds[1:]): # iterate over every detected bond except shortest
                angles = np.array([])
                for j, bond2 in enumerate(keptBonds): # compare to bonds we've already seen
#                     angles = np.append(angles, atomAnalysis.get_angle_value(image, [bond1, idx, bond2]))
                    angles = np.append(angles, atoms.get_angle(bond1, idx, bond2, mic = True))

                if np.all(angles > minAngle):
                    # keep if the angle of the new bond is large enough
                    # wrt the bonds we've decided to keep
                    keptBonds.add(bond1)
        coordinations[idx] = len(keptBonds)
        if coordinations[idx] < fullCoordinations[atom.symbol]:
            relativeCoordinations[idx] = -1
        elif coordinations[idx] > fullCoordinations[atom.symbol]:
            relativeCoordinations[idx] = 1
        else:
            relativeCoordinations[idx] = 0
        newBonds[idx] = keptBonds
    return relativeCoordinations, newBonds

def readStructs(datadir, slabEnergy, adsorbateEnergy):
    """
        Currently designed for output from single layer directory trees
        Reads in final adsorption geometries and energy data,
            returns dataframe with geometry and energy data

        Input:
            datadir: string that points to directory containing the following:
                - convergence: each line i has convergence status of run i
                - energies: each line i has final total energy from run i
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
    convergence = pd.read_csv(datadir + "convergence", header = None)
    energies = pd.read_csv(datadir + "energies", header = None)
    output =  pd.concat([energies, convergence], axis = 1)
    output.columns = ["E", "conv"]

    for i in os.listdir(datadir):
        key = re.search(r"output(\d+).gen", i)
        if key:
            key = int(key.group(1))
            geometries[key] =  gen.read_gen(datadir + i)
    output['geom'] = pd.Series(geometries)

    output = output[output['conv'] == "Geometry converged"]
    output = output.drop("conv", axis = 1)

    output['E_ads'] = output["E"] - (adsorbateEnergy + slabEnergy)
    return output


def convertAdsorbateToHe(struct, centerIndex, molIndices, height = None):
    """
    Preprocess final relaxed adsorption structures; replace adsorbate with He
    Input:
        struct: total structure (Atoms object)
        centerIndex: index of central atom (where He will be) (int)
        molIndices: list of indices to delete from the slab
        height(float) : height of He to be placed
    Output:
        output: Atoms object with He representing the location of the adsorbate
    """
    x, y, z = struct[centerIndex].position
    output = struct.copy()
    del output[[atom.index for atom in output if atom.index in molIndices]]
    if height:
        add_adsorbate(output, "He", height = height, position = (x, y))
    else:
        output.append(Atom("He", position=[x,y,z])) # adds to exact position of centeratom
    return output


def getSOAPs(geometries, rcut = 5, nmax = 10, lmax = 9, sigma = 0.1,
             periodic = True, crossover = True, sparse = False):
    """
    Takes a Series of geometries with one He present,
        returns SOAP representation of the chemical environment of He for each item
    Assumes any given structure in `geometries` has the same collection of elements
        as all the other structures
    Assumes any given structure in `geometries` has the same number of atoms as all
        the other structures

    Input:
        geometries: Series of Atoms objects; each must contain exactly 1 He atom
        rcut, nmax, lmax, sigma, periodic, crossover, sparse: SOAP parameters
    Output:
        output: Series of SOAP matrices, each corresponding to the appropriate index
    """
    refgeom = geometries.iloc[0] #use the first geometry as a reference geometry

    ## set up descriptor
    species = np.unique([i.symbol for i in refgeom])
    desc = SOAP(species=species, rcut = rcut, nmax = nmax, lmax = lmax,
                sigma = sigma, periodic = periodic, crossover = crossover, sparse = sparse)
    ## apply descriptor
    soaps = {}
    HeLoc = len(refgeom) - 1  # assume He atom is last one in Atoms list
    for i, geom in geometries.iteritems():
        tempSOAP = preprocessing.normalize(
            desc.create(geom, positions = [HeLoc], n_jobs = 4)) # SOAP representation of temp
        soaps[i] = tempSOAP[0]
    return pd.Series(soaps,name = 'SOAP')

# Python program to print connected
# components in an undirected graph
# https://www.geeksforgeeks.org/connected-components-in-an-undirected-graph/
 
 
class Graph:
 
    # init function to declare class variables
    def __init__(self, V):
        self.V = V
        self.adj = [[] for i in range(V)]
 
    def DFSUtil(self, temp, v, visited):
 
        # Mark the current vertex as visited
        visited[v] = True
 
        # Store the vertex to list
        temp.append(v)
 
        # Repeat for all vertices adjacent
        # to this vertex v
        for i in self.adj[v]:
            if visited[i] == False:
 
                # Update the list
                temp = self.DFSUtil(temp, i, visited)
        return temp
 
    # method to add an undirected edge
    def addEdge(self, v, w):
        self.adj[v].append(w)
        self.adj[w].append(v)
 
    # Method to retrieve connected components
    # in an undirected graph
    def connectedComponents(self):
        visited = []
        cc = []
        for i in range(self.V):
            visited.append(False)
        for v in range(self.V):
            if visited[v] == False:
                temp = []
                cc.append(self.DFSUtil(temp, v, visited))
        return cc

# This code is contributed by Abhishek Valsan

##############
# structures #
##############
mef = vasp.read_vasp("reference_files/CONTCAR_mef")
cf4 = vasp.read_vasp("reference_files/CONTCAR_cf4")
amorphous = gen.read_gen("reference_files/amorphous_base.gen")
xtl_n = vasp.read_vasp("reference_files/CONTCAR_nrich")
xtl_si = vasp.read_vasp("reference_files/CONTCAR_sirich")
xtl2x2 = gen.read_gen("reference_files/2x2xtl.gen")
xtl2x2_sifterm = gen.read_gen("reference_files/2x2xtl_sifterm.gen")
heavy_bomb = vasp.read_vasp("reference_files/CONTCAR_heavy_bombard")
bulk222 = vasp.read_vasp("reference_files/CONTCAR_222bulk")
annealed = vasp.read_vasp("reference_files/CONTCAR_annealed_unitcell")

#############
# constants #
#############
amu2g = 1.66054e-24 #multiply by this
a2cm = 1e-8 # multiply by this
