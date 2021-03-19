###########
# imports # 
###########
# base python
import os
import pickle
import copy
import re
from sys import getsizeof
import sys
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
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.kernel_ridge import KernelRidge
from sklearn import preprocessing
import matplotlib.tri as tri

#plt.rcParams["figure.figsize"] = (20,7)

#ase
from ase.io import gen, vasp, xyz, extxyz, dftb
from ase.io.dftb import read_dftb_velocities, write_dftb_velocities
from ase.calculators.dftb import Dftb
from ase import Atoms, Atom
from ase.constraints import FixAtoms
from ase.formula import Formula
from ase.visualize import view
from ase.build import make_supercell
from ase.visualize.plot import plot_atoms
from ase.build import add_adsorbate, add_vacuum
from ase.cell import Cell
from ase.neighborlist import NewPrimitiveNeighborList, natural_cutoffs
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

def postprocessResults(directory = "../"):
    """
        Takes in a list of indices, corresponding to the bombardment trials to analyze
        Looks for files named `results$i{bomb,quench,eq}.csv` in directory specified. 
        Returns list of 3 dfs; each one has elements and keys
    """

    subdirs = np.arange(10)
    bombdata = {i :pd.read_csv(directory + "results%dbomb.csv" % i, index_col=0) 
            for i in subdirs}
    quenchdata = {i : pd.read_csv(directory + "results%dquench.csv" % i, index_col=0) 
            for i in subdirs}
    eqdata = { i: pd.read_csv(directory + "results%deq.csv" % i, index_col=0) 
            for i in subdirs}

    return [bombdata, quenchdata, eqdata]



def postprocessAggregated(simindices, directory = "../"):
    """
        Takes in a list of indices, corresponding to the bombardment trials to analyze
        Looks for files named `aggregated_{bomb,quench,eq}$i` in directory specified. 
    """
    bombdata = {i :pd.read_csv(directory + "aggregated_bomb%d.csv" % i, index_col=0) 
            for i in simindices}
    quenchdata = {i : pd.read_csv(directory + "aggregated_quench%d.csv" % i, index_col=0) 
            for i in simindices}
    eqdata = { i: pd.read_csv(directory + "aggregated_eq%d.csv" % i, index_col=0) 
            for i in simindices}
    data = {"bomb": bombdata, "quench": quenchdata, "eq":eqdata}


    aggregated = {}
    for i in simindices:
        for step in ["bomb", "quench", "eq"]:
            aggregated["%i-%s" % (i, step)] = data[step][i].sum(axis = 1)
    return pd.DataFrame(aggregated)


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


def getSOAPs(geometries, species,
        rcut = 5, nmax = 10, lmax = 9, sigma = 0.1,
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
#   refgeom = geometries.iloc[0] #use the first geometry as a reference geometry

    ## set up descriptor
#   species = np.unique([i.symbol for i in refgeom])
    desc = SOAP(species=species, rcut = rcut, nmax = nmax, lmax = lmax,
                sigma = sigma, periodic = periodic, crossover = crossover, sparse = sparse)
    ## apply descriptor
    soaps = {}
    for i, geom in geometries.iteritems():
        HeLoc = len(geom) - 1  # assume He atom is last one in Atoms list
        tempSOAP = preprocessing.normalize(
            desc.create(geom, positions = [HeLoc], n_jobs = 4)) # SOAP representation of temp
        soaps[i] = tempSOAP[0]
    return pd.Series(soaps,name = 'SOAP')


def predictz(surf, x, y, zmodel, species):
    """
    surf: bare substrate
    x, y: position at which to place adsorbate
    zmodel: model object

    returns predicted z value
    """
    searchR = 2.2
    surf = surf.copy()
    add_adsorbate(surf, 'He', height = 0, position = (x, y))

    maxz = 0
    for atom in surf:
        if atom.symbol == "He": # don't use He position to determine max Z position
            continue
        _x, _y, _z = atom.position
        if ((x - _x)**2 + (y - _y)**2) ** 0.5 < searchR:
            if _z > maxz:
                maxz = _z + 2.5

    surf[-1].position[2] = maxz

    X = getSOAPs(pd.Series({0: surf}), species = species)[0].reshape(1, -1) #reshape because just one sample
    print(X.shape)
    if zmodel:
        predz = zmodel.predict(X)
#       print(maxz, predz)
    else:
#         TODO: implement me hehe
        predz = maxz + 2.5
    return predz

def getslab(struct):
    """
    Input: 
        struct: structre from which we will trim unbound species (.gen file)
    Output:
        baseslab: structure with unbound species trimmed
    """
    adjmat = Analysis(struct).adjacency_matrix[0]
    numnodes = adjmat.shape[0]
    g = Graph(numnodes)
    for i in range(numnodes):
        for j in range(numnodes):
            if adjmat[i,j]:
                g.addEdge(i,j)
    cc = g.connectedComponents()
    maingraph = np.array([i for i in cc if 0 in i][0])
    return struct[[atom.index for atom in struct if atom.index in maingraph]]
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
