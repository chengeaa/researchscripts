# imports 
# base python
import os, copy, re, random, pickle

# scipy
import numpy as np
import pandas as pd

#ase
from ase.io import gen, vasp, xyz, extxyz, dftb
from ase import Atoms, Atom
from ase.constraints import FixAtoms
from ase.visualize import view
from ase.build import make_supercell, add_adsorbate

#dscribe
from dscribe.descriptors import SOAP, MBTR
from dscribe.kernels import REMatchKernel, AverageKernel

from sklearn import preprocessing


def predictz(surf, x, y, zmodel = 'zmodel.pkl'):
    """
    surf: bare substrate
    x, y: position at which to place adsorbate
    zmodel: Path object pointing to the model pickle
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
    
    X = getSOAPs(pd.Series({0: surf}))[0].reshape(1, -1) #reshape because just one sample
    if zmodel: 
        with open('zmodel.pkl', 'rb') as f:
            zmodel = pickle.load(f)
        
        predz = zmodel.predict(X)
    else:
        # TODO: implement me hehe (someday, this will train a new model if one isn't provided?)
        predz = maxz
    return predz


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



def main(basename):
    """
    Perform ML-based isotherm seeding.

    Args:
        basename: name of base slab
    """

    # load z prediction and E_ads prediction models (pickled KRR models)
    with open('models/zmodel.pkl', 'rb') as f:
        zmodel = pickle.load(f)

    with open('models/Emodel.pkl', 'rb') as f:
        Emodel = pickle.load(f)


    # load base slab, remove extraneous atoms, and wrap
    base = gen.read_gen(basename) 
    del base[[atom.index for atom in base if atom.symbol in ['He', 'Ar']]]
    base.wrap()


    # generate regular grid based on cell parameters of slab
    a,b,c = base.cell
    a,b,c = np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c)
    npoints = 20
    apoints = np.linspace(0, a, npoints) # regular spacing
    bpoints = np.linspace(0, b, npoints) # regular spacing

    # place He atoms in grid points
    gridpoints = []
    for apoint in apoints:
        for bpoint in bpoints:
            newstruct = base.copy()
            zhat = predictz(newstruct, apoint, bpoint)
            newstruct.append(Atom('He', position = (apoint, bpoint, zhat)))
            gridpoints += [newstruct]


    # generate pd df with data
    gridpoints = pd.Series(gridpoints)
    gridpoints = pd.DataFrame({'geom': gridpoints})
    gridpoints = pd.concat([gridpoints, getSOAPs(gridpoints['geom'])], axis = 1)

    # data matrix for ML
    X = pd.DataFrame(gridpoints['SOAP'].to_list(), index = gridpoints.index)

    gridpoints['predE'] = Emodel.predict(X)

    charges = np.append(np.zeros(len(base)), gridpoints['predE'])
    base.set_initial_charges(charges)
    for geom in gridpoints['geom']:
        base.append(Atom("He", position = geom[-1].position))

    # TODO adaptive sampling portion 
    if visualize:
        view(visbase)

        print("pearson r:", 
        pearsonr([geom[-1].position[2] for geom in gridpoints['geom']],
                 gridpoints['predE']
                )
        )


if __name__ == "__main__":
    """
    Takes in  arguments:
    1. relative path for the base substrate
    """
    args = sys.argv[1:]
    if len(args) > 3:
        print(args)
        raise Exception("No more than 3 arguments allowed")
    main(*args)



## utilities below 

