from ase.geometry.analysis import Analysis
from ase.formula import Formula
from ase import io
import numpy as np 
import os
import matplotlib.pyplot as plt
import pandas as pd

def bondAnalysis(data, focusElement = "C", bondelems = ["C", "F", "H", "Si", "N"], verbose = False):
    """
    `data` should be a pd Series consisting of (structure id: Atoms object) pairs 
    Length of returned values reflects only # of atoms of focusElement that have at least one bond to an 
    atom in bondelems. 
    
    """
    analyses = {key:Analysis(value) for key, value in data.iteritems()}
    cbonds = {key: {
        i: a.get_bonds(focusElement, i)[0] for i in bondelems}
     for key, a in analyses.items()}
    cIdxs = {key: [atom.index for atom in value if atom.symbol == focusElement]
             for key, value in data.iteritems()}
    # construct cbonds 
    cbonds = {}
    for key, lst in cIdxs.items():
        _bonds = analyses[key].all_bonds[0]
        _struct = data[key]
        for idx in lst:
            mybonds = _bonds[idx]
            mybondDict = {}
            for bondelem in bondelems:
                mybondDict[bondelem] = sum(_struct[i].symbol == bondelem for i in mybonds)
            if np.sum(pd.Series(mybondDict)) == 0:
                if verbose:
                    print("no bonds between focusElement and bondelems detected")
            else:
                cbonds[(key, idx)] =  mybondDict
    cbonds = pd.DataFrame(cbonds).T
    # construct combos
    combos = {}
    combolists = {}
    for key, value in cbonds.iterrows():
        newkey = "".join([key*value for key,value in value.iteritems() if value > 0])
        newkey = Formula(newkey).format('hill')
        combos[newkey] = combos.get(newkey,0) + 1
        combolists[newkey] = combolists.get(newkey,[]) + [key]
    combos = pd.Series(combos)
    combolists = pd.Series(combolists)
    return cbonds, combos, combolists

##########################
## Set up analysis here ##
##########################

#read data
s = io.read("POSCAR")
data = pd.Series({"0-0": s})

#specify parameters here 
cbonds, combos, combolists = bondAnalysis(data, focusElement = "F", bondelems = ["N", "Si"])

#show output
print('Group lists:')
print(combolists)
print("count # F shared by NSi, assuming any exist (errors if none)")
print(len(combolists['NSi']))

