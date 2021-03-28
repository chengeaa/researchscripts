from ase.neighborlist import NewPrimitiveNeighborList, natural_cutoffs
import numpy as np


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
import numpy as np

