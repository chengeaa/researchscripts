from ase.neighborlist import NewPrimitiveNeighborList, natural_cutoffs
from utils import readStructs
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


def analyzeFragments(datadir, **kwargs):
    """
        Pass in `name` and `shallow` kwargs if needed for utils.readStruct function
    """
    geometries = readStructs(datadir, **kwargs)
    analyses = {key: Analysis(item) for key, item in geometries.items()}
    analyses = pd.Series(analyses)

    #####################
    ### fragmentation ###
    #####################

    fragmentLists = []
    for struct in geometries:
        adjmat = Analysis(struct).adjacency_matrix[0]
        numnodes = adjmat.shape[0]
        g = Graph(numnodes)
        for i in range(numnodes):
            for j in range(numnodes):
                if adjmat[i,j]:
                    g.addEdge(i,j)
        cc = g.connectedComponents()
        isSmallgraph = np.array([len(i) for i in cc]) < 10
        smallgraphs = []
        for i, subgraph in enumerate(cc):
            if isSmallgraph[i]:
                smallgraphs += [struct[[atom.index for atom in struct if atom.index in subgraph]]]
        fragmentLists += [smallgraphs]

    flatten = lambda t: [item for sublist in t for item in sublist]
    fragmentTypes = np.unique([i.symbols.get_chemical_formula() for i in flatten(fragmentLists)])
    fragdict = {i:j for i, j in zip(geometries.keys(), fragmentLists)}

    fragmentData = pd.DataFrame({key: [0] * len(geometries) for key in fragmentTypes})

    fragmentData.index = fragdict.keys()

    for key, fragmentList in fragdict.items():
        for fragment in fragmentList:
            _symbol = fragment.symbols.get_chemical_formula()
            fragmentData[_symbol].loc[key] += 1
    fragmentData.to_csv(datadir + "fragdata.csv")
    print(fragmentData.sum(axis = 0))


    ###################  
    ### bond counts ###
    ###################
    
    totalbonds = []
    bondcounts = {}
    e1 = "Si"
    e2 = "N"
    form = False
    for key, analysis in analyses.items(): 
        try:
            totalbonds += [len(analysis.get_bonds(e1, e2, unique = True)[0])]
            bondcounts[key] =  len(analysis.get_bonds(e1, e2, unique = True)[0])
        except:
            print('error on {}'.format(key))
    totalbonds = np.array(totalbonds)
    if form:
        print('percent runs with {}-{} bond formation = {}'.format(e1, e2, np.sum(totalbonds > 0)/170))
    else:
        print('average number of final {}-{} bonds = {}'.format(e1, e2, np.sum(totalbonds)/170))
    # plt.hist(totalbonds, bins = np.arange(0, np.max(totalbonds) + 1))
    # plt.hist(totalbonds, bins = np.arange(5, 14))
    if form:
        plt.title('distribution of # of {}-{} bonds formed'.format(e1, e2));
    else:
        plt.title('distribution of # of {}-{} bonds count'.format(e1, e2));
    plt.show()   

    from itertools import combinations
    elems = ["Si", "F", "N", "C", "H", "Ar"]

    result = {}
    for e1, e2 in combinations(elems, 2):
        bondcounts = {}
        for key, analysis in analyses.items(): 
            bondcounts[key] =  len(analysis.get_bonds(e1, e2, unique = True)[0])
        bondcounts = pd.Series(bondcounts)
        result["{}-{}".format(e1, e2)] = bondcounts
    pd.DataFrame(result).to_csv(datadir+"bondcounts.csv")
fragmentData