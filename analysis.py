from ase.neighborlist import NewPrimitiveNeighborList, natural_cutoffs
from researchscripts.utils import readStructs, KE, vAngle
from researchscripts.structure import Graph, getFragIndices
from ase.geometry.analysis import Analysis
from researchscripts.utils import readStructs
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd


def coordLabeller(atoms, fullCoordinations = {"Si": 4, "N":3, "H":1, "F": 1},
                  minAngles = {"Si": 90, "N": 109.5, "H": 360, "F": 360}, # 360 for atoms with max 1 bond
                  maxBonds_per_element = {"Si": 6, "N":4, "H":1, "F":1},
                  angle_tolerance = 0.15, #tolerance for valid bond angles
                  bond_tolerance = 0.25, #tolerance for valid bond angles
                  minz = 0 #minimum height above which to compute
                 ):
    """
    Takes a structure, returns two dictionaries, the keys of which are identical: 
    the index of the atom for which the statistic is calculated
    relativeCoordinations: -1 if atom i is undercoordinated, 1 if overcoordinated, 0 else
    bonds: list of bonds for atom i
    Defaults for angle and bond tolerance based on amorphous SiN
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


def analyzeFragments(datadir, **kwargs):
    """
        Pass in `name` and `shallow` kwargs if needed for utils.readStruct function
    """
    geometries = readStructs(datadir, **kwargs)['geom']
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
    
    # totalbonds = []
    # bondcounts = {}
    # for key, analysis in analyses.items(): 
        # try:
            # totalbonds += [len(analysis.get_bonds(e1, e2, unique = True)[0])]
            # bondcounts[key] =  len(analysis.get_bonds(e1, e2, unique = True)[0])
        # except:
            # print('error on {}'.format(key))
    # totalbonds = np.array(totalbonds)
    # if form:
        # print('percent runs with {}-{} bond formation = {}'.format(e1, e2, np.sum(totalbonds > 0)/170))
    # else:
        # print('average number of final {}-{} bonds = {}'.format(e1, e2, np.sum(totalbonds)/170))
    # # plt.hist(totalbonds, bins = np.arange(0, np.max(totalbonds) + 1))
    # # plt.hist(totalbonds, bins = np.arange(5, 14))
    # if form:
        # plt.title('distribution of # of {}-{} bonds formed'.format(e1, e2));
    # else:
        # plt.title('distribution of # of {}-{} bonds count'.format(e1, e2));
    # plt.show()   

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



def plotBondDistributionComposition(data, a='Si', b='N'):
    """
    Requires ``data``, a pd df with 'geom', 'coordlabels', and (optionally) 'wantedIndices'
    Gives overall distribution and percent makeup of elements A and B in that distribution
    """
    if 'wantedIndices' in data.columns:
        sibonds = pd.Series({idx:
            [(
                np.array([data['geom'][idx][i].symbol for i in value])
            ) for key, value in
             pd.Series(data['coordlabels'][idx][1])[data['wantedIndices'][idx]].items()
            ]
         for idx in data.index
        })
    else:
        sibonds = pd.Series({idx:
            [(
                np.array([data['geom'][idx][i].symbol for i in value])
            ) for key, value in
             pd.Series(data['coordlabels'][idx][1]).items()
            ]
         for idx in data.index
        })

    nfrac = np.zeros(7)
    sifrac = np.zeros(7)
    ndatapoints = np.zeros(7)

    for key, item in sibonds.items():
        for i in item:
            idx = len(i)
            ndatapoints[idx] += 1
            nfrac[idx] += np.sum(np.array(i) == a)/idx
            sifrac[idx] += np.sum(np.array(i) == b)/idx
    ndatapoints = np.array([i if i > 0 else 1 for i in ndatapoints])
    nfrac /= ndatapoints
    sifrac /= ndatapoints

    mask = np.array([i > 10 for i in ndatapoints])
    nfrac, sifrac, ndatapoints = nfrac[mask], sifrac[mask], ndatapoints[mask]


    plt.figure(figsize = (10,5))
    plt.hist(
        [i for sublist in
            pd.Series({idx:
                [len(
                    np.array([data['geom'][idx][i].symbol for i in value])
                ) for key, value in
                 pd.Series(data['coordlabels'][idx][1])[data['wantedIndices'][idx]].items()
                ]
             for idx in data.index
            })
         for i in sublist], bins = np.arange(7), density = True, label = 'total density'
    );
    plt.plot(np.arange(7)[mask], nfrac, marker = 'o', label = '{} fraction'.format(a))
    plt.plot(np.arange(7)[mask], sifrac, marker = 'o', label = '{} fraction'.format(b))
    plt.title("Normed hist of bond counts with {}/{} fraction".format(a,b))
    plt.legend()


def getBondSubset(data, elem):
    """
    Requires ``data`` pd df that has columns 'geom', 'coordlabels', 'wantedIndices'
    Returns pd Series with the subsetted bond counts of (A, elem) where A should be set by the input wantedIndices
    """
    result = pd.Series({idx: 
        [np.sum(
            np.array([data['geom'][idx][i].symbol for i in value]) == elem
        ) for key, value in 
         pd.Series(data['coordlabels'][idx][1])[data['wantedIndices'][idx]].items()
        ]
     for idx in data.index
    })
    return result

def getFragStats(traj, fragments, indices, verbose = False):
    trajVelos = np.array([i.arrays['vel'] for i in traj]) # Read velocities; shape is (T, N, 3)
    # Get the velocities of each fragment at each time point
    fragVelosByPoint = {}
    for key, frags in fragments.items():
        pointVelos = []
        for frag in frags:
            velo = trajVelos[indices[key], frag, :]
            pointVelos += [velo]
        fragVelosByPoint[key] = pointVelos

    # Collect the relevant statistics: KE, p, angles
    fragAnglesByPoint = {}
    fragPositionsByPoint = {}
    fragMomentaByPoint = {}
    fragKEByPoint = {}
    fragLabelByPoint = {}
#     print(fragVelosByPoint)
#     print(fragsByPoint)
    for key, velos in fragVelosByPoint.items(): # `velos` at each point is a list of velocities for each fragment
        pointAngles = []
        pointMomenta = []
        pointKE = []
        pointLabel = []
        pointPos = []
        for velo, frag in zip(velos, fragments[key]):
            _velo = velo/1000 # convert to Å/fs

            # Get label
            label = traj[0][frag].symbols.get_chemical_formula()

            # Get average fragment position (x,y,z)
            _positions = traj[indices[key]][frag].positions
            position = np.mean(_positions, axis = 0)

            # Masses should be the same for all timepoints within a trajectory,
            # but I'll just use the Point for temporary convenience
            masses = np.array([atom.mass for atom in traj[indices[key]][frag]])

            # Calculate average velocities by component and by atom
            _fragVelos = np.mean(_velo, axis = 0) # average velocity across each dimension
            _atomVelos = np.sqrt(np.sum(_velo**2, axis = 1)) #v_i = \sqrt{vx_i**2 + vy_i**2 + vz_i**2}

            # Calculate kinetic energies (per atom, per fragment)
            atomKEs = KE(_atomVelos, masses) #<1/2 * \sum m_i(vx^2+vy^2+vz^2)^(1/2)
            avgKE = np.mean(atomKEs)

            # Calculate momenta (per atom, per fragment) - magnitude only
            atomMomenta = np.abs(_atomVelos * masses)
            avgMomentum = np.mean(atomMomenta)

            # Calculate fragment angle
            angle = vAngle(_fragVelos)
            if verbose:
                print('point, fragment: {} (index {}) - {} (indices: {})'.format(
                    key, indices[key], label, frag))
                print('velocities: \n{}'.format(_velo))
                print('average velocities (componentwise): {}'.format(_fragVelos))
                print('average velocities (atomwise): {}'.format(_atomVelos))
                print("T_i = {}".format(atomKEs))
                print("<T> = {:.2E}".format(avgKE))
                print("p_i: {}".format(atomMomenta))
                print("<p>: {}".format(avgMomentum))
                print("theta = {:.2f}".format(angle))
                print()
            pointAngles += [angle]
            pointMomenta += [avgMomentum]
            pointKE += [avgKE]
            pointLabel += [label]
            pointPos += [position]
        fragAnglesByPoint[key] = pointAngles
        fragPositionsByPoint[key] = pointPos
        fragMomentaByPoint[key] = pointMomenta
        fragKEByPoint[key] = pointKE
    dfPoint = pd.DataFrame({
        "timeIndex" : indices,
        "fragIndices" : fragments,
        "fragPositions": fragPositionsByPoint,
        "fragVelos": fragVelosByPoint,
        "fragLabels": fragLabelByPoint,
        "angle": fragAnglesByPoint,
        "KE": fragKEByPoint,
        "momentum": fragMomentaByPoint
    })

    dfTraj = None
    return dfPoint, dfTraj

def getFragsTraj(traj, exclusions = ['CF4', 'Ar', 'CH3F']):
    """
    Returns two lists of lists:
    - list of lists, each list is a fragment
    - list of labels corresponding to each fragment
    """
    fragIndices = []
    fragLabels = []
    for frame in traj:
        frags = getFragIndices(frame, True)
        labels = [frame[frag].symbols.get_chemical_formula() for frag in frags]
        fragIndices += [[i for i, j in zip(frags, labels) if j not in exclusions]]
        fragLabels += [[j for i, j in zip(frags, labels) if j not in exclusions]]
    return fragIndices, fragLabels


def getFragsByPoint(traj, points, ArIndex = -1, veloConversionFactor = 1/1000, exclusions = ['CF4', 'Ar', 'CH3F'], 
        KEcutoff = 0.05, numPointstoCheck = 10, keyPoints = ["start", "maxfrag", "end"], verbose = False):
    """
    Produces a list of lists for each point, each parent being list of fragments and child being list of
    atom indices
    """
    pointIndices = {}
    trajVelos = np.array([i.arrays['vel'] for i in traj]) # Read velocities; shape is (T, N, 3)

    # Get argon KEs over traj
    arVelos = [i[ArIndex] for i in trajVelos]
    arVeloTotals = [np.sqrt(np.sum(i**2)) for i in arVelos]
    arKEs = [KE(i * veloConversionFactor) for i in arVeloTotals]

    # Determine indices for each key point
    pointIndices['start'] = 0
    pointIndices['end'] = -1

    # Check for index with maximum fragmentation post-collision
    checkedPoints = 0
    maxFrags = 0
    maxFragIndex = 0

    for i, arKE in enumerate(arKEs):
        if arKE < arKEs[0] * KEcutoff: # Fragment count begins when Ar KE drops below 0.05 of starting KE
            if checkedPoints < numPointstoCheck:
                checkedPoints += 1

                # This is kind of inefficient but we're gonna do it anyway
                fragIndices = getFragIndices(traj[i])
                frags = [f for f in fragIndices if
                                   traj[i][f].symbols.get_chemical_formula() not in exclusions]
                numFrags = len(frags)

                if verbose:
                    print('checking point ', i)
                    print('fragment indices: ', frags)
                    print('fragment labels: ', [traj[i][f].symbols.get_chemical_formula() for f in frags])
                    print('number of fragments: ', numFrags)

                # Stop or update
                if numFrags < maxFrags:
                    break # Stop searching once recombination begins happening
                if numFrags > maxFrags:
                    maxFragIndex = i
                    maxFrags = numFrags

    # Update maxfrag entry after checking the postcollision points for max fragmentation

    pointIndices['maxfrag'] = maxFragIndex

    # Check that all keyPoints have an index
    for key in keyPoints:
        if key not in pointIndices.keys():
            raise ValueError("Not all keys in keyPoints have corresponding value in pointIndices")

    # Obtain fragment indices
    fragsByPoint = {}
    for P in keyPoints:
        pointIndex = pointIndices[P]
        fragIndices = getFragIndices(traj[pointIndex])

        fragsByPoint[P] = [f for f in fragIndices if
                           traj[pointIndex][f].symbols.get_chemical_formula() not in exclusions]
    if verbose:
        print("indices considered: {}".format(pointIndices))
        print("fragment indices:")
        for key, item in fragsByPoint.items():
            print(key, item)
    return pointIndices, fragsByPoint
