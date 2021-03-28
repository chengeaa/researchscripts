#!/usr/bin/env python3
'''
Designed to work on the cluster, removing 'ejected species' after each equilibration step
'''
from dependencies import *
import time

def main(cutoff, 
        datapath, 
        zmodelPath, 
        EmodelPath, 
        smallset,
        adsorbate = 'mef',
        npoints = 20, 
        outputpath = 'input.gen'):
    last = time.time()
    npoints = int(npoints)
    cutoff = float(cutoff)
    smallset = bool(int(smallset))
    adsorbate_types = {'mef': mef, 'cf4': cf4}
    ads = adsorbate_types[adsorbate]

    with open(zmodelPath, 'rb') as f:
        zmodel = pickle.load(f)
    with open(EmodelPath, 'rb') as f:
        Emodel = pickle.load(f)

    # read in calculated structure
    if "gen" in datapath:
        data = gen.read_gen(datapath)
    elif "CAR" in datapath:
        data = vasp.read_vasp(datapath)

    print('maxz: ', max([i.position[2] for i in data]))
    print('data read')
    now = time.time()
    print(now - last)
    last = now

    # obtain base slab
    base = getslab(data)
    # assume any adsorption influence enters via config, independent of Ar
    del base[[atom.index for atom in base if atom.symbol in ['He', 'Ar']]]
    base.wrap()

    print('base obtained')
    now = time.time()
    print(now - last)
    last = now
    # set up gridpoints with predicted z heights

    a,b,c = base.cell
    a,b,c = np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c)
    apoints = np.linspace(0, a, npoints)
    bpoints = np.linspace(0, b, npoints)

    if smallset:
        species = ['Si', 'N', 'H', 'He']
    else:
        species = ["Si", "N", "H", "C", "F", "Ar", "He"]

    print(smallset, species)

    gridpoints = []
    for apoint in apoints:
        for bpoint in bpoints:
            newstruct = base.copy()
            print(newstruct)
            zhat = predictz(newstruct, apoint, bpoint, zmodel, species)
            newstruct.append(Atom('He', position = (apoint, bpoint, zhat)))
            gridpoints += [newstruct]

    print('gridpoints done')
    now = time.time()
    print(now - last)
    last = now

    gridpoints = pd.Series(gridpoints)
    gridpoints = pd.DataFrame({'geom': gridpoints})

    # add SOAP representation for gridpoint structs
    gridpoints = pd.concat([gridpoints, getSOAPs(gridpoints['geom'],
        species = species
        )], axis = 1)
    
    # create prediction matrix
    X = pd.DataFrame(gridpoints['SOAP'].to_list(), index = gridpoints.index)

    # predict energies, append to original df
    gridpoints['predE'] = Emodel.predict(X)

    # create 'visbase': struct with all He points included in one struct
    charges = np.append(np.zeros(len(base)), gridpoints['predE'])
    visbase = base.copy()
    for geom in gridpoints['geom']:
        visbase.append(Atom("He", position = geom[-1].position)) 
    visbase.set_initial_charges(charges)

    view(visbase)

    print('energy prediction done')
    now = time.time()
    print(now - last)
    last = now
    
    # assess gridpoints and place adsorbates
    gridpoints = gridpoints.sort_values(by = 'predE')
    gridpoints['xpos'] = [geom[-1].position[0] for geom in gridpoints['geom']]
    gridpoints['ypos'] = [geom[-1].position[1] for geom in gridpoints['geom']]
    gridpoints['zpos'] = [geom[-1].position[2] for geom in gridpoints['geom']]

    adsorbatePoints = []
    a = visbase.cell[0]
    b = visbase.cell[1]
    for _, row in gridpoints.iterrows():
        isclose = False
        point1 = np.array([row['xpos'], row['ypos']])
        for x, y, z in adsorbatePoints:
            for dispx in [-a, a*0, a]:
                for dispy in [-b, b*0, b]:
                    point2 = np.array([x, y])
                    point2 = point2 + dispx[:2] + dispy[:2]

                    if np.linalg.norm(point1 - point2) < cutoff:
                        isclose = True
        if not isclose:
            adsorbatePoints.append(np.append(point1, row['zpos']))

    print('placement done')
    now = time.time()
    print(now - last)
    last = now

    adsvisbase = base.copy()
    maxz = np.max([atom.position[2] for atom in adsvisbase])

    for point in adsorbatePoints:
        print(point[2])
        add_adsorbate(adsvisbase, ads, height = point[2] - maxz + 1, position = (point[0], point[1]))

    gen.write_gen(outputpath, adsvisbase)

    view([data,base, adsvisbase])
if __name__ == "__main__":
    """
    Takes in some structure with some stuff on surface maybe
    Removes that stuff
    Adds the molecules that I want adsorbed according to the predicted isotherm
    """

    args = sys.argv[1:]
    nargs = 8
    if len(args) > nargs:
        print(args)
        raise Exception("No more than {} arguments allowed".format(nargs))
    main(*args)
