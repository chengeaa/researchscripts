from dscribe.descriptors import SOAP

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

def convertAdsorbateToHe(struct, centerIndex, molIndices, height = None):
    """
    Preprocess final relaxed adsorption structures; replace adsorbate with He

    Args:
        - struct: total structure (Atoms object)
        - centerIndex: index of central atom (where He will be) (int)
        - molIndices: list of indices to delete from the slab
        - height(float) : height of He to be placed
    Returns:
        - output: Atoms object with He representing the location of the adsorbate
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
    Assumes any given structure in ``geometries`` has the same collection of elements
        as all the other structures
    Assumes any given structure in ``geometries`` has the same number of atoms as all
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

