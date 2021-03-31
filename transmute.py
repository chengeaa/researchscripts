from ase.geometry.analysis import Analysis
from utils import readStructs
from analysis import coordLabeller

def transmute(threshold, 
         transmute = True,
         target = "C",
         slabElems = ["Si", "N"], 
         surfDepth = 6,
#          simple = True,
         numStructs = -1,
         numOut = -1,
         **kwargs):
    """
        Find undercoordinated surface/near-surface atoms and transmute them to C
        Inputs:
            threshold (k): maxiumum number of bonds to be considered "undercoordinated" (should be 1 or 2)
            transmute: whether to transmute elements or not 
            slabElems: elements to consider as part of the 'core slab'; only these will be transmuted
            surfDepth: depth of atoms to consider below surface; surface defined as 
                max z coord among slabElems
            simple: Use basic bond counting as implemented in ASE, 
                or advanced with coordLabeller from `dependencies` (simple=False is UNTESTED!!!)
            numStructs: Number of structures to consider
            numOut: Number of structures to output
            **kwargs: to be passed to the readStructs function from `dependencies`
    """
    simple = True # only support simple bond count for now 
    if numStructs > 0:
        data = readStructs(**kwargs)['geom'][:numStructs]
    else:
        data = readStructs(**kwargs)['geom']
    
    ## make bond counts
    
    if simple:
        analyses = {key: Analysis(value) for key, value in data.items()}
        bonds = {key: value.all_bonds[0] for key, value in analyses.items()}
        numBonds = {key: [len(i) for i in value] for key, value in bonds.items()}
        for key, b in bonds.items():
            for idx, atomBonds in enumerate(b):
                for bond in atomBonds:
                    e1, e2 = data[key][idx].symbol, data[key][bond].symbol
                    if e1 not in slabElems or e2 not in slabElems:
                        numBonds[key][idx] -= 1
            
    else:
        bonds = {key: coordLabeller(value, 0, angle_tolerance = .25, bond_tolerance = .15, minz = minz)[1]
             for key, value in data.items()}
        numBonds = {key: [len(value[i]) for i in range(len(value))] for key, value in bonds.items()}
        for key, b in bonds.items():
            for idx, bond in b.items(): #dictionary
                for bond in atomBonds:
                    e1, e2 = data[key][idx].symbol, data[key][bond].symbol
                    if e1 not in slabElems or e2 not in slabElems:
                        numBonds[key][idx] -= 1
    
    ## look for criteria: surface depth, core slab elements, low coordination
    toReplace = {}
    for key, numBondList in numBonds.items():
        minz = np.max([atom.position[2] for atom in data[key] if atom.symbol in slabElems]) - surfDepth
        toReplace[key] = [int(nBonds <= threshold and 
                              data[key][i].symbol in slabElems and
                              data[key][i].position[2] > minz
                             ) 
                         for i, nBonds in enumerate(numBondList)]
    
    ## set labels for atoms to be transmuted, and generate counts
    numTransmuted = {}
    for key, value in data.items():
        data[key].set_tags(toReplace[key])
        data[key].wrap()
        numTransmuted[key] = np.sum(toReplace[key])
    data = pd.DataFrame(data)
    data['n'] = pd.Series(numTransmuted)
    data = data.sort_values("n", ascending = False)
    if numOut > 0:
        data = data.iloc[:numOut, :]
    
    ## transmute atoms
    if transmute:
        for key in data.index:
            for atom in data['geom'][key]:
                if toReplace[key][atom.index]:
                    atom.symbol = target
        
    return data

def write_transmuted(outdir, numout, numrep):
    """
        Write results of transmute function into some directory
        numout corresponds the top `numout` structures to output
        numrep is the number of times to replicate each of the above structures,
            with random adjustment of last atom position (should be Ar for bomb)
        files written in input$i-$j.gen format
    """
    np.random.seed(429)
    for i, key in enumerate(output.index[:10]):
        for j in range(17):
            temp = output.loc[key, 'geom'].copy()
            print(len(temp))
            x, y, _ = temp.cell.cartesian_positions(np.random.random(size = 3))
            temp[-1].position = [x, y, temp[-1].position[2]]
            print(len(temp))
            gen.write_gen(outdir + "input{}-{}.gen".format(i,j), temp)

if __name__ == "__main__":
    """
    Direct transumation of atoms based on coordination
    Arg 1: threshold (maximum number of NNs to select for transmutation)
    Arg 2-n: see main docstring
    """
    args = sys.argv[1:]
    
    if len(args) < 1:
        print(args)
        raise Exception("Too few args, need at least the NN threshold")
    transmute(*args)

