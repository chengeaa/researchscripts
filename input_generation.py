# imports 

# scipy
import numpy as np

#ase
from ase.io import gen, vasp
from ase.build import make_supercell
from ase.build import add_adsorbate


# function definitions

def gridVasp(adsorbate, slab, directory = 'tempout/', spacing = 0.2):
    """
    Takes in slab, adsorbate, and desired output directory.
    Writes numbered POSCAR files to output directory with adsorbate 
    placed in grid pattern

    Args:
        adsorbate: either path (str) for POS/CONTCAR or Atoms obj
        slab: either path (str) for POS/CONTCAR or Atoms obj
        directory: path (str) for desired output location
        spacing: spacing between grid points (eg, 0.2 is a 5x5 grid)
    Returns:
        None
    """
    if type(adsorbate) == str:
        adsorbate = vasp.read_vasp(adsorbate)
        
    if type(slab) == str:
        slab = vasp.read_vasp(slab)
    # TODO: make this work with path objects?  

    i = 0 # 0 index written POSCAR{i} files
    for _x in np.arange(0, 1 - spacing, spacing):
        for _y in np.arange(0, 1 - spacing, spacing):
            s = slab.copy()
            x, y, z = slab.cell.cartesian_positions([_x, _y, 0])
            add_adsorbate(s, adsorbate, height = 2, position = (x, y))
            s = Atoms(sorted(s, key = lambda x: x.symbol))
            s.cell = slab.cell
            vasp.write_vasp(directory + "POSCAR" + str(i), s, 
                            label = '{} on {}'.format(
                                adsorbate.get_chemical_formula(),
                                slab.get_chemical_formula()
                                ), 
                            vasp5=True)
            
            i += 1
def randomGrid(adsorbate, slab, h = 2, outputDir = "tempout/",
        numDirs = 10, runsPerDir = 17, shallow = False, kind = 'gen'):
    """
    Produces a collection of structures with adsorbate randomly placed on 
    given slab. 
    Writes (.gen or vasp) outputs to desired directory
    
    Args:
        adsorbate: Atoms obj of adsorbate
        slab: Atoms obj of slab
        h: height of adsorbate (above max slab position)
        outputDir: Path (str) for desired output location.
        numDirs: Number of batches. Defaults to 10.
        runsPerDir: Number of sims per batch. Defaults to 17.
        shallow: If all runs to be at one directory level. Defaults to False.
        kind: output desired, from (gen, vasp). Defaults to gen. 
    Returns:
        None
    """
    np.random.seed(429)

    for d in range(numDirs):
        for run in range(runsPerDir):
            s = slab.copy()

            # generate random positions with dummy z required (3)
            p = s.cell.cartesian_positions(np.random.random(3)) 

            # construct and write
            add_adsorbate(s, adsorbate, height =  h, position = p[:2])

            if kind == 'gen':
                outname = "input{}-{}.gen".format(d, run) if not shallow else "input{}.gen".format(
                        d * runsPerDir + run)
                gen.write_gen(outputDir + outname, s)
            elif kind == 'vasp':
                outname = "POSCAR{}-{}".format(d, run) if not shallow else "POSCAR{}".format(
                        d * runsPerDir + run)
                vasp.write_vasp(outputDir + outname, s,
                        sort = True, vasp5 = True
                        )
            else:
                raise AssertionError("kind should be from (vasp, gen)")

def randomGridMultiple(n, adsorbate, slab, h = 2, outputDir = "tempout/", 
        numDirs = 10, runsPerDir = 17, shallow = False, kind = 'gen', minDist = 1):
    """
    Produces a collection of structures (of size ``numDirs`` * ``runsperDir``;
    each has n adsorbates randomly placed on given slab. 
    Writes (.gen or vasp) outputs to desired directory
    
    Args:
        n: number of adsorbates per slab
        adsorbate: Atoms obj of adsorbate
        slab: Atoms obj of slab
        h: height of adsorbate (above max slab position)
        outputDir: Path (str) for desired output location.
        numDirs: Number of batches. Defaults to 10.
        runsPerDir: Number of sims per batch. Defaults to 17.
        shallow: If all runs to be at one directory level. Defaults to False.
        kind: output desired, from (gen, vasp). Defaults to gen. 
        minDist: minimum distance between placed adsorbates (in Å)
    Returns:
        None
    """

    np.random.seed(429)

    for d in range(numDirs):
        for run in range(runsPerDir):
            s = slab.copy()

            positions = []
            while len(positions) < n: 
                # generate random positions with dummy z required (3)
                pnew = s.cell.cartesian_positions(np.random.random(3)) 
                if positions:
                    dists = np.array([np.sqrt((pnew[0] - p[0])**2 + (pnew[1] - p[1])**2)
                        for p in positions])
                    if np.any(dists < minDist):
                        continue
                positions += [pnew]
                    
            # construct and write
            for p in positions:
                add_adsorbate(s, adsorbate, height =  h, position = p[:2])

            if kind == 'gen':
                outname = "input{}-{}.gen".format(d, run) if not shallow else "input{}.gen".format(
                        d * runsPerDir + run)
                gen.write_gen(outputDir + outname, s)
            elif kind == 'vasp':
                outname = "POSCAR{}-{}".format(d, run) if not shallow else "POSCAR{}".format(
                        d * runsPerDir + run)
                vasp.write_vasp(outputDir + outname, s,
                        sort = True, vasp5 = True
                        )
            else:
                raise AssertionError("kind should be from (vasp, gen)")
