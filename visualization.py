import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms

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

def plotElemDist(data, targetElem = "C", latticeElems = ["Si", "N", "H"], nbins = 25, stacked = False):
    """
    Plot distribution of element within slab, data should be arraylike collection of stuctures    
    """
    targetZs = []
    latticeZs = []

    # populate a cZs list of hists, latticeZs list of hists
    for key, value in data.iteritems():
        targetZs += [atom.position[2] for atom in value if atom.symbol == targetElem]
        latticeZs += [atom.position[2] for atom in value if atom.symbol in latticeElems]


    minZ, maxZ = np.min(latticeZs), np.max(latticeZs)
    bins = np.linspace(minZ, maxZ, nbins)
    width = (maxZ-minZ)/nbins

    if stacked:
        h = plt.hist([targetZs, latticeZs], bins = bins, density = True, alpha = 1, 
                 label = "stacked {} and {} distributions".format(targetElem, latticeElems), stacked = True)
        plt.vlines([minZ, maxZ], 0, np.max(h[:1]), label = "min and max Z positions")
    else:
        h1 = plt.hist(targetZs, bins = bins, density = True, alpha = 0.8, 
                 label = "{} distribution".format(targetElem))
        h2 = plt.hist(latticeZs, bins = bins, density = True, alpha = 0.2, 
                 label = "{} distribution".format(latticeElems))
        plt.vlines([minZ, maxZ], 0, np.max([h1[:1], h2[:1]]), label = "min and max Z positions")


    plt.legend()
    plt.show()