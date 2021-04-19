"""
Input: a list of directory of training datasets
Each directory containins:
    -.gen structure files (each with a single adsorbate) 
    - an ``energies`` file containing the calculated adsorption energies
Returns: 
    trained KRR models for z prediction and E prediction

Also writes .pkl of each model to current directory
"""
from dscribe.descriptors import SOAP


soap = SOAP(
    species=species,
    periodic=True,
    rcut=5,
    nmax=10,
    lmax=9,
)

print("number of soap features: {}".format(
    soap.get_number_of_features())
    )


# globals
jitter = True
jitterscale = 0.25
adslen = 5
# adsIdx = 370

species = ["Si", "N", "H", "C", "F", "Ar", "He"]

# specifics
datadirs = ["../mef.noH.bulk/"]

datalist = []
for datadir in datadirs:
    data = readStructs(datadir) #MeF, amorphous slab, DFTB
    data['processed'] = [
        convertAdsorbateToHe(
            i, len(i) - adslen, np.arange(len(i) - adslen, len(i))
        ) for i in data['geom']
    ]

    if jitter:
        np.random.seed(429)
        for struct in data['processed']:
            struct[-1].position[2] += np.random.normal(scale = jitterscale) 
            # add a bit of noise to z coord; reduce overfitting

    data = pd.concat([data, getSOAPs(data['processed'], species = species)], axis = 1)

    originalColumns = data.columns

    datalist += [pd.concat([pd.DataFrame(data['SOAP'].to_list(), index = data.index), data], axis = 1)]
data = pd.concat(datalist).reset_index(drop = True).fillna(0)
