"""
Input: a list of directory of training datasets
Each directory containins:
    -.gen structure files (each with a single adsorbate) 
    - an ``energies`` file containing the calculated adsorption energies
Returns: 
    trained KRR models for z prediction and E prediction

Also writes .pkl of each model to current directory
"""
#############
## imports ##
#############

# custom
from dscribe.descriptors import SOAP
from krr_utils import convertAdsorbateToHe, getSOAPs
from utils import readStructs

#scipy
import pandas as pd
import numpy as np

rseed = 429
np.random.seed(rseed)

#sklearn 
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.kernel_ridge import KernelRidge

############
## params ##
############
jitter = True
subset = 10 #nrows//subset = sample size used
jitterscale = 0.25
adslen = 5 #length of adsorbate molecule in dataset
species = ["Si", "N", "H", "C", "F", "Ar", "He"] #include ALL elements possible in dataset
datadirs = ["../mef.noH.bulk/"]

class zmodel(KernelRidge

datalist = []
for datadir in datadirs:
    data = readStructs(datadir) #MeF, amorphous slab, DFTB
    data = data.loc[np.random.choice(data.index, size = data.shape[0]//subset, replace = False)]
    # below, we assume the adsorbate is the last ``adslen`` atoms in the file
    data['processed'] = [
        convertAdsorbateToHe(
            i, len(i) - adslen, np.arange(len(i) - adslen, len(i))
        ) for i in data['geom']
    ]

    if jitter:
        for struct in data['processed']:
            struct[-1].position[2] += np.random.normal(scale = jitterscale) 
            # add a bit of noise to z coord; reduce overfitting

    data = pd.concat([data, getSOAPs(data['processed'], rcut = 5, sigma = 0.1, species = species)], axis = 1)

    originalColumns = data.columns

    datalist += [pd.concat([pd.DataFrame(data['SOAP'].to_list(), index = data.index), data], axis = 1)]
data = pd.concat(datalist).reset_index(drop = True).fillna(0)
print(data.head())
nfeatures = data.shape - originalColumns.shape[0]
print(nfeatures)


# initial training
X_train, X_test, y_train, y_test = train_test_split(
    data["processed"],
    data["E_ads"], random_state = rseed)
print("# points total: %d; #train points: %d; #test points: %d" % 
      (sum(validData), len(X_train), len(X_test)))


alphas = np.logspace(-10, -1, 25)
gammas = np.logspace(-10, -1, 25)

krr = KernelRidge()  # need to set gamma by grid search.. 
Emodel = GridSearchCV(krr, [{"alpha":alphas, "gamma":gammas}], cv = 5)
Emodel.fit(np.array(X_train), np.array(y_train))


# post-test, train on whole dataset for maximum accuracy for production

# X_train = data.loc[validData, ~np.isin(data.columns, originalColumns)]
# y_train = data.loc[validData, 'E_ads']
