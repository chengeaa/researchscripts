import pandas as pd
import numpy as np

def postprocessResults(directory = "../"):
    """
        Takes in a list of indices, corresponding to the bombardment trials to analyze
        Looks for files named ``results$i{bomb,quench,eq}.csv`` in directory specified. 
        Returns list of 3 dfs; each one has elements and keys
    """

    subdirs = np.arange(10)
    bombdata = {i :pd.read_csv(directory + "results%dbomb.csv" % i, index_col=0) 
            for i in subdirs}
    quenchdata = {i : pd.read_csv(directory + "results%dquench.csv" % i, index_col=0) 
            for i in subdirs}
    eqdata = { i: pd.read_csv(directory + "results%deq.csv" % i, index_col=0) 
            for i in subdirs}

    return [bombdata, quenchdata, eqdata]



def postprocessAggregated(simindices, directory = "../"):
    """
        Takes in a list of indices, corresponding to the bombardment trials to analyze
        Looks for files named ``aggregated_{bomb,quench,eq}$i`` in directory specified. 
    """
    bombdata = {i :pd.read_csv(directory + "aggregated_bomb%d.csv" % i, index_col=0) 
            for i in simindices}
    quenchdata = {i : pd.read_csv(directory + "aggregated_quench%d.csv" % i, index_col=0) 
            for i in simindices}
    eqdata = { i: pd.read_csv(directory + "aggregated_eq%d.csv" % i, index_col=0) 
            for i in simindices}
    data = {"bomb": bombdata, "quench": quenchdata, "eq":eqdata}


    aggregated = {}
    for i in simindices:
        for step in ["bomb", "quench", "eq"]:
            aggregated["%i-%s" % (i, step)] = data[step][i].sum(axis = 1)
    return pd.DataFrame(aggregated)
