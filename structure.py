from ase.io import gen
import pandas as pd
import numpy as np
from ase.geometry.analysis import Analysis

def getslab(struct):
    """
    Input: 
        struct: structre from which we will trim unbound species (.gen file)
    Output:
        baseslab: structure with unbound species trimmed
    """
    adjmat = Analysis(struct).adjacency_matrix[0]
    numnodes = adjmat.shape[0]
    g = Graph(numnodes)
    for i in range(numnodes):
        for j in range(numnodes):
            if adjmat[i,j]:
                g.addEdge(i,j)
    cc = g.connectedComponents()
    maingraph = np.array([i for i in cc if 0 in i][0])
    return struct[[atom.index for atom in struct if atom.index in maingraph]]

def getslab(struct):
    """
    Input: 
        struct: structre from which we will trim unbound species (.gen file)
    Output:
        baseslab: structure with unbound species trimmed
    """
    adjmat = Analysis(struct).adjacency_matrix[0]
    numnodes = adjmat.shape[0]
    g = Graph(numnodes)
    for i in range(numnodes):
        for j in range(numnodes):
            if adjmat[i,j]:
                g.addEdge(i,j)
    cc = g.connectedComponents()
    maingraph = np.array([i for i in cc if 0 in i][0])
    return struct[[atom.index for atom in struct if atom.index in maingraph]]
# Python program to print connected
# components in an undirected graph
# https://www.geeksforgeeks.org/connected-components-in-an-undirected-graph/

def getslabs(data, directory, useInputs = False):
    """
    Utility for getting and writing slab files from readData (utils.py) function
    data is the df from readData function or any df with (struct, in) and (struct, out) columns
    """

    if useInputs:
        slabSource = data['struct']['in']
    else:
        slabSource = data['struct']['out']

    dataDir = directory
    slabs = {}

    # # to generate slabs

    for key, value in slabSource.iteritems():
        slabs[key] = getslab(value)

    for key, value in slabs.items():
        gen.write_gen(dataDir + "slab{}.gen".format(key), value)

    # to read slabs
    for key in data.index:
        slabs[key] = gen.read_gen(dataDir + "slab{}.gen".format(key))

    if useInputs:
        data.loc[:, ('struct', 'inslab')] = pd.Series(slabs)
    else:
        data.loc[:, ('struct', 'outslab')] = pd.Series(slabs)

class Graph:
 
    # init function to declare class variables
    def __init__(self, V):
        self.V = V
        self.adj = [[] for i in range(V)]
 
    def DFSUtil(self, temp, v, visited):
 
        # Mark the current vertex as visited
        visited[v] = True
 
        # Store the vertex to list
        temp.append(v)
 
        # Repeat for all vertices adjacent
        # to this vertex v
        for i in self.adj[v]:
            if visited[i] == False:
 
                # Update the list
                temp = self.DFSUtil(temp, i, visited)
        return temp
 
    # method to add an undirected edge
    def addEdge(self, v, w):
        self.adj[v].append(w)
        self.adj[w].append(v)
 
    # Method to retrieve connected components
    # in an undirected graph
    def connectedComponents(self):
        visited = []
        cc = []
        for i in range(self.V):
            visited.append(False)
        for v in range(self.V):
            if visited[v] == False:
                temp = []
                cc.append(self.DFSUtil(temp, v, visited))
        return cc

# This code is contributed by Abhishek Valsan

