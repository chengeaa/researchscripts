#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import os
import matplotlib.animation as animation
import sys

def eigenval_analyze(filename='OUTCAR', unit = 'THz'):
    """
    Searches an OUTCAR for the specified unit, and returns a
    numpy array of the normal modes of that unit 
    """
    eigenval_pattern = r"-?(\d+\.\d+) " + unit
    collect = []
    with open(filename) as c:
        for line in c:
            m = re.search(eigenval_pattern, line) 
            if m is not None:
                result = float(m.groups()[0])
                collect += [(result)]
    return np.array(collect)

if __name__ == "__main__":
    """
    Arg 1: OUTCAR
    Arg 2: term
    """
    args = sys.argv[1:]
    if len(args) > 2:
        print(args)
        raise Exception("Too many args")
    print(eigenval_analyze(*args))
