import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import os
import matplotlib.animation as animation
import sys

def oszicarSearch(filename="OSZICAR", term = "E0"):
    """
    Searches an OSZICAR for the specified statistic, and returns a numpy 
    array of values for that statistic from each ionic loop.
    """
    numberSearch = r"\s*-?\d*.\d+E?\+?\d+"
    searchTerm = term + "= "  + numberSearch
    collect = []
    with open(filename) as c:
        for line in c:
            m = re.search(searchTerm, line) 
            if m is not None:
                result = float(re.search(numberSearch, m[0])[0])
                collect += [result]
    return np.array(collect)

def main(filename="OSZICAR", outputname="energy", interval=16+2/3):
    """
    Calls OSZICARsearch to generate the potential and kinetic energy arrays;
    sums the two, and creates an animation 
    """
    interval = float(interval)
    E0, EK = oszicarSearch(filename, "E0"), oszicarSearch(filename, "EK")
    l = E0 + EK #sum of these numpy arrays
    lmax, lmin = max(l), min(l)
    lrange = lmax - lmin

    def data_gen(t=0):
        while t < len(l):
            t += 1
            yield t, l[t]


    def init():
        ax.set_ylim(lmin - .1 * lrange, lmax + .1 * lrange)
        ax.set_xlim(0, len(l))
        del xdata[:]
        del ydata[:]
        line.set_data(xdata, ydata)
        return line,

    fig, ax = plt.subplots()
    line,  = ax.plot([], [])
    ax.grid()
    xdata, ydata = [], []


    def run(data):
        # update the data
        t, y = data
        xdata.append(t)
        ydata.append(y)
        xmin, xmax = ax.get_xlim()
        line.set_data(xdata, ydata)

        return line,

    ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=interval,
                                  repeat=False, init_func=init)
    plt.xlabel("Timestep (fs)")
    plt.ylabel("Energy (eV)")
    ani.save_count = len(l) - 1
    ani.save(outputname+".mp4", writer = "ffmpeg")


if __name__ == "__main__":
    """
    Takes in two arguments:
    1. OSZICAR name
    2. Desired movie name
    3. Desired interval between frames (default 16+2/3)
    """
    args = sys.argv[1:]
    if len(args) > 3:
        print(args)
        raise Exception("No more than 3 arguments allowed")
    main(*args)