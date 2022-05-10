from ase import io 
from ase.visualize import view
from ase.cell import Cell
from ase.geometry.analysis import Analysis
import matplotlib.pyplot as plt
from researchscripts.analysis import coordLabeller   
import pandas as pd

def main(filepath, a1, a2, x, y, z, outputname)

    #parameters for orthorhombic cell with periodic boundaries 

    s = io.read(filepath)
    c = Cell([
        [x, 0, 0],
        [0, y, 0],
        [0, 0, z],
    ])
    s.cell = c
    s.pbc = True


    ##########################################
    ## TUNE THESE VALUES FOR YOUR SYSTEM!!! ##
    ##########################################

    rel, b = coordLabeller(s, 
                  fullCoordinations={"Li": 1, "F": 1, "O": 2},
                  minAngles = {"Li": 90, "F": 120, "O": 105},
                  maxBonds_per_element={"Li": 6, "F": 3, "O": 4},
                  angle_tolerance=0,
                  bond_tolerance=0,
                  minz = -1
                 )

    _values = {f"F{i}-Li{k}": a.get_bond_value(0,(i, k))
         for i, j in b.items() for k in j
         if s[i].symbol == a1 and s[k].symbol == a2
        }


    plt.hist(_values.values())
    plt.savefig(outputname)

    pd.Series(_values).to_csv(outputname)

    print("number of bonds detected: ", len(_values.values()))



if __name__ == "__main__":
    filepath = "../O64_2X2X2_last.xyz"
    a1, a2 = "Li", "F"
    x, y, z = 28.735385, 26.490013, 34.503785
    outputname = "bonds"
        """
    Arg 1: filepath
    Arg 2: a1, element 1 of interest
    Arg 3: a2, element 2 of interest
    Arg 3, 4, 5: x, y, z dimensions of orthorhombic cell
    Arg 6: name for output files
    """
    args = sys.argv[1:]
    if len(args) > 6:
        print(args)
        raise Exception("Too many args")
    if len(args == 0):
        print("using default parameters")
        args = [filepath, a1, a2, x, y, z]
    main(*args)
