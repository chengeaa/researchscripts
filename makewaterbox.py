#!/usr/bin/env python
'''Generates water box

Writes text files for O and H atoms that can be pasted into POSCAR to fill in 
water box 

Args:
    n - number of layers
    num_per_layer - number of waters per layer
    correction - 0 if a layer is missing in default output

Returns:
    None; writes files o.txt and h.txt 

'''

import numpy as np

def main(n, num_per_layer, correction = 0):

    #reference water location
    initial_O = np.array([0.01, 0.01, 0.01]) 
    initial_H1 = np.array([0.01, 0.01, 0.18]) 
    initial_H2 = np.array([0.01, 0.14, 0.01]) 

    #initialize ouptut lists
    O_output = np.array([])
    H_output = np.array([])

    #reference to a vector
    lowerlim = 0.16
    upperlim = 0.60 
    vacuum_height = 1 - (upperlim - lowerlim)

    #set spacings
    a_space = vacuum_height/(n / num_per_layer) #assume even n and 2 waters per layer
    b_space = 1/2 
    c_space = 1/4 

    for i in np.arange(int(1/a_space) + correction):
        for j in np.arange(num_per_layer):
            #offsets a and c locations in step, b location within the layer
            if i * a_space < lowerlim or i * a_space > upperlim:
                O_output = np.append(O_output, np.array(initial_O + np.array([i * a_space, j * b_space + i*c_space, i*c_space])))
                O_output = np.append(O_output, np.array(["  T"]*3 ))
                H_output = np.append(H_output, np.array(initial_H1 + np.array([i * a_space, j * b_space+ i*c_space, i*c_space])))
                H_output = np.append(H_output, np.array(["  T"]*3 ))
                H_output = np.append(H_output, np.array(initial_H2 + np.array([i * a_space, j * b_space+ i*c_space, i*c_space])))
                H_output = np.append(H_output, np.array(["  T"]*3 ))
    O_output = O_output.reshape(-1, 6)
    H_output = H_output.reshape(-1, 6)

    np.savetxt("../o_test.txt", O_output, fmt = "%s")
    np.savetxt("../h_test.txt", H_output, fmt = "%s")

if __name__ == "__main__":
    """
    Arg 1: number of waters
    Arg 2: number of waters per layer 
    Arg 3: correction (1 if default run has one layer too few)
    """
    args = sys.argv[1:]
    if len(args) > 3:
        print(args)
        raise Exception("Too many args")
    main(*args)


