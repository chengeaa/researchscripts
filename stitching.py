from ase.io import xyz, gen, extxyz
from researchscripts.analysis import getFragsTraj
import time

##################################################################################
## in this file, I define various utilities to address the fact that deleting   ##
## atoms will lead to reindexing by ASE (no gaps in index sequences are allowed ##
##################################################################################

def reindexTrajBreak(tempIn, tempOut, oldIndices):
    """
    Get aligned indices for tempOut given also tempIn and oldIndices (adjusted indices for tempIn), 
    adjusting for loss of deleted or added atoms
    """
    tempIn = tempIn.copy()
    tempOut = tempOut.copy()
    newIndices = []
    
    for index, atom in zip(oldIndices, tempIn):
        # create a copy of the output geometry with the addition of the input atom of interest
        tempRef = tempOut.copy()
        tempRef.append(atom)
        dists = tempRef.get_all_distances(mic = True)
        
        # now we check if this atom still exists in the output
        
        #not 0, numerical tolerance
        if np.any(dists[-1, :-1] < 1e-5): 
            newIndices += [index]
            
    
    padLen = len(tempOut) - len(newIndices)
    output = np.concatenate((newIndices, np.zeros(padLen, dtype = int)))
    return output
def stichTrajectories(bombID, bombardments = 1, prefix = "geom.out"):
    """
    bombID (string): a string in the format of "{batch}-{sample}", eg "5-3"
    bombardments (arraylike): result of range() or np.arange() specifiying which bombardment events to use
    otherwise, should be an int, and the range (0, bombardments) will be used 
    filename (string): name of .xyz trajectory file in each replicate
    """
    startTime = time.time()
    bombardments = bombardments if hasattr(bombardments, '__iter__') else np.arange(bombardments)
    batch, sample = bombID.split("-")
    
    # Create an empty array to populate with the trajectory
    
    print("bombardments:", bombardments)
    nMaxAtoms = np.max(
        [len(
            gen.read_gen("{}/{}/{}/{}/geom.out.gen".format(_b, step, batch, sample, prefix))
        )
         for _b in bombardments for step in ['bomb', 'quench', 'eq']
        ]
    )
    
    # nMaxAtoms is guaranteed to correspond to some frame with an Ar in it
    # then, I basically generate nBombardments slots at the bottom of the df for the Ar atoms introduced
    nMaxAtoms = nMaxAtoms - 1 + len(bombardments) 
    print("nMaxAtoms: ", nMaxAtoms)
    
    #create the final matrix that will actually represent stitched traj
    trajFrame = np.zeros((nMaxAtoms, 0), dtype = object)
    
    frameIdx = 0 #initialize global frame count

    for _b in bombardments:
        with open("{}/bomb/{}/{}/{}.xyz".format(_b, batch, sample, prefix)) as f1:
            with open("{}/quench/{}/{}/{}.xyz".format(_b, batch, sample, prefix)) as f2:
                with open("{}/eq/{}/{}/{}.xyz".format(_b, batch, sample, prefix)) as f3:
#                     print("{}/bomb/{}/{}/{}.xyz".format(_b, batch, sample, prefix))
                    _btemp = list(extxyz.read_extxyz(f1, index = slice(0, None)))
                    _qtemp = list(extxyz.read_extxyz(f2, index = slice(0, None)))
                    _etemp = list(extxyz.read_extxyz(f3, index = slice(0, None)))
                    trajList = _btemp + _qtemp + _etemp #list form of this thing
                    
                    if _b == bombardments[0]: 
                        lastLen = len(_btemp[0]) # initialize lastLen in the very first frame
#                     fragNames, fragIdxs = getFragsTraj(trajList)                   
                        newIndices = np.arange(lastLen)
                    for frame in trajList:
#                         try:
#                         print(frame)
                        frame.set_momenta(frame.get_masses().reshape(-1,1) * frame.arrays['vel'])
#                         except:
#                             pass
#                             print("no velocity data")
                        trajFrame = np.hstack((trajFrame, np.zeros((nMaxAtoms, 1), dtype = object)))
                        if len(frame) != lastLen: # check for changes in nAtoms present in frame; indicative of step change
#                             print("entered if block for ", frameIdx)
                            newIndices = reindexTrajBreak(prevFrame, frame, newIndices)
                            ArAdded = (
                                frame[-1].symbol == "Ar" and
                                np.all(frame[-1].position != prevFrame[-1].position)
                            )
                            if ArAdded:
                                newIndices[-1] = nMaxAtoms - (len(bombardments) - _b)
                                
                            lastLen = len(frame)
                        trajFrame[newIndices, frameIdx] = frame
                        prevFrame = frame
                        frameIdx += 1
#                     trajFrame[trajFrame == 0] = np.nan
    endTime = time.time()
    print("execution time (s) = {}".format((endTime - startTime)))
    return pd.DataFrame(trajFrame)

