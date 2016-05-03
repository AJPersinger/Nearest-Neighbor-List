#  Title: Nearest Neigbhor List
#  Purpose: Implement nearest neigbhor searches and creates a list of nearest atoms
#             in a structure
#  Author: Axel-Jose Persinger

from scipy.spatial import cKDTree
import scipy as sp
import numpy as np
import time

# --------------------------------------------------------------------------------------------- #
#  Gobal Variables:
#    numberDim is the number of dimensions
#    t is the initial time for tic-tock functionality
# --------------------------------------------------------------------------------------------- #
numberDim = 3

# Tick!
t = time.time()



# --------------------------------------------------------------------------------------------- #
#  Takes a file as an input and creates a list of all coordinates returns a list of atom
# --------------------------------------------------------------------------------------------- #
def atomLister(file):
    # Opens a LAMMPS dump file named  'dump' and reads into a variable
    with open(file) as file:
        fileList = file.read().splitlines()

    # Removes uncessesary header information
    fileList = fileList[9:]

    # Initializes the list to store atoms
    atomList = [line.split()[2:-6] for line in fileList]

    return atomList
# --------------------------------------------------------------------------------------------- #



# --------------------------------------------------------------------------------------------- #
#  Takes the atom list and creates a KDTree. Returns the nearest neighbors for all atoms
# --------------------------------------------------------------------------------------------- #
def kdTree(atomList, distance):
    # Converts the atomList into an array
    atomArray = np.array(atomList).astype(np.float)

    # Changes the shape of the array to handle the dimensions
    atomArray.shape = atomArray.size / numberDim, numberDim

    # Creates ckdTree
    tree = cKDTree(atomArray, leafsize=atomArray.shape[0]+1)

    # Uses list comprehension to store all nearest neighbors
    ndxList = [(str(tree.query([atom], k = distance))) for atom in atomList]
    # print ndxList

    # Creates final comprehension to store all nearest neighbors in the form:
    #  neighborList:[atom][neighbors]
    # neighborList = [[atom for atom in atomArray] for string in ndxList]

    #neighborList = [[ndxList[i]] for i in range(10)]
    #del neighborList[::2]
    #neighborList.pop(0)
    #print ndxList

    return ndxList

# --------------------------------------------------------------------------------------------- #



# Replace 'out.Pt_3nm.20000' with a LAMMPS dump file
listOfAtoms = atomLister('out.Pt.0')
# Replace 10 with a distance
print kdTree(listOfAtoms, 3.332)

# Tock!
print time.time() - t
