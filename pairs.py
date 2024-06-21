# finds all possible pairs for a set of vectors*, with no duplicate pairs (eg, [B, A] when [A, B] already exists) or pairs of the same vector (eg, [A, A])
# *(unit-size vectors, with x coord, y coord and angle)

import numpy as np
import math

# create a random array of vectors to demonstrate the function:

def fill(size):
    array = np.zeros((size, 3))
    i = 0
    while i < array.shape[0]:
        array[i,0] = np.random.randint(1, 10)
        array[i,1] = np.random.randint(1, 10)
        array[i,2] = np.random.randint(0, 360)
        i = i + 1
    return array
        
a = fill(167)

# pair finding function:

def pair(array):
    length = array.shape[0]
    width = array.shape[1]
    pairs = np.zeros(((math.comb(length, 2)), 2*width))
    counter = 1
    j = 0
    k = 1
    jmax = length - counter

    while counter < (length): 
        while j < jmax:
            pairs[j,0] = array[(counter-1), 0]
            pairs[j,1] = array[(counter-1), 1]
            pairs[j,2] = array[(counter-1), 2]
            pairs[j,3] = array[k, 0]
            pairs[j,4] = array[k, 1]
            pairs[j,5] = array[k, 2]
            j = j + 1
            k = k + 1
        counter = counter + 1
        k = counter
        jmax = jmax+(length - counter)
    return pairs

pairs = pair(a)
