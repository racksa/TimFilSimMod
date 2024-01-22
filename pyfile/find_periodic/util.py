import numpy as np

def box(x, box_size):
    return x - np.floor(x/box_size)*box_size