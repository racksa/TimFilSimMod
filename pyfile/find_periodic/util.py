import numpy as np
import math
import pandas as pd

def box(x, box_size):
    return x - np.floor(x/box_size)*box_size

def cartesian_to_spherical(x):
    """
    Convert Cartesian coordinates to spherical polar coordinates.
    
    Args:
        x (float, float, float): cartesian-coordinate.
    
    Returns:
        tuple: (r, theta, phi), where r is the radial distance, theta is the polar angle (azimuthal angle),
               and phi is the elevation angle (zenith angle).
    """
    r = math.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
    theta = math.atan2(x[1], x[0])
    phi = math.acos(x[2] / r)
    
    return r, theta, phi

def read_fil_references(fileName):
    ret = pd.read_csv(fileName, sep = ' ', header=None)
    return ret.iloc[:,:-1].to_numpy().reshape(-1)