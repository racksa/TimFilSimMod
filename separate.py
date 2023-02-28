from cProfile import label
from re import A
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import erf

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import subprocess
import os
import random
import pandas as pd

filePath = 'simulation_data.dat'
infoFile = open(filePath, 'r')
lines = infoFile.readlines()
N = 3920
infoFile.close()

read_frame = pd.read_csv(filePath, delimiter=' ', header=None, skiprows=int(0*(N+1)), nrows=N)

x = read_frame[0]
y = read_frame[1]
z = read_frame[2]

infoFile = open(filePath, 'r')
pos_lines = list()
for row in range(N):
    pos_string = str(x[row]) + ' ' + str(y[row]) + ' ' + str(z[row]) + ' 0.7184' + ' 0' + '\n'
    pos_lines.append(pos_string)

infoFile = open("pos_data.dat", 'w')
infoFile.writelines(pos_lines)
infoFile.close()