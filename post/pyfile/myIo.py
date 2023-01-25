import os
import pandas as pd

def read_blob_references(fileName):
    ret = pd.read_csv(fileName, sep = ' ', header=None)
    return ret.iloc[:,:-1].to_numpy().reshape(-1)

def read_fil_references(fileName):
    ret = pd.read_csv(fileName, sep = ' ', header=None)
    return ret.iloc[:,:-1].to_numpy().reshape(-1)

def read_body_states(fileName):
    ret = pd.read_csv(fileName, sep=' ', header=None)
    ret = ret.iloc[:, 1:-1]
    return ret.to_numpy()

def read_seg_states(fileName):
    ret = pd.read_csv(fileName, sep = ' ', header=None)
    ret = ret.iloc[:, 1:-1]
    return ret.to_numpy()
    
def read_pars(fileName):
    ret_pardict = {}
    df = pd.read_csv(fileName, sep=' %% ', header=None, engine='python')
    for i in range(len(df)):
        ret_pardict[df.iloc[i, 1]] = df.iloc[i, 0]
    return ret_pardict

def write_line(text, fileName):
    with open(fileName, 'a') as the_file:
        the_file.write(text + '\n')
    
def clean_file(fileName):
    open(fileName, 'w')