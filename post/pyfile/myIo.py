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

def get_boxsize_from_name(filename):
    str_list = filename.split('_')
    try:
        Lx, Ly, Lz = [float(s) for s in str_list[-3:]]
        return Lx, Ly, Lz
    except:
        print("WARNING: Filename not supported for auto boxing.")
        return (float('inf'), float('inf'), float('inf'))
    
def get_ciliate_data_from_name(filename):
    str_list = filename.split('_')
    try:
        R, Tor = 0, 0
        for s in str_list[-2:]:
            print(s)
            if(s.endswith('R')):
                R = float(s[:-1])
            if(s.endswith('torsion')):
                Tor = float(s[:-7])
        return R, Tor
    except:
        print("WARNING: Filename not supported for auto ciliating.")
        return (float('inf'), float('inf'), float('inf'))