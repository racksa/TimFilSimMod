#!/usr/bin/python3

import driver
import arnoldi
import numpy as np
import time
import sys
import os
import util

def main():
    # Initialisation
    NSEG = 20      # Number of segments
    NFIL = 159       # Number of filaments
    NBLOB = 9000
    AR = 8
    
    # Number of time steps (ndts) and fixT
    ndts = 300

    Num_evals = 6
    Ustar = np.zeros(2*NFIL)

    # Initialise the driver
    d = driver.DRIVER()
    d.cuda_device = 1
    d.date = 'soln'
    d.category = 'JFNK_sims/'
    d.dir = f"data/{d.category}{d.date}/"
    eval_filename =  d.dir + f"eigenvalues_pb.dat"
    evec_filename =  d.dir + f"eigenvectors_pb.dat"

    output_filename = d.dir + f"psi_guess{NFIL}.dat"

    with open(output_filename, 'r') as file:
        num_lines = sum(1 for line in file)
        solns = np.loadtxt(output_filename)

    num_lines = 1

    for i in range(num_lines):
        k = solns[i,0]
        T = solns[i,1]
        dt = T/ndts
        Ustar = solns[i,2:]
        Ustar[:NFIL] = util.box(Ustar[:NFIL], 2*np.pi)
        print(f'k={k} T={T}')

        a = arnoldi.ARNOLDI(NSEG, NFIL, NBLOB, AR, T,\
                            Ustar,Num_evals,k,d, eval_filename, evec_filename)

        # Begin arnoldi iterations
        a.arnoldi_for_eigenvalues(T)

main()