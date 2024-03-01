#!/usr/bin/python3

import driver
import arnoldi
import numpy as np
import time
import sys
import os

def main():
    # Initialisation
    NSEG = 20      # Number of segments
    NFIL = 159       # Number of filaments
    NBLOB = 9000
    AR = 8
    
    folder = f"data/JFNK_sims/20240214_periodic_s"
    output_filename = folder + "/" + f"psi_guess{NFIL}.dat"

    # Number of time steps (ndts) and fixT
    ndts = 300

    nseg = 20
    nfil = 1
    f = 5.0
    k = 1000.0
    ndts = 300
    

    Num_evals = 6


    Ustar = np.zeros(2*nfil)

    with open(output_filename, 'r') as file:
        num_lines = sum(1 for line in file)
        solns = np.loadtxt(output_filename)

    num_lines = 10

    for i in range(num_lines):
        k = solns[i,0]
        T = solns[i,1]
        dt = T/ndts
        Ustar = solns[i,2:]
        print(f'k={k} T={T}')

        eval_filename =  folder + f"/eigenvalues_pb.dat"
        evec_filename =  folder + f"/eigenvectors_pb.dat"

        d = driver.DRIVER()
        a = arnoldi.ARNOLDI(nseg,nfil,Ustar,Num_evals,d)

        # These are always changed at the beginning of each Arnoldi iteration
        d.change_variables(NFIL, NSEG, NBLOB, AR, k, T, 1.)
        d.update_globals_file()


        # Begin arnoldi iterations
        # a.arnoldi_for_eigenvalues(Num_evals, simulation_file, d.dir, T, f, eval_filename, evec_filename)

main()