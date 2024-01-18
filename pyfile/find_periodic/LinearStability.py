#!/usr/bin/python3

import driver
import arnoldi
import numpy as np
import time
import sys
import os

def main():
    # Initialisation
    nseg = int(20)
    nfil = int(1)
    f = float(5.0)
    k = float(1000.0)
    y_spacing_const = float(5) # Note: this is applied as a spacing of L/y_spacing_const
    x_spacing_const = float(5) # Note: this is applied as a spacing of L/x_spacing_const
    ndts = int(10)
    dt = float(0.01)

    sim_3D = 1 #0 if 2D, otherwise 3D
    T = dt*ndts
    Num_evals = 6

    #y_spacing_const_to_cycle = [0.5,1,2,4,8,10,20,30,40]


    Ustar = np.zeros(3*(nseg-1)*nfil)

    print(os.getcwd())
    new_x = np.loadtxt("PB_sol.dat")
    ndts = 200

    for i in range(37):
        f = new_x[i,0]
        dt = new_x[i,1]/ndts
        Ustar = new_x[i,2:]
        print("f = " + str(f) + ", T = " + str(new_x[i,1]) + ", Ustar = " +str(Ustar))
        #eval_filename = f"eigenvalues_{nfil}fil_{nseg}seg_{k}spring_{y_spacing_const}spacing_varying_f.dat"
        #evec_filename = f"eigenvectors_{nfil}fil_{nseg}seg_{k}spring_{y_spacing_const}spacing_varying_f.dat"

        eval_filename = f"eigenvalues_pb.dat"
        evec_filename = f"eigenvectors_pb.dat"

        d = driver.DRIVER()
        a = arnoldi.ARNOLDI(nseg,nfil,Ustar,Num_evals)

        # These are always changed at the beginning of each Arnoldi iteration 
        d.change_variables(nseg, nfil, f, k, y_spacing_const, x_spacing_const, ndts, dt)
        simulation_file = d.update_globals_file()
        #d.run()

        #TODO: delete pre-existing data files?

        # Begin arnoldi iterations
        a.arnoldi_for_eigenvalues(sim_3D, Num_evals, simulation_file, d.dir, T, f, eval_filename, evec_filename)

main()