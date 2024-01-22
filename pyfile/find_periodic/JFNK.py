import numpy as np
import newton_solver
import os
import time
import util

def main():
    # Define 'global variables'
    epsJ = 1e-5 # 1e-6  # epsilon used in Jacobian approximation
    #ndts = 200   # Number of timesteps taken in period T
    #fixT = 1   # Fix T for equilibrium, rather than PO solution
    #follower_force = 50      # Follower force (parameter of dynamical system)
    NSEG = 20      # Number of segments
    NFIL = 159       # Number of filaments
    NBLOB = 5000
    AR = 6

    output_filename = "data/expr_sims/20240118_periodic/psi_guess.dat"

    # Number of time steps (ndts) and fixT
    ndts = 300
    fixT = 0
 
    # n = 3*(NSEG-1)*NFIL+1 #4#3 * (N - 1) * Nf + 1  # Dimension of system, including unknown params
    n = 2*NFIL+1
    mgmres = 2  # 10  # max GMRES iterations
    nits = 150  # max Newton iterations
    rel_err = 1e-3  # 1e-8 Relative error |F|/|x|
    del_value = -1  # These rarely need changing for any problem
    mndl = 1e-20
    mxdl = 1e20
    gtol = 1e-3  # 1e-4

    f_range = np.array([0.08])
    for k in f_range:

        print('-----------Spring constant = ' + str(k))

        new_x = find_new_x(fixT,NSEG,NFIL,output_filename)
        newton = newton_solver.NEWTON_SOLVER(new_x,epsJ,ndts,fixT,k,NFIL, NSEG, NBLOB, AR)

       # Scale parameters by |x| then call black box
        d = np.sqrt(np.sum(newton.new_x[1:] * newton.new_x[1:]))
        # Do we want rel_err?
        tol = rel_err * d
        print(f"Requested tol={tol}")
        del_value = del_value * d
        mndl = mndl * d
        mxdl = mxdl * d

        info = 1

        info = newton.NewtonHook(mgmres, n, gtol, tol, del_value, mndl, mxdl, nits, info)
        
    #     save_solution(np.concatenate(([follower_force],newton.new_x)),output_filename)

def find_new_x(fixT,NSEG,NFIL,input_filename):

    if (fixT == 1):
        
        return np.concatenate(([0.2],1e-8 * np.random.standard_normal(3*(NSEG-1)*NFIL))) # Current best x
    
    else:

        full_input = np.loadtxt(input_filename) #TODO: fix this!
        
        # full_input[2:2+NFIL] = np.sin(full_input[2:2+NFIL])
        # full_input[2:2+NFIL] = np.exp(1j*full_input[2:2+NFIL])
        full_input[2:2+NFIL] = util.box(full_input[2:2+NFIL], 2*np.pi)

        return full_input[1:]

def save_solution(data,filename):
    print(data)
    os.chdir("../stability/")
    input_filename = filename
    print(input_filename)

    with open(input_filename, "ab") as  f:
        f.write(b"\n")
        np.savetxt(f, data, newline = " ")

    return

# Run code
main()
