import numpy as np
import newton_solver
import os
import time

def main():
    # Define 'global variables'
    epsJ = 1e-5 # 1e-6  # epsilon used in Jacobian approximation
    #ndts = 200   # Number of timesteps taken in period T
    #fixT = 1   # Fix T for equilibrium, rather than PO solution
    #follower_force = 50      # Follower force (parameter of dynamical system)
    NSEG = 20      # Number of segments
    NFIL = 1       # Number of filaments
    y_spacing_const = 5
    x_spacing_const = y_spacing_const
    output_filename = "psi.dat"

    # Number of time steps (ndts) and fixT
    ndts = 300
    fixT = 0
 
    # n = 3*(NSEG-1)*NFIL+1 #4#3 * (N - 1) * Nf + 1  # Dimension of system, including unknown params
    n = 2*639+1
    mgmres = 5  # 10  # max GMRES iterations
    nits = 150  # max Newton iterations
    rel_err = 1e-6  # 1e-8 Relative error |F|/|x|
    del_value = -1  # These rarely need changing for any problem
    mndl = 1e-20
    mxdl = 1e20
    gtol = 1e-3  # 1e-4

    f_range = np.array([0.06])
    for k in f_range:

        print('-----------Spring constant = ' + str(k))
    
        new_x = find_new_x(fixT,NSEG,NFIL,output_filename)
        newton = newton_solver.NEWTON_SOLVER(new_x,epsJ,ndts,fixT,k,NSEG,NFIL, y_spacing_const, x_spacing_const)

       # Scale parameters by |x| then call black box
        d = np.sqrt(np.sum(newton.new_x[1:] * newton.new_x[1:]))
        # Do we want rel_err?
        tol = rel_err * d
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

        # print(os.getcwd())
        # print(input_filename)
        new_x = np.loadtxt(input_filename) #TODO: fix this!
        
        new_x[-1,2:2+639] = new_x[-1,2:2+639] - np.floor(new_x[-1,2:2+639]/(2*np.pi))*2*np.pi
        # print(np.shape(new_x))
        # print(new_x[-1,1:])
        return new_x[-1,1:]

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
