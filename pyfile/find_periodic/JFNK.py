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
    k = 1000                 # Spring constant
    NSEG = 20      # Number of segments
    NFIL = 1       # Number of filaments
    y_spacing_const = 5
    x_spacing_const = y_spacing_const
    #input_filename = "beating_f_150_liealgebra_configuration.dat"
    output_filename = "psi15.dat"

    # Number of time steps (ndts) and fixT
    ndts = 200
    fixT = 0
 
    n = 3*(NSEG-1)*NFIL+1 #4#3 * (N - 1) * Nf + 1  # Dimension of system, including unknown params
    mgmres = 5  # 10  # max GMRES iterations
    nits = 150  # max Newton iterations
    rel_err = 1e-6  # 1e-8 Relative error |F|/|x|
    del_value = -1  # These rarely need changing for any problem
    mndl = 1e-20
    mxdl = 1e20
    gtol = 1e-3  # 1e-4

    #f_range = np.array([150,140,130,120,110,100])
    f_range = np.array([90,80,70,60,55,50,45,40])
    #f_range = np.array([150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,400])
    f_range = np.array([200])
    for follower_force in f_range:

        print('Follower force = ' + str(follower_force))
    
        new_x = find_new_x(fixT,NSEG,NFIL,output_filename)
    #     newton = newton_solver.NEWTON_SOLVER(new_x,epsJ,ndts,fixT,k,follower_force,NSEG,NFIL, y_spacing_const, x_spacing_const)

    #    # Scale parameters by |x| then call black box
    #     d = np.sqrt(np.sum(newton.new_x[1:] * newton.new_x[1:]))
    #     # Do we want rel_err?
    #     tol = rel_err * d
    #     del_value = del_value * d
    #     mndl = mndl * d
    #     mxdl = mxdl * d

    #     info = 1

    #     info = newton.NewtonHook(mgmres, n, gtol, tol, del_value, mndl, mxdl, nits, info)

    #     print(newton.new_x)

    #     save_solution(np.concatenate(([follower_force],newton.new_x)),output_filename)

def find_new_x(fixT,NSEG,NFIL,input_filename):

    if (fixT == 1):
        
        return np.concatenate(([0.2],1e-8 * np.random.standard_normal(3*(NSEG-1)*NFIL))) # Current best x
    
    else:

        print(os.getcwd())
        print(input_filename)
        new_x = np.loadtxt(input_filename) #TODO: fix this!
        print(np.shape(new_x))
        print(new_x[-1,1:])
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
