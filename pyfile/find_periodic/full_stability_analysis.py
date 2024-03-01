import numpy as np
import driver
import arnoldi
import newton_solver
import os

def main():
    
    k = 100            # Spring constant
    NSEG = 20      # Number of segments
    NFIL = 2       # Number of filaments
    y_spacing_const = 0.5 # Note: this is applied as a spacing of L/y_spacing_const
    x_spacing_const = y_spacing_const # Note: this is applied as a spacing of L/x_spacing_const
    #follower_force = 50
    fixT = 0 # Steady solution (fixT = 1) or periodic (fixT = 0)
#    input_filename = "FF_Length_20NSEG_Whirling_solutions.dat"
#    JFNK_output_filename =  "FF_Length_20NSEG_Whirling_solutions.dat"
#    eigenvalue_output_filename = f"FFLength_{NSEG}NSEG_Whirling_eigenvalues.dat"
#    eigenvector_output_filename = f"FFLength_{NSEG}NSEG_Whirling_eigenvectors.dat"

    input_filename = "psi.dat"
    JFNK_output_filename =  "test_JFNK.dat"
    eigenvalue_output_filename = "TwoFil_axonemal_20NSEG_eigenvalues.dat"
    eigenvector_output_filename = "TwoFil_axonemal_20NSEG_eigenvectors.dat"
    
    #eigenvalue_output_filename = f"FFLength_{NSEG}NSEG_2D_eigenvalues.dat"
    #eigenvector_output_filename = f"FFLength_{NSEG}NSEG_2D_eigenvectors.dat"
    new_x = find_new_x(fixT,NSEG,NFIL,input_filename)
    print(new_x)
    ndts = find_ndts(fixT) # Number of timesteps

    #stability_analysis_singlepoint(k,NSEG,NFIL,y_spacing_const,x_spacing_const,follower_force,fixT,new_x,ndts,JFNK_output_filename,eigenvalue_output_filename,eigenvector_output_filename)
    stability_analysis_continuation(k,NSEG,NFIL,y_spacing_const,x_spacing_const,fixT,new_x,ndts,JFNK_output_filename,eigenvalue_output_filename,eigenvector_output_filename)
    #stability_analysis_verification(k,NSEG,NFIL,y_spacing_const,x_spacing_const,fixT,ndts,input_filename,JFNK_output_filename,eigenvalue_output_filename,eigenvector_output_filename)
    #stability_analysis_vertical_equilib(k,NSEG,NFIL,y_spacing_const,x_spacing_const,fixT,ndts,eigenvalue_output_filename,eigenvector_output_filename)
    return 0


def stability_analysis_singlepoint(k,NSEG,NFIL,y_spacing_const,x_spacing_const,follower_force,fixT,new_x,ndts,JFNK_output_filename,eigenvalue_output_filename,eigenvector_output_filename):
    # If trying one value (steady or periodic)
    new_x = track_solution(k,NSEG,NFIL,y_spacing_const,x_spacing_const,follower_force,ndts,fixT,new_x,JFNK_output_filename)
    find_eigenmodes(k,NSEG,NFIL,y_spacing_const,x_spacing_const,follower_force,ndts,fixT,new_x,eigenvalue_output_filename,eigenvector_output_filename)

def stability_analysis_continuation(k,NSEG,NFIL,y_spacing_const,x_spacing_const,fixT,new_x,ndts,JFNK_output_filename,eigenvalue_output_filename,eigenvector_output_filename):
    
    # Continuation
    f_range = np.array([10,11,12,13,14,15,16,17,18,19,20])
    f_range = np.array([0.03])
    #f_range = np.linspace(20,150,131)
    #f_range = np.array([103,104])
    print(ndts)
    print(new_x[0]/ndts)

    for follower_force in f_range:
        print("Follower force = " + str(follower_force) + ", new_x = " + str(new_x))
        new_x = track_solution(k,NSEG,NFIL,y_spacing_const,x_spacing_const,follower_force,ndts,fixT,new_x,JFNK_output_filename)
        print('JFNK analysis complete')
        # find_eigenmodes(k,NSEG,NFIL,y_spacing_const,x_spacing_const,follower_force,ndts,fixT,new_x,eigenvalue_output_filename,eigenvector_output_filename)
        # print('Stability analysis complete')

def stability_analysis_verification(k,NSEG,NFIL,y_spacing_const,x_spacing_const,fixT,ndts,input_filename,JFNK_output_filename,eigenvalue_output_filename,eigenvector_output_filename):
    dat = np.loadtxt(input_filename); 
    sz = dat.shape
    print(sz)
    for i in range(sz[0]):
        follower_force = dat[i,0]
        new_x = dat[i,1:]
        print('follower force = ' + str(follower_force))
        new_x = track_solution(k,NSEG,NFIL,y_spacing_const,x_spacing_const,follower_force,ndts,fixT,new_x,JFNK_output_filename)
        print('JFNK analysis complete')
        find_eigenmodes(k,NSEG,NFIL,y_spacing_const,x_spacing_const,follower_force,ndts,fixT,new_x,eigenvalue_output_filename,eigenvector_output_filename)
        print('Stability analysis complete')

def stability_analysis_vertical_equilib(k,NSEG,NFIL,y_spacing_const,x_spacing_const,fixT,ndts,eigenvalue_output_filename,eigenvector_output_filename):

    new_x = np.concatenate(([0.4],np.zeros(3*(NSEG-1)*NFIL)))
    f = np.linspace(0,10)
    for i in range(f.shape[0]):
        follower_force = f[i]
        find_eigenmodes(k,NSEG,NFIL,y_spacing_const,x_spacing_const,follower_force,ndts,fixT,new_x,eigenvalue_output_filename,eigenvector_output_filename)

def find_new_x(fixT,NSEG,NFIL,filename):
    if fixT == 1: return np.concatenate(([0.2],1e-2 * np.random.standard_normal(3*(NSEG-1)*NFIL))) # Current best x and ndts
    else: new_x = np.loadtxt(filename); return new_x[-1,1:]

def find_ndts(fixT):
    if fixT == 1: return 10
    else: return 100#200

def track_solution(k,NSEG,NFIL,y_spacing_const,x_spacing_const,follower_force,ndts,fixT,new_x,output_filename):

    # Define 'global variables'
    epsJ = 1e-5 # 1e-6  # epsilon used in Jacobian approximation
    n = 3*(NSEG-1)*NFIL+1  # Dimension of system, including unknown params
    mgmres = 5  # 10  # max GMRES iterations
    nits = 150  # max Newton iterations
    rel_err = 1e-6  # 1e-8 Relative error |F|/|x|

    # These rarely need changing for any problem:
    del_value = -1  
    mndl = 1e-20
    mxdl = 1e20
    gtol = 1e-3  # 1e-4

    newton = newton_solver.NEWTON_SOLVER(new_x,epsJ,ndts,fixT,k,follower_force,NSEG,NFIL, y_spacing_const, x_spacing_const)

    # Scale parameters by |x| then call black box
    d = np.sqrt(np.sum(newton.new_x[1:] * newton.new_x[1:]))

    if fixT == 1:
        tol = rel_err
        del_value = del_value
        mndl = mndl 
        mxdl = mxdl 
    else:
        tol = rel_err * d
        del_value = del_value * d
        mndl = mndl * d
        mxdl = mxdl * d

    info = 1

    info = newton.NewtonHook(mgmres, n, gtol, tol, del_value, mndl, mxdl, nits, info)

    print(newton.new_x)

    save_solution(np.concatenate(([follower_force],newton.new_x)),output_filename)

    return newton.new_x


def save_solution(data,filename):
    os.chdir("../stability/")

    with open(filename, "ab") as  f:
        f.write(b"\n")
        np.savetxt(f, data, newline = " ")

    return


def find_eigenmodes(k,nseg,nfil,y_spacing_const,x_spacing_const,f,ndts,fixT,new_x,eval_filename,evec_filename):
    # Initialisation
    sim_3D = 0 #0 if 2D, otherwise 3D
    Num_evals = 6

    if (fixT == 1): # Steady state (vertical equilibrium
        Ustar = np.zeros(3*(nseg-1)*nfil)
    else:           # Periodic solution, given by fed-in data
        Ustar = new_x[1:]
    
    T = new_x[0]
    dt = T/ndts
    
    d = driver.DRIVER()
    a = arnoldi.ARNOLDI(nseg,nfil,Ustar,Num_evals)

    # These are always changed at the beginning of each Arnoldi iteration 
    d.change_variables(nseg, nfil, f, k, y_spacing_const, x_spacing_const, ndts, dt)
    simulation_file = d.update_globals_file()

    # Begin arnoldi iterations
    a.arnoldi_for_eigenvalues(sim_3D, Num_evals, simulation_file, d.dir, T, f, eval_filename, evec_filename)

main()