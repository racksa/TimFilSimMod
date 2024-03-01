import numpy as np
import gmresc
import driver
import os
import subprocess
import time
import util
import configparser
from scipy.sparse.linalg import gmres
from scipy.sparse.linalg import LinearOperator

class NEWTON_SOLVER:
    
    def __init__(self, new_x, epsJ, ndts, fixT, k, NFIL, NSEG, NBLOB, AR):
        self.new_x = new_x
        self.ndts = ndts
        self.k = k
        self.NSEG = NSEG
        self.NFIL = NFIL
        self.NBLOB = NBLOB
        self.AR = AR
        self.epsJ = epsJ
        self.fixT = fixT
        self.new_fx = []
        self.new_tol = []
        self.new_del = []
        self.new_nits = 0
        self.new_gits = 0
        self.dt = []

        self.d = driver.DRIVER()


    def dotprd(self, n_, a, b):
        # dot product.  n_==-1 is a flag to exclude parameter T.
        n1 = 0
        if n_ == -1:
            n1 = 1
        d = np.sum(a[n1:] * b[n1:])
        return d
    
    def Lorenz_f(self, x):
        dx = np.zeros_like(x)
        dx[0] = self.p[0] * (-x[0] + x[1])
        dx[1] = x[0] * (-x[2] + self.p[1]) - x[1]
        dx[2] = x[0] * x[1] - self.p[2] * x[2]
        return dx
    
    def run_filament_code(self, ndts, x):

        # Insert the spring constant to the left of the array
        x = np.insert(x, 0, self.k)

        # Save states to a file to read and run
        np.savetxt(self.d.dir + "psi.dat", x, newline = " ")

        # Change globals.ini file
        self.d.change_variables(self.NFIL, self.NSEG, self.NBLOB, self.AR, self.k, x[1], ndts/self.ndts)
        self.d.update_globals_file()

        # with open(self.d.dir + "pars.dat", "ab") as f:
        #     f.write(b"\n")
        #     np.savetxt(f, [x[1], ndts/self.ndts], newline = " ")

        # Run code
        # print(f"running new sim with T = {x[1]}/ steps = {ndts}")
        # config = configparser.ConfigParser()
        # config.read(self.d.globals_name)
        # print(config['Parameters']['sim_length'])
        self.d.run()

        # Return output
        output_filename = self.d.dir + self.d.simName + "_true_states.dat"
        return np.loadtxt(output_filename)

    def steporbit(self, ndts, x):
        
        if ndts != 1:
            self.dt = x[0] / self.ndts

        # passing the states as the initial condition
        a = self.run_filament_code(ndts, x)
        
        a = a[-1][2:]

        y = np.zeros_like(x)
        
        y[1:] = a
        
        return y

    def getrhs(self, x):
        # Function to be minimised
        y_ = self.steporbit(self.ndts, x)

        # Calculate the difference
        y = y_ - x 

        # Take into account the periodicity and remove fils near the pole
        y[1:self.NFIL+1] -= 2*np.pi

        # Remove the fils near the poles
        fil_references = util.read_fil_references(self.d.dir + self.d.simName + '_fil_references.dat')
        fil_references_sphpolar = np.zeros((self.NFIL,3))
        for i in range(self.NFIL):
            fil_references_sphpolar[i] = util.cartesian_to_spherical(fil_references[3*i: 3*i+3])

        near_pole_ind = np.where(np.sin(fil_references_sphpolar[:,2]) < 0.0 )
        # print(near_pole_ind)
        # print('-------------')
        
        # print(y[1:][near_pole_ind])
        # print(y[(1+self.NFIL):][near_pole_ind])
        # print('-------------')
        
        
        y[1:][near_pole_ind] = 0
        y[(1+self.NFIL):][near_pole_ind] = 0

        # print(y[1:][near_pole_ind])
        # print(y[(1+self.NFIL):][near_pole_ind])
        # print('-------------')

        # Set the period to 0 (constraints, rhs=0)
        y[0] = 0.0  

        return y

    def saveorbit(self, save):
        # called each Newton iteration
        print(f'[\033[32mnewton\033[m]: iteration {self.new_nits}')

        aux = np.copy(self.new_x[1:]) # using the F(x) norm instead of F(x+T)
        # detach error from norm
        aux[:self.NFIL] = np.ones(self.NFIL)*np.pi
        norm_x = np.linalg.norm(aux)

        # norm_x = np.sqrt(self.dotprd(-1, self.new_x, self.new_x))

        relative_err = self.new_tol / norm_x

        print(f'[\033[32mnewton\033[m]: T = {self.new_x[0]}')
        print(f'[\033[32mnewton\033[m]: norm x = {norm_x}, new tolerance = {self.new_tol}')
        print(f'[\033[32mnewton\033[m]: relative error: {relative_err}')

        # save_x = np.insert(self.new_x, 0, self.k)
        if (save==1):
            with open(self.d.dir + f"errors.dat", "ab") as f:
                f.write(b"\n")
                np.savetxt(f, [self.k, relative_err], newline=" ")
                
        # SAVE current solution, new_x (You can add your saving logic here)
    
    def multJp(self, x):
        # preconditioner for multJ.  Empty - no preconditioner required
        return x

    def multJ(self, dx):
        # Action of Jacobian on update dx, and evaluation of constraint on update.
        # Use approximation    dF(x_n)/dx . dx = (F(x_n+eps.dx)-F(x_n))/eps

        # print(f"////////////dx = {dx[0]}, {dx[1]}, {dx[2]}")
        # (F(x0+eps*x) - F(x0))/eps
        eps = np.sqrt(self.dotprd(-1, dx, dx)) 
        eps = self.epsJ * np.sqrt(self.dotprd(-1, self.new_x, self.new_x)) / eps
        # should eps include \del period?
        y = self.new_x + eps * dx 
        # print(f"computing F(x_i + eps * dx) with {y}")
        # print(f'[\033[34mgmres\033[m]: computing F(x_i + eps * dx) with eps={eps} T={self.new_x[0]} eps*dT={eps * dx[0]}')
        s = self.getrhs(y) # F(x_i + eps * dx)
        y = (s - self.new_fx) / eps 
    
        if self.fixT: # No extra constraint if T fixed
            y[0] = 0.0
        else: # Constrant: dx . \dot{x} = 0
            # print(f"computing \dot(x) with {self.new_x}")
            s = self.steporbit(1, self.new_x) 
            self.dt = self.new_x[0] / self.ndts
            s = (s - self.new_x) / self.dt
            # s = (s - self.new_x) / self.dt - 2*np.pi/self.new_x[0]
            y[0] = self.dotprd(-1, s, dx) # s = \dot{x}
    
        return y


    def NewtonHook(self, m, n, gtol, tol, del_value, mndl, mxdl, nits, info):
        g = gmresc.GMRES()
        # print(f"computing F(x_i) with {self.new_x}")
        self.new_x[1:self.NFIL+1] = util.box(self.new_x[1:self.NFIL+1], 2*np.pi)
        self.new_fx = self.getrhs(self.new_x)
        self.new_tol = np.sqrt(self.dotprd(-1, self.new_fx, self.new_fx))
        self.new_del = del_value
        self.dt = self.new_x[0] / self.ndts

        v = np.zeros((n, m + 1))
        mxdl_ = mxdl
        ginfo = info

        if del_value < 0.0:
            self.new_del = self.new_tol / 10.0
            mxdl_ = 1e99

        if info == 1:
            print(f'[\033[32mnewton\033[m]: nits={self.new_nits}  res={self.new_tol}')

        # stt = 0
        # print('new_x (x(0))', self.new_x[stt:stt+10])
        # print('new_fx (x(T)-x(0))', self.new_fx[stt:stt+10])

        self.saveorbit(0)
        x_ = np.copy(self.new_x)
        fx_ = np.copy(self.new_fx)
        tol_ = self.new_tol
        tol__ = 1e99

        if self.new_tol < tol:
            if info == 1:
                print('[\033[32mnewton\033[m]: input already converged')
                self.saveorbit(1)
            info = 0
            return info

        # main loop:
        while True:
            if self.new_del < mndl:
                if info == 1:
                    print('[\033[32mnewton\033[m]: trust region too small')
                info = 3
                self.saveorbit(1)
                return info

            # Find hookstep s and update x
            s = np.zeros(n)
            gres = gtol * self.new_tol
            gdel = self.new_del

            if ginfo != 2:
                self.new_gits = m

            if del_value == 0.0:
                self.new_gits = 9999

            s, gres, gdel, self.new_gits, ginfo = g.GMRESm(m, n, s, fx_, self.multJ, self.multJp, self.dotprd, gres, gdel, self.new_gits, ginfo)
            # print(f'[\033[35mDEBUG\033[m]: del x from GMRES is {-s}')
            ginfo = info
            self.new_x = x_ - s
            self.new_x[1:self.NFIL+1] = util.box(self.new_x[1:self.NFIL+1], 2*np.pi)
            # print(f'[\033[35mDEBUG\033[m]: new T is {self.new_x[0]}')

            # print(f"computing F(x_i) with {self.new_x}")
            self.new_fx = self.getrhs(self.new_x)
            self.new_tol = np.sqrt(self.dotprd(-1, self.new_fx, self.new_fx))
            snrm = np.sqrt(self.dotprd(-1, s, s))
            ared = tol_ - self.new_tol
            pred = tol_ - gdel

            if info == 1:
                print(f'[\033[32mnewton\033[m]: nits={self.new_nits}  res={self.new_tol}')
                print(f'[\033[32mnewton\033[m]: gits={self.new_gits}  del={self.new_del}')
                print(f'[\033[32mnewton\033[m]: |s|={snrm}  pred={pred}')
                print(f'[\033[32mnewton\033[m]: ared/pred={ared/pred}')

            if del_value == 0.0:
                if info == 1:
                    print('[\033[32mnewton\033[m]: took full newton step')

            elif self.new_tol > tol__:
                if info == 1:
                    print('[\033[32mnewton\033[m]: accepting the previous step')
                self.new_x = x__
                self.new_fx = fx__
                self.new_tol = tol__
                self.new_del = del__
            elif ared < 0.0:
                if info == 1:
                    print('[\033[32mnewton\033[m]: norm increased, try a smaller step')
                self.new_del = snrm * 0.5
                ginfo = 2

                if ared/pred < -50:
                    print('[\033[32mnewton\033[m]: no valid update found')
                    info = 3
                    self.saveorbit(1)
                    return info
            
            elif ared / pred < 0.75:
                if info == 1:
                    print('[\033[32mnewton\033[m]: step is okay, trying a smaller step')
                x__ = self.new_x
                fx__ = self.new_fx
                tol__ = self.new_tol
                if ared / pred > 0.1:
                    del__ = snrm
                if ared / pred <= 0.1:
                    del__ = snrm * 0.5
                self.new_del = snrm * 0.7
                ginfo = 2
            elif snrm < self.new_del * 0.9:
                if info == 1:
                    print('[\033[32mnewton\033[m]: step is good, took a full newton step')
                self.new_del = min(mxdl_, snrm * 2.0)
            elif self.new_del < mxdl_ * 0.9:
                if info == 1:
                    print('[\033[32mnewton\033[m]: step is good, trying a larger step')
                x__ = self.new_x
                fx__ = self.new_fx
                tol__ = self.new_tol
                del__ = self.new_del
                self.new_del = min(mxdl_, snrm * 2.0)
                ginfo = 2

            # Check if need to try another s
            if ginfo == 2:
                continue

            # End of iteration
            self.saveorbit(0)
            self.new_nits += 1
            x_ = self.new_x
            fx_ = self.new_fx
            tol_ = self.new_tol
            tol__ = 1e99

            if self.new_tol < tol:
                if info == 1:
                    print('[\033[32mnewton\033[m]: converged')
                info = 0
                self.saveorbit(1)
                return info
            elif self.new_nits == nits:
                if info == 1:
                    print('[\033[32mnewton\033[m]: reached max its')
                info = 2
                return info
