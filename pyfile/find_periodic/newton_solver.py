import numpy as np
import gmres
import driver
import os
import subprocess
import time

class NEWTON_SOLVER:
    
    def __init__(self, new_x, epsJ, ndts, fixT, k, NSEG, NFIL, y_spacing_const, x_spacing_const):
        self.new_x = new_x
        self.ndts = ndts
        self.k = k
        self.NSEG = NSEG
        self.NFIL = NFIL
        self.y_spacing_const = y_spacing_const 
        self.x_spacing_const = x_spacing_const
        self.epsJ = epsJ
        self.fixT = fixT
        self.new_fx = []
        self.new_tol = []
        self.new_del = []
        self.new_nits = 0
        self.new_gits = 0
        self.dt = []
        self.period = 1.0

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
    
    def save_configuration_to_file(self,initial_condition):
         
        simulation_dir = "data/find_periodic/" 
        input_filename = simulation_dir + "input_liealgebra_configuration.dat"

        np.savetxt(input_filename, initial_condition, newline = " ")

        return

    def change_variables(self,ndts):

        # These are always changed at the beginning of each Arnoldi iteration 
        self.d.change_variables(self.k, self.period)
        simulation_file = self.d.update_globals_file()

        return simulation_file
    
    def read_data_from_file(self,sim_dir_and_name):
        
        output_filename = sim_dir_and_name + "_true_states.dat"
        return np.loadtxt(output_filename)

    def run_filament_code(self,a,ndts):

        # Save Lie alg elements
        self.save_configuration_to_file(a)

        # Debugging
        # filename = "data/find_periodic/saving_intermidate_liealg"
        # with open(filename, "ab") as f:
        #     np.savetxt(f, a, newline = " ")
        #     f.write(b"\n")
        
        # Change globals.ini file
        sim_name = self.change_variables(ndts)

        # Run code
        # os.chdir('../filament')
        # command = "./run_cilia_sim.sh"

        # subprocess.run(['sh', './run_cilia_sim.sh'])
        #os.system(command)

        self.d.run()

        # Return output
        return self.read_data_from_file(self.d.dir + sim_name)


    def steporbit(self, ndts, x):
        
        if ndts != 1:
            self.dt = x[0] / self.ndts

        # x[1:] is the phases
        a = self.run_filament_code(x[1:],ndts)
        
        a = a[-1][2:]

        y = np.zeros_like(x)
        

        y[1:] = a
        
        return y
        

    def getrhs(self, n_, x):
        # function to be minimised
        y_ = self.steporbit(self.ndts, x)
        y_[:639] = (y_[:639] - np.floor(y_[:639]/(2*np.pi))*(2*np.pi)) 
        x[:639] = (x[:639] - np.floor(x[:639]/(2*np.pi))*(2*np.pi))

        # print('y', y_)
        # print('x', x)

        y = y_ - x # Calculate the difference
        y[0] = 0.0  # Set the first element to 0 (constraints, rhs=0)
        return y

    def saveorbit(self):
        # called each Newton iteration
        print(f'newton: iteration {self.new_nits}')

        norm_x = np.sqrt(self.dotprd(-1, self.new_x, self.new_x))
        if self.fixT == 1:
            relative_err = self.new_tol 
        else:
            relative_err = self.new_tol / norm_x

        print(f'norm x = {norm_x}, new tolerance = {self.new_tol}')
        print(f'relative error: {relative_err}')
#        time.sleep(1)
        
        self.d.save_orbit()

        # SAVE current solution, new_x (You can add your saving logic here)
    
    def multJp(self, n, x):
        # preconditioner for multJ.  Empty - no preconditioner required
        return x

    def multJ(self, n_, dx):
        # Action of Jacobian on update dx, and evaluation of constraint on update.
        # Use approximation    dF(x_n)/dx . dx = (F(x_n+eps.dx)-F(x_n))/eps

        # (F(x0+eps*x) - F(x0))/eps
        eps = np.sqrt(self.dotprd(1, dx, dx))
        eps = self.epsJ * np.sqrt(self.dotprd(1,self.new_x, self.new_x)) / eps
        y = self.new_x + eps * dx
        s = self.getrhs(n_, y)
        y = (s - self.new_fx) / eps
    
        if self.fixT: # No extra constraint if T fixed
            y[0] = 0.0
        else: # Constrant: dx . \dot{x} = 0
            s = self.steporbit(1,self.new_x)
            self.dt = self.new_x[0] / self.ndts
            s = (s - self.new_x) / self.dt
            y[0] = self.dotprd(-1, s, dx)
    
        return y


    def NewtonHook(self, m, n, gtol, tol, del_value, mndl, mxdl, nits, info):
        g = gmres.GMRES()
        self.new_fx = self.getrhs(n, self.new_x)
        self.new_tol= np.sqrt(self.dotprd(n, self.new_fx, self.new_fx))
        self.new_del = del_value
        self.dt = self.new_x[0] / self.ndts

        v = np.zeros((n, m + 1))
        mxdl_ = mxdl
        ginfo = info

        if del_value < 0.0:
            self.new_del = self.new_tol / 10.0
            mxdl_ = 1e99

        if info == 1:
            print(f'newton: nits={self.new_nits}  res={self.new_tol}')

        self.saveorbit()
        x_ = np.copy(self.new_x)
        fx_ = np.copy(self.new_fx)
        tol_ = self.new_tol
        tol__ = 1e99

        if self.new_tol < tol:
            if info == 1:
                print('newton: input already converged')
            info = 0
            return info

        # main loop:
        while True:
            if self.new_del < mndl:
                if info == 1:
                    print('newton: trust region too small')
                info = 3
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
            ginfo = info
            self.new_x = x_ - s

            self.new_fx = self.getrhs(n, self.new_x)
            self.new_tol = np.sqrt(self.dotprd(n, self.new_fx, self.new_fx))
            snrm = np.sqrt(self.dotprd(n, s, s))
            ared = tol_ - self.new_tol
            pred = tol_ - gdel

            if info == 1:
                print(f'newton: nits={self.new_nits}  res={self.new_tol}')
                print(f'newton: gits={self.new_gits}  del={self.new_del}')
                print(f'newton: |s|={snrm}  pred={pred}')
                print(f'newton: ared/pred={ared/pred}')

            if del_value == 0.0:
                if info == 1:
                    print('newton: took full newton step')

            elif self.new_tol > tol__:
                if info == 1:
                    print('newton: accepting the previous step')
                self.new_x = x__
                self.new_fx = fx__
                self.new_tol = tol__
                self.new_del = del__
            elif ared < 0.0:
                if info == 1:
                    print('newton: norm increased, try a smaller step')
                self.new_del = snrm * 0.5
                ginfo = 2
            elif ared / pred < 0.75:
                if info == 1:
                    print('newton: step is okay, trying a smaller step')
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
                    print('newton: step is good, took a full newton step')
                self.new_del = min(mxdl_, snrm * 2.0)
            elif self.new_del < mxdl_ * 0.9:
                if info == 1:
                    print('newton: step is good, trying a larger step')
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
            self.new_nits += 1
            self.saveorbit()
            x_ = self.new_x
            fx_ = self.new_fx
            tol_ = self.new_tol
            tol__ = 1e99

            if self.new_tol < tol:
                if info == 1:
                    print('newton: converged')
                info = 0
                return info
            elif self.new_nits == nits:
                if info == 1:
                    print('newton: reached max its')
                info = 2
                return info
