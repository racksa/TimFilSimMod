import numpy as np
import gmres
import driver
import os
import subprocess
import time
import util

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
    
    def run_filament_code(self, x):

        # Save states to a file to read and run
        x = np.insert(x, 0, self.k)
        # x[2:self.NFIL+2] = np.arcsin(x[2:self.NFIL+2]) + 2*np.pi
        # x[2:self.NFIL+2] = np.angle(x[2:self.NFIL+2])
        x[2:self.NFIL+2] = util.box(x[2:self.NFIL+2], 2*np.pi)
        np.savetxt(self.d.dir + "psi.dat", x, newline = " ")

        # Change globals.ini file
        self.d.change_variables(self.NFIL, self.NSEG, self.NBLOB, self.AR, self.k, x[1], 1.0)
        sim_name = self.d.update_globals_file()

        # Run code
        print(f"Running new sim with T={x[1]}")
        self.d.run()

        # Return output
        output_filename = self.d.dir + sim_name + "_true_states.dat"
        return np.loadtxt(output_filename)

    def steporbit(self, ndts, x):
        
        if ndts != 1:
            self.dt = x[0] / self.ndts

        # x[1:] is the phases
        a = self.run_filament_code(x)
        
        a = a[-1][2:]

        y = np.zeros_like(x)
        
        y[1:] = a
        
        return y
        

    def getrhs(self, n_, x):
        # function to be minimised
        y_ = self.steporbit(self.ndts, x)
        # y_[1:self.NFIL+1] = np.sin(y_[1:self.NFIL+1])
        # y_[1:self.NFIL+1] = np.exp(1j*y_[1:self.NFIL+1])
        # y_[1:self.NFIL+1] = util.box(y_[1:self.NFIL+1], 2*np.pi)

        y = y_ - x # Calculate the difference
        y[1:self.NFIL+1] -= 2*np.pi # Take into account the periodicity of the phase

        y[0] = 0.0  # Set the first element to 0 (constraints, rhs=0)

        # stt = 0
        # print('x(T)', y_[stt:stt+10])
        # print('x(0)', x[stt:stt+10])
        # print('x(T)-x(0)', y[stt:stt+10])
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
        self.new_tol = np.sqrt(self.dotprd(n, self.new_fx, self.new_fx))
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

        # stt = 0
        # print('new_x (x(0))', self.new_x[stt:stt+10])
        # print('new_fx (x(T)-x(0))', self.new_fx[stt:stt+10])

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