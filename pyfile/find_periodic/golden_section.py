import math
import driver
import numpy as np
import util


invphi = (math.sqrt(5) - 1) / 2  # 1 / phi
invphi2 = (3 - math.sqrt(5)) / 2  # 1 / phi^2

def main():
    tol = 1e-4

    NSEG = 20      # Number of segments
    NFIL = 159       # Number of filaments
    NBLOB = 9000
    AR = 8

    d = driver.DRIVER()
    d.cuda_device = 2
    d.date = 'gss'
    d.category = 'JFNK_sims/'
    d.dir = f"data/{d.category}{d.date}/"

    output_filename = d.dir + f"psi_guess{NFIL}.dat"
    with open(output_filename, 'r') as file:
        num_lines = sum(1 for line in file)
        solns = np.loadtxt(output_filename)

    num_lines = 1
    for i in range(num_lines):
        k = solns[i,0]
        T = solns[i,1]
        T0 = 0.97
        T3 = 1.
        Ustar = solns[i,2:]
        Ustar[:NFIL] = util.box(Ustar[:NFIL], 2*np.pi)
        print(f'k={k} T={T}')

        wrapper = gss_wrapper(NFIL, NSEG, NBLOB, AR, k, Ustar,\
                              d, T0, T3, tol,\
                                )
        wrapper.gss()


class gss_wrapper:

    def __init__(self, NFIL, NSEG, NBLOB, AR, ct_var, Ustar, \
                 driver_, T0, T3, tol=1e-5,\
                 ):

        self.driver_ = driver_
        self.T0 = T0
        self.T3 = T3
        self.tol = tol

        self.NFIL = NFIL
        self.NSEG = NSEG
        self.NBLOB = NBLOB
        self.AR = AR
        self.ct_var = ct_var
        self.Ustar = Ustar
        
    def gss(self):
        """Golden-section search.
        """
        a = self.T0
        b = self.T3
        tol = self.tol
        h = b - a

        # Required steps to achieve tolerance
        n = int(math.ceil(math.log(tol / h) / math.log(invphi)))
        print(f"\033[36mRequire {n} steps for tol={tol}\033[m")

        c = a + invphi2 * h
        d = a + invphi * h
        yc = self.f(c)
        yd = self.f(d)

        for k in range(n - 1):
            
            if yc < yd:  # yc > yd to find the maximum
                b = d
                d = c
                yd = yc
                h = invphi * h
                c = a + invphi2 * h
                yc = self.f(c)
            else:
                a = c
                c = d
                yc = yd
                h = invphi * h
                d = a + invphi * h
                yd = self.f(d)

            print(f"\033[37mIteration {k}: error = yc\033[m")
            with open(self.driver_.dir + "Ts.dat", "ab") as fi:
                fi.write(b"\n")
                np.savetxt(fi, [a, d, yc, yd], newline = " ")

        if yc < yd:
            return (a, d)
        else:
            return (c, b)
        
    def f(self, T):
        input_filename = self.driver_.dir + "psi.dat"
        x = np.insert( self.Ustar, 0, [self.ct_var, T ])

        np.savetxt(input_filename, x, newline = " ")

        self.driver_.change_variables(self.NFIL, self.NSEG, self.NBLOB, self.AR, self.ct_var, T, 1.)
        self.driver_.update_globals_file()
        
        self.driver_.run()

        output_filename = self.driver_.dir + self.driver_.simName + "_true_states.dat"
        U = np.loadtxt(output_filename)[-1][2:]
        U[:self.NFIL] -= 2*np.pi

        return np.linalg.norm((U - self.Ustar))

main()