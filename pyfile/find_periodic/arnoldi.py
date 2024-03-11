import numpy as np
import numpy.linalg as la
import os
import time

class ARNOLDI:

    def __init__(self, NSEG, NFIL, NBLOB, AR, T, Ustar, Num_evals, ct_var, d, evalf, evecf):
        # Num_evals can be assigned as a memeber
        self.Num_evals = Num_evals
        self.NSEG = NSEG
        self.NFIL = NFIL
        self.NBLOB = NBLOB
        self.AR = AR
        self.NTOTAL = 2*NFIL
        self.T = T
        self.Q = np.zeros((self.NTOTAL,self.NTOTAL))
        self.H = np.zeros((self.NTOTAL,self.NTOTAL))
        self.epsilon = 5e-1
        self.tol = 1e-3
        self.old_evalue = 10*np.ones(self.Num_evals,dtype = np.complex )
        self.difference = 2*self.tol
        self.Ustar = Ustar
        self.ct_var = ct_var

        self.d = d

        self.evalf = evalf
        self.evecf = evecf
    

    def generate_initial_condition(self):
        
        b = np.random.rand(self.NTOTAL)-0.5
        b[:self.NFIL] *= 50
        # b[self.NFIL:] = 0
        
        b = b/la.norm(b)
        print(b[:10]*self.epsilon)
        print(self.Ustar[:10])
        self.Q[:,0] = b
        
            
    def save_initial_condition(self, k):

        input_filename = self.d.dir + "psi.dat"

        initial_condition = self.Ustar + self.epsilon * self.Q[:,k]
        x = np.insert( initial_condition, 0, [self.ct_var, self.T ])

        np.savetxt(input_filename, x, newline = " ")

    def Ax(self):

        self.d.change_variables(self.NFIL, self.NSEG, self.NBLOB, self.AR, self.ct_var, self.T, 1.)
        self.d.update_globals_file()
        self.d.run()

        # Read output data
        output_filename = self.d.dir + self.d.simName + "_true_states.dat"
        U = np.loadtxt(output_filename)[-1][2:]
        U[:self.NFIL] -= 2*np.pi

        # print(U - self.Ustar)
        print(f'norm of dU(T) = {la.norm(U-self.Ustar)}, eps*dU = {la.norm(U-self.Ustar)/self.epsilon}')
        print((U-self.Ustar)[:10])
        

        # Normalise data
        return (U - self.Ustar)/self.epsilon


    def gramschmidt_iteration(self,k,U):
        
        # Carry out Gram-Schmidt on U against Q
        for i in range(k+1):
            qi = self.Q[:, i]
            self.H[i,k] = qi @ U
            U = U - self.H[i,k]*qi

        if k+1 < self.NTOTAL:
            self.H[k+1, k] = la.norm(U)
            self.Q[:, k+1] = U/self.H[k+1, k]

    def log_over_T_and_sort(self,evals,T):
        evals = np.log(np.array(evals,dtype = np.complex_))/T
        idx = evals.argsort()[::-1]   
        evals = evals[idx]

        return evals, idx

    def find_difference(self, k, T):

        new_evals = la.eig(self.H[0:k+1, 0:k+1])[0]

        new_evals, idx = self.log_over_T_and_sort(new_evals,T)

        temp = min(k+1, self.Num_evals)
        print(f'\033[33m[Old eval\033[m] = {self.old_evalue[0:temp]}')
        print(f'\033[33m[New eval\033[m] = {new_evals[0:temp]}')
        
        self.difference = la.norm(np.subtract(new_evals[0:temp], self.old_evalue[0:temp]))

        self.old_evalue[0:temp] = new_evals[0:temp]
        # print("Max eigenvalues" + str(self.old_evalue))


    def save_evals_and_evecs_to_file(self, k, T):

        eigenValues, eigenVectors = la.eig(self.H[0:k,0:k])

        eigenValues, idx = self.log_over_T_and_sort(eigenValues,T)
        eigenValues = np.append([self.ct_var],eigenValues)

        eigenVectors = self.Q[:,0:k] @ eigenVectors[:,idx]

        with open(self.evalf, "ab") as f:
            f.write(b"\n")
            np.savetxt(f, eigenValues[0:self.Num_evals+1], newline = " ")
            

        with open(self.evecf, "ab") as f:
            f.write(b"\n")
            np.savetxt(f, eigenVectors[:,0:self.Num_evals], newline = " ")
            

    def arnoldi_for_eigenvalues(self, T):
        self.generate_initial_condition()
        k = 0

        while np.absolute(self.difference) > self.tol or k < self.Num_evals:

            self.save_initial_condition(k)

            U = self.Ax()

            self.gramschmidt_iteration(k,U)

            self.find_difference(k, T)

            print(f"[\033[32mArnoldi\033[m]Difference at step " + str(k) + " is " + str(self.difference))
            k += 1    

        self.save_evals_and_evecs_to_file(k, T)
