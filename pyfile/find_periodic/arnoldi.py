import numpy as np
import numpy.linalg as la
import os
import time

class ARNOLDI:

    def __init__(self,NSEG,NFIL,Ustar,Num_evals):
        NTOTAL = 3*(NSEG-1)*NFIL
        self.NTOTAL = NTOTAL
        self.Q = np.zeros((NTOTAL,NTOTAL))
        self.H = np.zeros((NTOTAL,NTOTAL))
        self.epsilon = 1e-3
        self.tol = 1e-8
        self.old_evalue = 10*np.ones(Num_evals,dtype = np.complex )
        self.difference = 2*self.tol
        self.Ustar = Ustar
    

    def generate_initial_condition(self,sim_3D):
        
        if sim_3D:
            b = np.random.rand(self.NTOTAL)
        else:
            b = np.zeros((int(self.NTOTAL/3),3))
            b[:,1] = np.random.rand(int(self.NTOTAL/3))
            b = np.reshape(b,self.NTOTAL)
        
        b = b/la.norm(b)
        self.Q[:,0] = b
        
        del b
            
    def save_initial_condition(self,k):

        simulation_dir = "../filament/" 
        input_filename = simulation_dir + "input_liealgebra_configuration.dat"

        initial_condition = self.Ustar + self.epsilon * self.Q[:,k]

        np.savetxt(input_filename, initial_condition, newline = " ")

    def Ax(self,sim_dir,sim_name):

        os.chdir('../filament')
        command = f"./run_cilia_sim.sh"
        os.system(command)

        # Read output data
        output_filename = sim_dir + sim_name + "_liealgebra_configuration.dat"
        U = np.loadtxt(output_filename)
        
        # Normalise data
        return np.subtract(U,self.Ustar)/self.epsilon


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

    def find_difference(self, k, Num_evals,T):

        new_evals = la.eig(self.H[0:k+1,0:k+1])[0]

        new_evals, idx = self.log_over_T_and_sort(new_evals,T)

        temp = min(k+1,Num_evals)
        print("New eval: " + str(new_evals[0]) + ", old eval: " + str(self.old_evalue[0]))

        self.difference = la.norm(np.subtract(new_evals[0:temp],self.old_evalue[0:temp]))
 #       self.difference = la.norm(np.subtract(new_evals[0],self.old_evalue[0]))
        #TODO: change to check first Num_evals eigenvalues

        self.old_evalue[0:temp] = new_evals[0:temp]
        print("Max eigenvalues" + str(self.old_evalue))


    def save_evals_and_evecs_to_file(self, k, T, Num_evals, varied_value, eval_filename, evec_filename):

        eigenValues, eigenVectors = la.eig(self.H[0:k,0:k])

        eigenValues, idx = self.log_over_T_and_sort(eigenValues,T)
        eigenValues = np.append([varied_value],eigenValues)

        eigenVectors = self.Q[:,0:k] @ eigenVectors[:,idx]

        os.chdir('../stability')

        with open(eval_filename, "ab") as f:
            np.savetxt(f, eigenValues[0:Num_evals+1], newline = " ")
            f.write(b"\n")

        with open(evec_filename, "ab") as f:
            np.savetxt(f, eigenVectors[:,0:Num_evals], newline = " ")
            f.write(b"\n")

    def arnoldi_for_eigenvalues(self, sim_3D, Num_evals, sim_name, sim_dir, T, varied_value, eval_filename, evec_filename):
        self.generate_initial_condition(sim_3D)
        k = 0

        while np.absolute(self.difference) > y  v or k < Num_evals:

            self.save_initial_condition(k)

            U = self.Ax(sim_dir,sim_name)
            self.gramschmidt_iteration(k,U)

            self.find_difference(k, Num_evals,T)

            print("Difference at step " + str(k) + " is " + str(self.difference))
            k += 1    

        self.save_evals_and_evecs_to_file(k, T, Num_evals, varied_value, eval_filename, evec_filename)
