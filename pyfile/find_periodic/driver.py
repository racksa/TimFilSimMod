import configparser
import subprocess
import os
import math
import numpy as np

class DRIVER:

    def __init__(self):
        self.globals_name = 'globals.ini'
        self.exe_name = 'cilia_periodic'
        self.date = '20240118_periodic'
        self.afix = ''
        self.dir = f"data/expr_sims/{self.date}{self.afix}/"
        self.pars_list = {
                     "nswim": [],
                     "nseg": [],
                     "nfil": [],
                     "nblob": [],
                     "ar": [],
                     "spring_factor": [],
                     "force_mag": [],
                     "seg_sep": [],
                     "period": [],
                     "sim_length": []}
        
        self.current_thread = 0
        self.num_thread = 1
        self.cuda_device = 1
    
    def create_ini(self):
        ini = configparser.ConfigParser()
        ini.add_section('Parameters')
        ini.add_section('Filenames')
        with open(self.globals_name, 'w') as configfile:
            ini.write(configfile, space_around_delimiters=False)
        

    def write_ini(self, section, variable, value):
        ini = configparser.ConfigParser()
        ini.read(self.globals_name)

        ini.set(section, variable, f'{value}')

        # Save the changes back to the file
        with open(self.globals_name, 'w') as configfile:
            ini.write(configfile, space_around_delimiters=False)


    def change_variables(self, nfil, nseg, nblob, ar, k, period, sim_length):
        nseg = nseg
        nfil = nfil
        nblob = nblob
        ar = ar
        sim_length = sim_length
        force_mag = 1
        seg_sep = 2.6
        
        self.pars_list["nswim"].append(1)
        self.pars_list["nseg"].append(nseg)
        self.pars_list["nfil"].append(nfil)
        self.pars_list["nblob"].append(nblob)
        self.pars_list["ar"].append(ar)
        self.pars_list["spring_factor"].append(k)
        self.pars_list["force_mag"].append(force_mag)
        self.pars_list["seg_sep"].append(seg_sep)
        self.pars_list["period"].append(period)
        self.pars_list["sim_length"].append(sim_length)


    def update_globals_file(self):
        self.create_ini()

        readphase_index = ''
        # Iterate through the sim list and write to .ini file and execute
        for key, value in self.pars_list.items():
            self.write_ini("Parameters", key, float(self.pars_list[key][0]))
        self.simName = f"ciliate_{self.pars_list['nfil'][0]:.0f}fil_{self.pars_list['nblob'][0]:.0f}blob_{self.pars_list['ar'][0]:.2f}R_{self.pars_list['spring_factor'][0]:.3f}torsion"
        simulation_file = self.simName
        self.write_ini("Filenames", "simulation_file", self.simName)
        self.write_ini("Filenames", "simulation_dir", self.dir)
        self.write_ini("Filenames", "simulation_readphase_name", f"phases.dat")
        self.write_ini("Filenames", "simulation_readangle_name", f"angles.dat")

        return simulation_file
    
    def save_orbit(self):
        # print("Saving orbits using the driver....")
        input_filenames = self.simName + '_true_states.dat'
        
        input_filename = self.dir + input_filenames

        # Open the input file in read mode
        with open(input_filename, 'r') as input_file:
            # Read all lines from the file
            lines = input_file.readlines()

            # Check if the file is not empty
            if lines:
                # Extract the last line
                first_line = lines[0]
                last_line = lines[-1]
                true_states_start = np.array(first_line.split(), dtype=float)
                true_states = np.array(last_line.split(), dtype=float)

                np.savetxt(self.dir + f"phases.dat", true_states[2:self.pars_list["nfil"][0]+2], delimiter=' ', newline=' ')
                np.savetxt(self.dir + f"angles.dat", true_states[self.pars_list["nfil"][0]+2:], delimiter=' ', newline=' ')
                
                # np.savetxt(self.dir + f"psi.dat", true_states, delimiter=' ', newline=' ')
                # np.savetxt(self.dir + f"psi_start.dat", true_states_start, delimiter=' ', newline=' ')

                # print(f"[SUCCESS]: last line copied from '{input_filename}' to '{output_filenames[0]}'.")
            else:
                print(f"The file '{input_filename}' is empty.")

    def run(self):
        
        # for key, value in self.pars_list.items():
        #     self.write_ini("Parameters", key, float(self.pars_list[key][0]))
        # self.simName = f"ciliate_{self.pars_list['nfil'][0]:.0f}fil_{self.pars_list['nblob'][0]:.0f}blob_{self.pars_list['ar'][0]:.2f}R_{self.pars_list['spring_factor'][0]:.3f}torsion"
        # self.write_ini("Filenames", "simulation_file", self.simName)
        # self.write_ini("Filenames", "simulation_dir", self.dir)

        # print(self.pars_list)

        command = f"export OPENBLAS_NUM_THREADS=1; \
                    export CUDA_VISIBLE_DEVICES={self.cuda_device}; \
                    ./bin/{self.exe_name} > terminal_outputs/output_{self.date}_{self.pars_list['nfil'][0]:.0f}fil_{self.pars_list['sim_length'][0]:.0f}period_{0}.out"
        # subprocess.run(command)
        os.system(command)