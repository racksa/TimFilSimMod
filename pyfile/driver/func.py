import configparser
import os
import math
import util

class DRIVER:

    def __init__(self):
        self.globals_name = 'globals.ini'
        self.dir = "data/expr_sims/20230814/"
        self.pars_list = {"nfil": [],
                     "nblob": [],
                     "ar": [],
                     "spring_factor": []}
        
        self.sweep_shape = (8, 8, 1, 1)

        self.num_sim = 0

        self.current_thread = 0
        self.num_thread = 1
        self.cuda_device = 0
    
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

    def create_rules(self):
        # Define the rule of sweeping simulations
        for i in range(self.sweep_shape[0]):
            for j in range(self.sweep_shape[1]):
                for k in range(self.sweep_shape[2]):
                    for l in range(self.sweep_shape[3]):
                        nfil = int(32*(i+1) )
                        nblob = int(2000*(1.2**j))
                        ar = round(3*(1.2**j), 2)
                        spring_factor = round(2, 2)

                        self.pars_list["nfil"].append(nfil)
                        self.pars_list["nblob"].append(nblob)
                        self.pars_list["ar"].append(ar)
                        self.pars_list["spring_factor"].append(spring_factor)
        # Write rules to sim list file
        self.write_rules()

    def delete_files(self):
        util.delete_files_in_directory(self.dir)

    def write_rules(self):
        sim = configparser.ConfigParser()
        sim.add_section('Parameter list')
        for key, value in self.pars_list.items():
            sim['Parameter list'][key] = ', '.join(map(str, value))
        with open(self.dir+"rules.ini", 'w') as configfile:
            sim.write(configfile, space_around_delimiters=False)

    def read_rules(self):
        sim = configparser.ConfigParser()
        sim.read(self.dir+"rules.ini")
        for key, value in self.pars_list.items():
            self.pars_list[key] = [float(x) for x in sim["Parameter list"][key].split(', ')]
        self.num_sim = len(self.pars_list["nfil"])

    def run(self):
        self.create_ini()
        self.write_ini("Filenames", "simulation_dir", self.dir)

        # Read rules from the sim list file
        self.read_rules()

        thread_list = util.even_list_index(self.num_sim, self.num_thread)
        sim_index_start = thread_list[self.current_thread]
        sim_index_end = thread_list[self.current_thread+1]

        print(f"Partitioning {self.num_sim} into {self.num_thread} threads\n" +\
              f"Partition index: {self.current_thread} / {self.num_thread-1} \n" + \
              f"[{sim_index_start} - {sim_index_end}] / {thread_list}\n" +\
              f"on GPU: {self.cuda_device}")

        # Iterate through the sim list and write to .ini file and execute
        for i in range(sim_index_start, sim_index_end):
            for key, value in self.pars_list.items():
                self.write_ini("Parameters", key, float(self.pars_list[key][i]))
            self.write_ini("Filenames", "simulation_file", f"ciliate_{self.pars_list['nfil'][i]:.0f}fil_{self.pars_list['nblob'][i]:.0f}blob_{self.pars_list['ar'][i]:.2f}R_{self.pars_list['spring_factor'][i]:.2f}torsion")
            self.write_ini("Filenames", "simulation_dir", self.dir)

            command = f"export OPENBLAS_NUM_THREADS=1; \
                        export CUDA_VISIBLE_DEVICES={self.cuda_device}; \
                        ./bin/cilia > terminal_outputs/output_{self.pars_list['nfil'][i]:.0f}fil_{i}.out"

                # command = f"export OPENBLAS_NUM_THREADS=1; \
                #             export CUDA_VISIBLE_DEVICES={self.cuda_device}; \
                #             ./bin/cilia"
            os.system(command)