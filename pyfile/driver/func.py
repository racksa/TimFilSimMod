import configparser
import os
import math
import util

class DRIVER:

    def __init__(self):
        self.globals_name = 'globals.ini'
        self.dir = "data/expr_sims/20230731/"
        self.pars = {"nfil": [],
                     "nblob": [],
                     "ar": [],
                     "spring_factor": []}
        
        self.sweep_shape = (12, 1, 1, 1)

        self.num_sim = 0

        self.current_thread = 0
        self.num_thread = 1
        self.cuda_device = 0

    def write_ini(self, section, variable, value):
        ini = configparser.ConfigParser()
        ini.read('globals.ini')

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
                        nfil = int(128)
                        nblob = int(1000*1.1**i)
                        ar = round(3*1.1**i)
                        spring_factor = round(2)

                        self.pars["nfil"].append(nfil)
                        self.pars["nblob"].append(nblob)
                        self.pars["ar"].append(ar)
                        self.pars["spring_factor"].append(spring_factor)
        # Write rules to sim list file
        self.write_rules()

    def delete_files(self):
        util.delete_files_in_directory(self.dir)

    def write_rules(self):
        sim = configparser.ConfigParser()
        sim.add_section('Parameter list')
        for key, value in self.pars.items():
            sim['Parameter list'][key] = ', '.join(map(str, value))
        with open(self.dir+"rules.ini", 'w') as configfile:
            sim.write(configfile, space_around_delimiters=False)

    def read_rules(self):
        sim = configparser.ConfigParser()
        sim.read(self.dir+"rules.ini")
        for key, value in self.pars.items():
            self.pars[key] = [float(x) for x in sim["Parameter list"][key].split(', ')]

    def run(self):
        
        self.write_ini("Filenames", "simulation_dir", self.dir)

        # Read rules from the sim list file
        self.read_rules()
        self.num_sim = len(self.pars["nfil"])

        thread_list = util.even_list_index(self.num_sim, self.num_thread)

        sim_index_start = thread_list[self.current_thread]
        sim_index_end = thread_list[self.current_thread+1]

        print(f"Partitioning {self.num_sim} into {self.num_thread} threads\n" +\
              f"Partition index: {self.current_thread} / {self.num_thread-1} \n" + \
              f"[{sim_index_start} - {sim_index_end}] / {thread_list}\n" +\
              f"on GPU: {self.cuda_device}")

        # Iterate through the sim list and write to .ini file and execute
        for i in range(sim_index_start, sim_index_end):
            for key, value in self.pars.items():
                self.write_ini("Parameters", key, float(self.pars[key][i]))
            self.write_ini("Filenames", "simulation_file", f"ciliate_{self.pars['nfil'][i]:.0f}fil_{self.pars['nblob'][i]:.0f}blob_{self.pars['ar'][i]:.2f}R_{self.pars['spring_factor'][i]:.2f}torsion")

            command = f"export OPENBLAS_NUM_THREADS=1; \
                        export CUDA_VISIBLE_DEVICES={self.cuda_device}; \
                        ./bin/cilia"
            os.system(command)