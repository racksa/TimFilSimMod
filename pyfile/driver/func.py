import configparser
import os
import math
import util

class DRIVER:

    def __init__(self):
        self.globals_name = 'globals.ini'
        self.exe_name = 'cilia_readphase'
        self.date = '20240104_readphase_hold'
        self.date = '20240112_readphase_free'
        self.afix = ''
        self.dir = f"data/expr_sims/{self.date}{self.afix}/"
        # self.dir = f"data/expr_sims/{self.date}{self.afix}/"
        self.pars_list = {
                     "nswim": [],
                     "nseg": [],
                     "nfil": [],
                     "nblob": [],
                     "ar": [],
                     "spring_factor": [],
                     "force_mag": [],
                     "seg_sep": []}

        # self.sweep_shape = (3, 8, 6, 1)
        self.sweep_shape = (1, 12, 4, 1)
        self.sweep_shape = (24, 1, 1, 1)

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

                        seg_sep = 2.6
                        nseg = 20
                        force_mag = 1

                        # fil_density = 3.0 - 0.6*k
                        # blob_density = 75
                        # nfil = int( 192 + 96*j )
                        # ar = round( (nfil/fil_density)**.5 ,2)
                        # nblob = int(ar**2*blob_density)
                        # spring_factor = round(0.5+ 0.25*i, 2)

                        # # planar hexagonal
                        # nfil = int(128)
                        # nblob = int(12800)
                        # ar = round(12.65 + 0.01*i, 2)
                        # spring_factor = round(0.005 + 0.00*i, 3)

                        # k-means
                        nfil = int(639 + 0*i)
                        nblob = int(40961 + 0*i)
                        ar = round(15.00, 2)
                        spring_factor = round(0.008 + 0.002*i*(i//6+1), 3)

                        # # icosahedral
                        # nfil = int(640)
                        # nblob = int(40962)
                        # ar = round(15.0, 2)
                        # spring_factor = round(0.005 + 0.008*i*(i//4+1), 3)

                        # # centric
                        # nfil = int(768)
                        # nblob = int(19200)
                        # ar = round(12.65, 2)
                        # spring_factor = round(0.005 + 0.008*i*(i//4+1), 3)

                        # # # ishikawa
                        # nfil = int(160)
                        # nblob = int(10242)
                        # ar = round(6.00, 2)
                        # spring_factor = round(0.005 + 0.00*i*(i//4+1), 3)
                        # nseg = 40


                        if(self.exe_name == 'cilia_ref'):
                            nfil = 1
                            nblob = 0
                        
                        self.pars_list["nswim"].append(1)
                        self.pars_list["nseg"].append(nseg)
                        self.pars_list["nfil"].append(nfil)
                        self.pars_list["nblob"].append(nblob)
                        self.pars_list["ar"].append(ar)
                        self.pars_list["spring_factor"].append(spring_factor)
                        self.pars_list["force_mag"].append(force_mag)
                        self.pars_list["seg_sep"].append(seg_sep)
        # Write rules to sim list file
        self.write_rules()

    def delete_files(self):
        util.delete_files_in_directory(self.dir)

    def view_files(self):
        util.view_files_in_directory(self.dir)

    def write_rules(self):
        sim = configparser.ConfigParser()
        sim.add_section('Parameter list')
        for key, value in self.pars_list.items():
            sim['Parameter list'][key] = ', '.join(map(str, value))
        with open(self.dir+"rules.ini", 'w') as configfile:
            sim.write(configfile, space_around_delimiters=False)

    def read_rules(self):
        sim = configparser.ConfigParser()
        try:
            sim.read(self.dir+"rules.ini")
            for key, value in self.pars_list.items():
                if(key in sim["Parameter list"]):
                    self.pars_list[key] = [float(x) for x in sim["Parameter list"][key].split(', ')][0::1]
            self.num_sim = len(self.pars_list["nfil"])            
        except:
            print("WARNING: " + self.dir + "rules.ini not found.")

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
            self.simName = f"ciliate_{self.pars_list['nfil'][i]:.0f}fil_{self.pars_list['nblob'][i]:.0f}blob_{self.pars_list['ar'][i]:.2f}R_{self.pars_list['spring_factor'][i]:.3f}torsion"
            self.write_ini("Filenames", "simulation_file", self.simName)
            self.write_ini("Filenames", "simulation_dir", self.dir)
            self.write_ini("Filenames", "simulation_readphase_name", f"phases{int(i)}.dat")
            self.write_ini("Filenames", "simulation_readangle_name", f"angles{int(i)}.dat")


            command = f"export OPENBLAS_NUM_THREADS=1; \
                        export CUDA_VISIBLE_DEVICES={self.cuda_device}; \
                        ./bin/{self.exe_name} > terminal_outputs/output_{self.date}_{self.pars_list['nfil'][i]:.0f}fil_{i}.out"

            os.system(command)
            try:
                util.copy_last_line(self.dir + self.simName + '_filament_phases.dat', self.dir + f"phases{int(i)}.dat")
            except:
                print('**********phase file non-exist*********')