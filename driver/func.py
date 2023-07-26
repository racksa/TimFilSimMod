import multiprocessing
import configparser
import os

class DRIVER:

    def __init__(self, fileName):
        self.globals_name = 'globals.ini'
        self.nfil = 16
        self.nblob = 3000
        self.ar = 5
        self.spring_factor = 2

    def read_pars(self):
        # Read the parameters from the INI file
        config = configparser.ConfigParser()
        config.read(self.globals_name)

        # Access the parameters
        nfil = config.getint('Parameters', 'nfil')
        nblob = config.getint('Parameters', 'nblob')

    def write_pars(self, section, variable, value):
        # Modify the parameters
        config = configparser.ConfigParser()
        config.read('globals.ini')

        config.set(section, variable, f'{value}')

        # Save the changes back to the file
        with open(self.globals_name, 'w') as configfile:
            config.write(configfile, space_around_delimiters=False)

    def run(self):
        self.cuda_device = 2

        for i in range(5):
            self.nfil = int(64 * (1+0.2*i)**2)
            self.nblob = int(1000 * (1+0.2*i)**2)
            self.ar = 3 * (1+0.2*i)

            self.write_pars("Parameters", "nfil", int(self.nfil))
            self.write_pars("Parameters", "nblob", int(self.nblob))
            self.write_pars("Parameters", "ar", int(self.ar))
            self.write_pars("Parameters", "spring_factor", int(self.spring_factor))
            self.write_pars("Filenames", "simulation_dir", "data/expr_sims/fixed_number_density/")
            self.write_pars("Filenames", "simulation_file", f"ciliate_{self.nfil}fil_{self.nblob}blob_{self.ar:.2f}R_{self.spring_factor:.2f}torsion")
            
            command = f"export OPENBLAS_NUM_THREADS=1; \
                        export CUDA_VISIBLE_DEVICES={self.cuda_device}; \
                        ./bin/cilia"
            
            os.system(command)

    # def modi_makefile(self, variable_name, new_value):
    #     # Read the content of the Makefile
    #     with open('makefile', 'r') as file:
    #         content = file.read().splitlines()
        
    #     for i, line in enumerate(content):
    #         if line.startswith(f"{variable_name} = "):
    #             content[i] = f"{variable_name} = {new_value}"
        
    #     # Write the updated content back to the Makefile
    #     updated_content = '\n'.join(content)        
    #     with open('makefile', 'w') as file:
    #         file.write(updated_content)
        
    # def modify_and_run(self, nfil, nblob):
    #     self.modi_makefile('CNFIL', f'{nfil}')
    #     self.modi_makefile('CNBLOB', f'{nblob}')
    #     os.system("make cilia_auto &")
    #     print(nfil, nblob)

    # def compile_all(self):
    #     nfils = [16, 24]
    #     nblobs = [4000, 5000]
    #     processes = [multiprocessing.Process(target=self.modify_and_run, args=(nfils[i], nblobs[i])) for i in range(len(nfils))]
    #     for process in processes:
    #         process.start()

    #     for process in processes:
    #         process.join()

    