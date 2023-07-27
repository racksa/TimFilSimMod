import configparser
import os

class DRIVER:

    def __init__(self, fileName):
        self.globals_name = 'globals.ini'
        self.nfil = 16
        self.nblob = 3000
        self.ar = 5
        self.spring_factor = 2

    def write_pars(self, section, variable, value):
        config = configparser.ConfigParser()
        config.read('globals.ini')

        config.set(section, variable, f'{value}')

        # Save the changes back to the file
        with open(self.globals_name, 'w') as configfile:
            config.write(configfile, space_around_delimiters=False)

    def run(self):
        self.cuda_device = 5

        for i in range(9, 10):
            self.nfil = int(64 * (1+0.2*i)**2)
            self.nblob = int(1000 * (1+0.2*i)**2)
            self.ar = 3 * (1+0.2*i)
            self.spring_factor = 2.0

            self.write_pars("Parameters", "nfil", int(self.nfil))
            self.write_pars("Parameters", "nblob", int(self.nblob))
            self.write_pars("Parameters", "ar", int(self.ar))
            self.write_pars("Parameters", "spring_factor", int(self.spring_factor))
            self.write_pars("Filenames", "simulation_dir", "data/expr_sims/global/")
            self.write_pars("Filenames", "simulation_file", f"ciliate_{self.nfil}fil_{self.nblob}blob_{self.ar:.2f}R_{self.spring_factor:.2f}torsion")
            
            command = f"export OPENBLAS_NUM_THREADS=1; \
                        export CUDA_VISIBLE_DEVICES={self.cuda_device}; \
                        ./bin/cilia"
            
            os.system(command)