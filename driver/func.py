import multiprocessing
import configparser
import os

class DRIVER:

    def __init__(self, fileName):
        self.fileName = fileName

    def read_pars(self):
        # Read the parameters from the INI file
        config = configparser.ConfigParser()
        config.read(self.fileName)

        # Access the parameters
        nfil = config.getint('Parameters', 'nfil')
        nblob = config.getint('Parameters', 'nblob')

        print("Current values:")
        print("nfil:", nfil)
        print("nblob:", nblob)

    def write_pars(self, nfil, nblob):
        # Modify the parameters
        config = configparser.ConfigParser()
        config.read(self.fileName)

        config.set('Parameters', 'nfil', f'{nfil}')
        config.set('Parameters', 'nblob', f'{nblob}')

        # Save the changes back to the file
        with open(self.fileName, 'w') as configfile:
            config.write(configfile)

    def modi_makefile(self, variable_name, new_value):
        # Read the content of the Makefile
        with open('makefile', 'r') as file:
            content = file.read().splitlines()
        
        for i, line in enumerate(content):
            if line.startswith(f"{variable_name} = "):
                content[i] = f"{variable_name} = {new_value}"
        
        # Write the updated content back to the Makefile
        updated_content = '\n'.join(content)        
        with open('makefile', 'w') as file:
            file.write(updated_content)
        
    def modify_and_run(self, nfil, nblob):
        self.modi_makefile('CNFIL', f'{nfil}')
        self.modi_makefile('CNBLOB', f'{nblob}')
        os.system("make cilia_auto &")
        print(nfil, nblob)

    def compile_all(self):
        nfils = [16, 24]
        nblobs = [4000, 5000]
        processes = [multiprocessing.Process(target=self.modify_and_run, args=(nfils[i], nblobs[i])) for i in range(len(nfils))]
        for process in processes:
            process.start()

        for process in processes:
            process.join()

        # for i in range(2):
        #     nfil = 16 + 8*i
        #     nblob = 4000 + 2000*i
        #     # self.write_pars(nfil, nblob)
        #     # self.read_pars()
        #     self.modi_makefile('CNFIL', f'{nfil}')
        #     self.modi_makefile('CNBLOB', f'{nblob}')
            
        #     multiprocessing.Process(target=self.run_command, args=()) 

    def run_all(self):
        # Command to be executed in the terminal
        command = f"export OPENBLAS_NUM_THREADS=1; \
                    export CUDA_VISIBLE_DEVICES=3; \
                    ./bin/cilia"