import configparser
import os
import math

class DRIVER:

    def __init__(self):
        self.globals_name = '../filament/globals.ini' # name and directory!
        self.dir = "data/"
        self.pars_list = {
                     "nseg": [],
                     "nfil": [],
                     "follower_force_nondim": [],
                     "spring_const": [],
                     "seeding_y_spacing_const": [],
                     "seeding_x_spacing_const": [],
                     "ndts": [],
                     "dt": []}
        
        # self.sweep_shape = (3, 8, 6, 1)
    
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

    def change_variables(self, nseg, nfil, f, k, y_spacing_const, x_spacing_const, ndts, dt):
        # nfil = int( 400+12*j)
        
        self.pars_list["nseg"].append(nseg)
        self.pars_list["nfil"].append(nfil)
        self.pars_list["follower_force_nondim"].append(f)
        self.pars_list["spring_const"].append(k)
        self.pars_list["seeding_y_spacing_const"].append(y_spacing_const)
        self.pars_list["seeding_x_spacing_const"].append(x_spacing_const)
        self.pars_list["ndts"].append(ndts)
        self.pars_list["dt"].append(dt)


    def update_globals_file(self):
        self.create_ini()

        # Iterate through the sim list and write to .ini file and execute
        for key, value in self.pars_list.items():
            self.write_ini("Parameters", key, float(self.pars_list[key][0]))
        simulation_file = f"FFLength_{self.pars_list['nfil'][0]:.0f}fil_{self.pars_list['nseg'][0]:.0f}seg_{self.pars_list['follower_force_nondim'][0]:.0f}f_{self.pars_list['spring_const'][0]:.0f}k"
        self.write_ini("Filenames", "simulation_file", simulation_file)
        self.write_ini("Filenames", "simulation_dir", self.dir)
        return simulation_file

    def run(self):
        self.update_globals_file()

        os.chdir('../filament')
        command = f"nohup ./run_cilia_sim.sh > nohup_FFLength_{self.pars_list['nfil'][0]:.0f}fil_{self.pars_list['nseg'][0]:.0f}seg_{self.pars_list['follower_force_nondim'][0]:.0f}f_{self.pars_list['spring_const'][0]:.0f}k.out &"
        
        os.system(command)