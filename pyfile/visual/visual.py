import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import myIo
import util
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import matplotlib.patches as patches
import matplotlib.animation as animation
import subprocess
import random
import time

overlap = {name for name in mcolors.CSS4_COLORS
           if f'xkcd:{name}' in mcolors.XKCD_COLORS}
color_list = list(overlap)
color_list.sort()
color_list.pop(2)
color_list.pop(2)

# color_list = ['#00ffff', '#faebD7', '#838bbb', '#0000ff', '	#8a2be2', '#ff4040', '#7fff00', '#ff6103', '#9932cc', '#ff1493', '#030303']

enable_periodic = False
big_sphere = False

simDir = 'data/box_height_sim/100fil_sims/'
simName = simDir + 'test_fil_1000_1000_1375'


selection = 9
for i in range(selection, selection+1):
    nfil = int(64 * (1+0.2*i)**2)
    nblob = int(1000 * (1+0.2*i)**2)
    ar = 3 * (1+0.2*i)
    spring_factor = 2
    print(f" nfil={nfil}\n nblob={nblob}\n ar={ar:.2f}")
    simDir = 'data/expr_sims/fixed_number_density/'
    simName = simDir + f'ciliate_{nfil}fil_{nblob}blob_{ar:.2f}R_{spring_factor:.2f}torsion'

# simDir = 'data/phase_model/fixed_density/'
# # simName = simDir + 'test_bab_64fil_2000blob_2R_2torsion'
# # simName = simDir + 'test_bab_256fil_8000blob_4R_2torsion'
# simName = simDir + 'test_bab_1024fil_32000blob_8R_2torsion'

# simDir = 'data/phase_model/'
# simName = simDir + 'test_bab_1024fil_40000blob_6R_2torsion'

# # simDir = 'data/phase_model/fixed_filament/'
# # simName = simDir + 'test_bab_128fil_6000blob_4R_2torsion'
# # simName = simDir + 'test_bab_128fil_12000blob_6R_2torsion'

# simDir = 'data/phase_model/single_fil/'
# simName = simDir + 'test_bab_1fil'
# simName = simDir + 'test_bab_1fil_white'

# simDir = 'data/4096fil_sims/'
# simName = simDir + 'test_fil_6400_6400_800'

# simDir = 'data/box_height_sim/1fil_sims/'
# simName = simDir + 'test_fil_100_100_1375'
# simDir = 'data/box_height_sim/256fil_sims/'
# simName = simDir + 'test_fil_1600_1600_1600'
# simDir = 'data/box_height_sim/1024fil_sims/'
# simName = simDir + 'test_fil_3200_3200_1000'

# simDir = 'data/rod_sims/2304rod_sims/'
# simName = simDir + 'test_rod_960_960_60'
# simName = simDir + 'test_rod_1920_1920_60'
# simName = simDir + 'test_rod_3840_3840_60'

# simDir = 'data/rod_sims/8100rod_sims/'
# simName = simDir + 'test_rod_1800_1800_75'
# simName = simDir + 'test_rod_3600_3600_75'
# simName = simDir + 'test_rod_7200_7200_75'

# simDir = 'data/rod_sims/rod7744/'
# simName = simDir + 'test_rod_7744'

patternDatafileName = simName + '_pattern.dat'
superpuntoDatafileName = simName + '_superpunto.dat'
fcmPosfileName = simName + '_flow_pos.dat'
fcmForcefileName = simName + '_flow_force.dat'
fcmTorquefileName = simName + '_flow_torque.dat'

Lx, Ly, Lz = myIo.get_boxsize_from_name(simName)

if(np.isinf(np.array([Lx, Ly, Lz])).any() and enable_periodic):
    Lx = 1000.
    Ly = 1000.
    Lz = 1000.
    print(f"Manually setting the boxsize to ({Lx}, {Ly}, {Lz}).")

# Lz *= 5

# fil_4096
# Lx = 3840.
# Ly = 3840.
# Lz = 240.
# 256*256*16

# rod_1024_2
# Lx = 512.
# Ly = 512.
# Lz = 16.*10
# 384*384*12

#rod_1024_l
# Lx = 1280.
# Ly = 1280.
# Lz = 40.*10
# 256*256*8

# rod_16384
# Lx = 2048
# Ly = 2048
# Lz = 32*2

# rod_7744
# Lx = 1760.
# Ly = 1760.
# Lz = 10.3125*10
# 1024*1024*6


class VISUAL:

    def __init__(self):
        self.pars = myIo.read_pars(simName + '.par')
        if(not 'PRESCRIBED_CILIA' in self.pars):
            self.pars['PRESCRIBED_CILIA'] = 0
        if(self.pars['NBLOB']>0):
            self.blob_references = myIo.read_blob_references(simName + '_blob_references.dat')
        if(self.pars['NFIL']>0):
            self.fil_references = myIo.read_fil_references(simName + '_fil_references.dat')
        self.dt = self.pars['DT']*self.pars['PLOT_FREQUENCY_IN_STEPS']
        self.L = 14.14*self.pars['NBLOB']/22.
        self.frames = min(300001, sum(1 for line in open(simName + '_body_states.dat')))

        self.plot_end_frame = self.frames
        self.plot_start_frame = max(0, self.plot_end_frame-30)
        self.plot_interval = 1

        self.plot_hist_frame = np.array([self.frames-1])
        self.plot_seg_frames = [self.plot_end_frame-1-2*i for i in range(14)]
        self.plot_rod_frame = np.array([self.plot_end_frame-1])
        self.plot_multi_rod_frames = np.array([100, 1000, 10000])
        self.plot_single_fil_frames = [self.plot_end_frame-1-2*i for i in range(15)]
        # self.plot_single_fil_frames = [self.plot_end_frame-1]
        self.plot_phase_frames = [self.plot_end_frame-1]
        self.fcm_frame = self.plot_end_frame-1

        self.phase_video = False
        self.ciliate_video = False

        self.output_to_superpunto = False
        self.output_to_fcm = False

    def set_plot_dim(self, dim):
        self.plot_dim = dim

    def enable_superpunto(self):
        self.output_to_superpunto = True
        myIo.clean_file(superpuntoDatafileName)
        myIo.clean_file(patternDatafileName)

    def enable_fcm(self):
        self.output_to_fcm = True
        myIo.clean_file(fcmPosfileName)
        myIo.clean_file(fcmForcefileName)
        myIo.clean_file(fcmTorquefileName)

    def write_data(self, x, r, filename, box=True, center=True, superpunto=True, color=0):
        plot_x = x.copy()
        if(box):
            plot_x[0] = util.box(plot_x[0], Lx)
            plot_x[1] = util.box(plot_x[1], Ly)
            plot_x[2] = util.box(plot_x[2], Lz)
            if(center):
                plot_x[0] -= 0.5*Lx
                plot_x[1] -= 0.5*Ly
        if superpunto:
            plot_x[1] = -plot_x[1]
            plot_x[2] = -plot_x[2]
            myIo.write_line(str(plot_x[0]) + ' ' +\
                            str(plot_x[1]) + ' ' +\
                            str(plot_x[2]) + ' ' +\
                            str(r) + ' ' +\
                            str(color),
                            filename)
        else:
            myIo.write_line(str(plot_x[0]) + ' ' +\
                            str(plot_x[1]) + ' ' +\
                            str(plot_x[2]),
                            filename)
      
    def plot(self):
        if (self.plot_dim == 2):
            ax = plt.figure().add_subplot(1,1,1)
            ax.axis('equal')
        if (self.plot_dim == 3):
            ax = plt.figure().add_subplot(projection='3d')

        start = time.time()
        seg_states_f = open(simName + '_seg_states.dat', "r")
        body_states_f = open(simName + '_body_states.dat', "r")
        if(self.pars['NBLOB']>0):
            blob_forces_f = open(simName + '_blob_forces.dat', "r")
        seg_forces_f = open(simName + '_seg_forces.dat', "r")
        print("open files time = ",(time.time()-start))
        if (self.pars['PRESCRIBED_CILIA'] == 1):
            fil_phases_f = open(simName + '_filament_phases.dat', "r")
            fil_angles_f = open(simName + '_filament_shape_rotation_angles.dat', "r")
            AR, torsion = myIo.get_ciliate_data_from_name(simName)
            radius = 0.5*AR*2.2*int(self.pars['NSEG'])

        for i in range(self.plot_end_frame):
            print(" frame ", i, "/", self.frames, "          ", end="\r")
            body_states_str = body_states_f.readline()
            if(self.pars['NFIL']>0):
                seg_states_str = seg_states_f.readline()
                if(self.output_to_fcm):
                    seg_forces_str = seg_forces_f.readline()
                if (self.pars['PRESCRIBED_CILIA'] == 1):
                    fil_phases_str = fil_phases_f.readline()
                    # fil_angles_str = fil_angles_f.readline()

                    fil_phases = np.array(fil_phases_str.split()[1:], dtype=float)
                    fil_phases = util.box(fil_phases, 2*np.pi)
            if(self.pars['NBLOB']>0 and self.output_to_fcm):
                blob_forces_str = blob_forces_f.readline()

            if(i%self.plot_interval==0 and i>=self.plot_start_frame):
                body_states = np.array(body_states_str.split()[1:], dtype=float)
                if(self.pars['NFIL']>0):
                    seg_states = np.array(seg_states_str.split()[1:], dtype=float)
                    if(self.output_to_fcm):
                        seg_forces = np.array(seg_forces_str.split()[1:], dtype=float)
                if(self.pars['NBLOB']>0 and self.output_to_fcm):
                    blob_forces = np.array(blob_forces_str.split()[1:], dtype=float)

                

                if(self.output_to_superpunto):
                    myIo.write_line('#', superpuntoDatafileName)
                if((self.output_to_fcm and i == self.fcm_frame) or (self.output_to_superpunto)):
                    for swim in range(int(self.pars['NSWIM'])):
                        body_pos = body_states[7*swim : 7*swim+3]
                        R = util.rot_mat(body_states[7*swim+3 : 7*swim+7])
                        R = np.linalg.inv(R)
                        R = np.eye(3)
                        # To find blob position
                        if(big_sphere):
                            self.write_data(body_pos, radius, superpuntoDatafileName, enable_periodic, color=16777215)
                        for blob in range(int(self.pars['NBLOB'])):
                            blob_x, blob_y, blob_z = util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[3*blob:3*blob+3])
                            if(self.output_to_superpunto):
                                if(not big_sphere):
                                    self.write_data(np.matmul(R, [blob_x, blob_y, blob_z]), float(self.pars['RBLOB']), superpuntoDatafileName, enable_periodic, color=16777215)
                            if(self.output_to_fcm and i == self.fcm_frame):
                                blob_fx, blob_fy, blob_fz = blob_forces[3*blob : 3*blob+3]
                                self.write_data([blob_x, blob_y, blob_z], 0, fcmPosfileName, box=True, center=False, superpunto=False)
                                self.write_data([blob_fx, blob_fy, blob_fz], 0, fcmForcefileName, box=False, center=False, superpunto=False)
                                self.write_data([0, 0, 0], 0, fcmTorquefileName, box=False, center=False, superpunto=False)

                        # Robot arm to find segment position (Ignored plane rotation!)
                        for fil in range(int(self.pars['NFIL'])):
                            fil_color = int("000000", base=16)
                            if (self.pars['PRESCRIBED_CILIA'] == 0):
                                fil_i = int(4*fil*self.pars['NSEG'])
                                fil_base_x, fil_base_y, fil_base_z = body_pos + np.matmul(R, self.fil_references[3*fil : 3*fil+3])
                                old_seg_pos = np.array([fil_base_x, fil_base_y, fil_base_z])
                            elif (self.pars['PRESCRIBED_CILIA'] == 1):
                                fil_i = int(3*fil*self.pars['NSEG'])
                                old_seg_pos = seg_states[fil_i : fil_i+3]
                                
                                # WRITE A FUNCTION FOR THIS!!
                                cmap_name = 'hsv'
                                cmap = plt.get_cmap(cmap_name)
                                rgb_color = cmap(fil_phases[fil]/(2*np.pi))[:3]  # Get the RGB color tuple
                                rgb_hex = mcolors.rgb2hex(rgb_color)[1:]  # Convert RGB to BGR hexadecimal format
                                bgr_hex = rgb_hex[4:]+rgb_hex[2:4]+rgb_hex[:2]
                                fil_color = int(bgr_hex, base=16)
                                # print("\n", bgr_hex, fil_color, "\t")

                            if(self.output_to_superpunto):
                                self.write_data(np.matmul(R, old_seg_pos-[50,50,50]), float(self.pars['RSEG']), superpuntoDatafileName, enable_periodic, True, True, color=fil_color)
                            if(self.output_to_fcm):
                                seg_fx, seg_fy, seg_fz, seg_tx, seg_ty, seg_tz = seg_forces[0 : 6]
                                self.write_data(old_seg_pos, 0, fcmPosfileName, box=True, center=False, superpunto=False)
                                self.write_data([seg_fx, seg_fy, seg_fz], 0, fcmForcefileName, box=False, center=False, superpunto=False)
                                self.write_data([seg_tx, seg_ty, seg_tz], 0, fcmTorquefileName, box=False, center=False, superpunto=False)

                            for seg in range(1, int(self.pars['NSEG'])):
                                if (self.pars['PRESCRIBED_CILIA'] == 0):
                                    q1 = seg_states[fil_i+4*(seg-1) : fil_i+4*seg]
                                    q2 = seg_states[fil_i+4*seg : fil_i+4*seg+4]
                                    
                                    t1 = util.find_t(q1)
                                    t2 = util.find_t(q2)
                                    
                                    seg_pos = old_seg_pos + 0.5*self.pars['DL']*(t1 + t2)
                                    old_seg_pos = seg_pos
                                elif (self.pars['PRESCRIBED_CILIA'] == 1):
                                    seg_pos = seg_states[fil_i+3*(seg-1) : fil_i+3*seg] 
                                if(self.output_to_superpunto):
                                    self.write_data(np.matmul(R, seg_pos-[50,50,50]), float(self.pars['RSEG']), superpuntoDatafileName, enable_periodic, True, True, color=fil_color)
                                if(self.output_to_fcm):
                                    seg_fx, seg_fy, seg_fz, seg_tx, seg_ty, seg_tz = seg_forces[6*seg : 6*seg+6]
                                    self.write_data(seg_pos, 0, fcmPosfileName, box=True, center=False, superpunto=False)
                                    self.write_data([seg_fx, seg_fy, seg_fz], 0, fcmForcefileName, box=False, center=False, superpunto=False)
                                    self.write_data([seg_tx, seg_ty, seg_tz], 0, fcmTorquefileName, box=False, center=False, superpunto=False)

        if(self.output_to_fcm):
            fcm_directory = "../CUFCM/data/flow_data/"
            subprocess.call("cp " + fcmPosfileName + " " + fcm_directory, shell=True)
            subprocess.call("cp " + fcmForcefileName + " " + fcm_directory, shell=True)
            subprocess.call("cp " + fcmTorquefileName + " " + fcm_directory, shell=True)

    # Filaments
    def plot_seg_vel(self):
        seg_states_f = open(simName + '_seg_states.dat', "r")
        seg_vels_f = open(simName + '_seg_vels.dat', "r")
        body_states_f = open(simName + '_body_states.dat', "r")

        ax = plt.figure().add_subplot(1,1,1)
        end_pos = np.zeros((int(self.pars['NFIL']), 3))
        end_vel = np.zeros((int(self.pars['NFIL']), 3))
        tip_traj = np.zeros((len(self.plot_seg_frames), int(self.pars['NFIL']), 3))
        current_time = 0

        for i in range(self.plot_end_frame):
            print(" frame ", i, "/", self.frames, "          ", end="\r")
            body_states_str = body_states_f.readline()
            seg_states_str = seg_states_f.readline()
            seg_vels_str = seg_vels_f.readline()

            if(i in self.plot_seg_frames):
                seg_states = np.array(seg_states_str.split()[1:], dtype=float)
                seg_vels = np.array(seg_vels_str.split()[1:], dtype=float)
                body_states = np.array(body_states_str.split()[1:], dtype=float)

                R = util.rot_mat(body_states[3 : 7])

                # Robot arm to find segment position (Ignored plane rotation!)
                for fil in range(int(self.pars['NFIL'])):
                    if (self.pars['PRESCRIBED_CILIA'] == 0):
                        fil_i = int(4*fil*self.pars['NSEG'])
                        fil_base_x, fil_base_y, fil_base_z = np.matmul(R, self.fil_references[3*fil : 3*fil+3])
                        old_seg_pos = np.array([fil_base_x, fil_base_y, fil_base_z])
                    elif (self.pars['PRESCRIBED_CILIA'] == 1):
                        fil_i = int(3*fil*self.pars['NSEG'])
                        old_seg_pos = seg_states[fil_i : fil_i+3]

                    for seg in range(1, int(self.pars['NSEG'])):
                        if (self.pars['PRESCRIBED_CILIA'] == 0):
                            q1 = seg_states[fil_i+4*(seg-1) : fil_i+4*seg]
                            q2 = seg_states[fil_i+4*seg : fil_i+4*seg+4]
                            
                            t1 = util.find_t(q1)
                            t2 = util.find_t(q2)
                            
                            seg_pos = old_seg_pos + 0.5*self.pars['DL']*(t1 + t2)
                            old_seg_pos = seg_pos
                        elif (self.pars['PRESCRIBED_CILIA'] == 1):
                            seg_pos = seg_states[fil_i+3*(seg-1) : fil_i+3*seg]

                    end_pos[fil] = seg_pos
                    end_vel[fil] = seg_vels[6*int(self.pars['NSEG'])*fil + 6*int(self.pars['NSEG']) - 6 : 6*int(self.pars['NSEG'])*fil + 6*int(self.pars['NSEG']) - 3]
                    tip_traj[current_time][fil] = seg_pos
                current_time += 1
    
        for fil in range(int(self.pars['NFIL'])):
            ax.plot(tip_traj[:,fil,0], tip_traj[:,fil,1], c='black')
            ax.arrow(end_pos[fil,0], end_pos[fil,1], end_vel[fil,0], end_vel[fil,1], color='r', lw=1, head_width=18)

        ax.set_ylabel(r"y")
        ax.set_xlabel(r"x")
        ax.set_xlim(0, Lx)
        ax.set_ylim(0, Ly)

        L = 2.2*self.pars['NSEG']
        x_ticks = np.linspace(0, Lx, 6)
        x_ticks_label = ["{:.2f}L".format(x/L) for x in x_ticks]
        ax.set_xticks(x_ticks, x_ticks_label)
        y_ticks = np.linspace(0, Ly, 6)
        y_ticks_label = ["{:.2f}L".format(y/L) for y in y_ticks]
        ax.set_yticks(y_ticks, y_ticks_label)

        ax.set_aspect('equal')
        nfil = int(self.pars['NFIL'])
        # ax.legend(title=rf'{int(Lx)}$\times${int(Ly)}$\times${int(Lz)}', loc='upper left')
        plt.savefig(f'fig/fil_vel_{nfil}fil_{int(Lx)}_{int(Ly)}_{int(Lz)}.eps', bbox_inches = 'tight', format='eps')
        plt.savefig(f'fig/fil_vel_{nfil}fil_{int(Lx)}_{int(Ly)}_{int(Lz)}.png', bbox_inches = 'tight', format='png')
        plt.savefig(f'fig/fil_vel_{nfil}fil_{int(Lx)}_{int(Ly)}_{int(Lz)}.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

    def plot_pattern(self):

        seg_states_f = open(simName + '_seg_states.dat', "r")
        body_states_f = open(simName + '_body_states.dat', "r")
        for i in range(self.plot_end_frame):
            print(" frame ", i, "/", self.frames, "          ", end="\r")
            body_states_str = body_states_f.readline()
            seg_states_str = seg_states_f.readline()

            myIo.write_line('#', superpuntoDatafileName)

            if(i%self.plot_interval==0 and i>=self.plot_start_frame):
                body_states = np.array(body_states_str.split()[1:], dtype=float)
                seg_states = np.array(seg_states_str.split()[1:], dtype=float)
                for swim in range(int(self.pars['NSWIM'])):
                    body_pos = body_states[7*swim : 7*swim+3]
                    R = util.rot_mat(body_states[7*swim+3 : 7*swim+7])
                    # Robot arm to find segment position (Ignored plane rotation!)
                    for fil in range(int(self.pars['NFIL'])):
                        
                        fil_i = int(4*fil*self.pars['NSEG'])
                        fil_base_x, fil_base_y, fil_base_z = body_pos + np.matmul(R, self.fil_references[3*fil : 3*fil+3])
                        old_seg_pos = np.array([fil_base_x, fil_base_y, fil_base_z])

                        self.write_data(old_seg_pos, float(self.pars['RSEG']), patternDatafileName, enable_periodic, True, True)

                        for seg in range(1, int(self.pars['NSEG'])):
                            q1 = seg_states[fil_i+4*(seg-1) : fil_i+4*seg]
                            q2 = seg_states[fil_i+4*seg : fil_i+4*seg+4]
                            
                            t1 = util.find_t(q1)
                            t2 = util.find_t(q2)
                            
                            seg_pos = old_seg_pos + 0.5*self.pars['DL']*(t1 + t2)
                            old_seg_pos = seg_pos
                            
                            self.write_data(seg_pos, float(self.pars['RSEG']), patternDatafileName, enable_periodic, True, True)
    
    def plot_fil3d(self):
        fig = plt.figure(figsize=(16,4),dpi=400)
        ax = fig.add_subplot(projection='3d')
        seg_states_f = open(simName + '_seg_states.dat', "r")
        body_states_f = open(simName + '_body_states.dat', "r")
        current_frame=0
        for i in range(self.plot_end_frame):
            print(" frame ", i, "/", self.frames, "          ", end="\r")
            body_states_str = body_states_f.readline()
            seg_states_str = seg_states_f.readline()

            if(i in self.plot_single_fil_frames):
                body_states = np.array(body_states_str.split()[1:], dtype=float)
                seg_states = np.array(seg_states_str.split()[1:], dtype=float)
                for swim in range(int(self.pars['NSWIM'])):
                    body_pos = body_states[7*swim : 7*swim+3]
                    R = util.rot_mat(body_states[7*swim+3 : 7*swim+7])
                    # Robot arm to find segment position (Ignored plane rotation!)
                    for fil in range(int(self.pars['NFIL'])):                    
                        segs_pos = np.zeros((int(self.pars['NSEG']), 3))
                        if (self.pars['PRESCRIBED_CILIA'] == 0):
                            fil_i = int(4*fil*self.pars['NSEG'])
                            fil_base_x, fil_base_y, fil_base_z = np.matmul(R, self.fil_references[3*fil : 3*fil+3])
                            fil_base_z -= 50
                            old_seg_pos = np.array([fil_base_x, fil_base_y, fil_base_z])
                        elif (self.pars['PRESCRIBED_CILIA'] == 1):
                            fil_i = int(3*fil*self.pars['NSEG'])
                            old_seg_pos = seg_states[fil_i : fil_i+3]
                        segs_pos[0] = old_seg_pos

                        for seg in range(1, int(self.pars['NSEG'])):
                            if (self.pars['PRESCRIBED_CILIA'] == 0):
                                q1 = seg_states[fil_i+4*(seg-1) : fil_i+4*seg]
                                q2 = seg_states[fil_i+4*seg : fil_i+4*seg+4]
                                
                                t1 = util.find_t(q1)
                                t2 = util.find_t(q2)
                                
                                seg_pos = old_seg_pos + 0.5*self.pars['DL']*(t1 + t2)
                                old_seg_pos = seg_pos
                            elif (self.pars['PRESCRIBED_CILIA'] == 1):
                                seg_pos = seg_states[fil_i+3*(seg-1) : fil_i+3*seg]
                            segs_pos[seg] = seg_pos
                        
                        ax.plot(segs_pos[:,0], segs_pos[:,1], segs_pos[:,2], c = 'black', alpha=0.1+current_frame*0.06, lw=0.3)
                current_frame += 1
        
        ax.set_proj_type('ortho')
        ax.set_proj_type('persp', 0.08)  # FOV = 157.4 deg
        ax.view_init(elev=15, azim=90)
        ax.dist=10
        ax.axis('off')
        ax.grid(False)

        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        

        L = 2.2*self.pars['NSEG']
        x_ticks = np.linspace(0, Lx, 6)
        x_ticks_label = ["{:.2f}L".format(x/L) for x in x_ticks]
        ax.set_xticks(x_ticks, x_ticks_label)
        y_ticks = np.linspace(0, Ly, 6)
        y_ticks_label = ["{:.2f}L".format(y/L) for y in y_ticks]
        ax.set_yticks(y_ticks, y_ticks_label)
        # z_ticks = np.arange(0, Lz, step=50)
        # z_ticks_label = ["{:.2f}L".format(z/L) for z in z_ticks]
        # ax.set_zticks(z_ticks, z_ticks_label)

        ax.text(100, 50, 80, "H={:.2f}L".format(Lz/L))
        ax.set_zticks([], [])

        ax.set_xlim(0, Lx)
        ax.set_ylim(0, Ly)
        ax.set_zlim(0, 50)
        ax.set_box_aspect(aspect = (1,1,50/Lx))        
        nfil = int(self.pars['NFIL'])
        # fig.savefig(f'fig/fil_{nfil}_fil_{int(Lx)}_{int(Ly)}_{int(Lz)}.eps', bbox_inches = 'tight',  format='eps')
        fig.savefig(f'fig/fil_{nfil}_fil_{int(Lx)}_{int(Ly)}_{int(Lz)}.png', bbox_inches = 'tight',  format='png')
        fig.savefig(f'fig/fil_{nfil}_fil_{int(Lx)}_{int(Ly)}_{int(Lz)}.pdf', bbox_inches = 'tight', format='pdf')
        plt.show()

    def plot_seg_force(self):
        seg_forces_f = open(simName + '_seg_forces.dat', "r")
        tether_force_f = open(simName + '_tether_force.dat', "r")

        NF0 = float(self.pars['DIMENSIONLESS_FORCE']) * float(self.pars['NFIL'])
        T0 = float(self.pars['STEPS_PER_PERIOD']) / float(self.pars['PLOT_FREQUENCY_IN_STEPS'])

        compute_start = self.plot_start_frame
        compute_start = 0

        # total_base_force = np.zeros((self.plot_end_frame - compute_start))
        total_tether_force = np.zeros((self.plot_end_frame - compute_start))

        ax = plt.figure().add_subplot(1,1,1)
        for i in range(self.plot_end_frame):
            print(" frame ", i, "/", self.frames, "          ", end="\r")
            # seg_forces_str = seg_forces_f.readline()
            tether_forces_str = tether_force_f.readline()

            if(i >= compute_start):
                # seg_forces = np.array(seg_forces_str.split()[1:], dtype=float)
                tether_force = np.array(tether_forces_str.split()[1:], dtype=float)

                for fil in range(int(self.pars['NFIL'])):

                    # total_base_force[i-compute_start] += np.linalg.norm(seg_forces[6*fil : 6*fil+3])
                    total_tether_force[i-compute_start] += np.linalg.norm(tether_force[3*fil : 3*fil+3])

        time_array = np.arange(compute_start, self.plot_end_frame )/T0
        ax.set_ylabel(r"$<\lambda>/F_0$")
        ax.set_xlabel(r"$t/T_0$")
        # ax.plot(time_array[1:], total_base_force[1:])
        ax.plot(time_array[1:], total_tether_force[1:]/NF0)
        ax.set_xlim(left=0)
        # ax.set_ylim(min(total_tether_force[min(self.frames, 100):])/NF0*0.8, \
        #             max(total_tether_force[min(self.frames, 100):])/NF0*1.1)
        nfil = int(self.pars['NFIL'])
        plt.savefig(f'fig/tether_force_{nfil}fil_{int(Lx)}_{int(Ly)}_{int(Lz)}.eps', format='eps')
        plt.savefig(f'fig/tether_force_{nfil}fil_{int(Lx)}_{int(Ly)}_{int(Lz)}.png', format='png')
        plt.show()

    # Phase model
    def plot_phase_heatmap(self):
        colormap = 'cividis'
        colormap = 'twilight_shifted'
        fil_phases_f = open(simName + '_filament_phases.dat', "r")
        fil_angles_f = open(simName + '_filament_shape_rotation_angles.dat', "r")

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        fil_references_sphpolar = np.zeros((int(self.pars['NFIL']),3))

        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable
        norm = Normalize(vmin=0, vmax=2*np.pi)
        sm = ScalarMappable(cmap=colormap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm)
        cbar.ax.set_yticks(np.linspace(0, 2*np.pi, 7), ['0', 'π/3', '2π/3', 'π', '4π/3', '5π/3', '2π'])
        cbar.set_label(r"phase")

        ax.set_ylabel(r"$\theta$")
        ax.set_xlabel(r"$\phi$")
        ax.set_xlim(-np.pi, np.pi)
        ax.set_ylim(0, np.pi)
        ax.set_xticks(np.linspace(-np.pi, np.pi, 5), ['-π', '-π/2', '0', 'π/2', 'π'])
        ax.set_yticks(np.linspace(0, np.pi, 5), ['0', 'π/4', 'π/2', '3π/4', 'π'])
        nfil = int(self.pars['NFIL'])

        def animation_func(t):
            print(t)
            ax.cla()
            fil_phases_str = fil_phases_f.readline()
            # fil_angles_str = fil_angles_f.readline()

            fil_phases = np.array(fil_phases_str.split()[1:], dtype=float)
            fil_phases = util.box(fil_phases, 2*np.pi)
            for i in range(int(self.pars['NFIL'])):
                fil_references_sphpolar[i] = util.cartesian_to_spherical(self.fil_references[3*i: 3*i+3])
                
            ax.scatter(fil_references_sphpolar[:,1], fil_references_sphpolar[:,2], c=fil_phases, cmap=colormap)

        if(self.phase_video):
            plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
            print("video")
            ani = animation.FuncAnimation(fig, animation_func, frames=500, interval=10, repeat=False)
            plt.show()
            # FFwriter = animation.FFMpegWriter(fps=10)
            # ani.save(f'fig/fil_phase_{nfil}fil_anim.mp4', writer=FFwriter)
            
            
        else:
            for t in range(self.plot_end_frame):
                print(" frame ", t, "/", self.frames, "          ", end="\r")
                fil_phases_str = fil_phases_f.readline()
                fil_angles_str = fil_angles_f.readline()

                if(t in self.plot_phase_frames):

                    fil_phases = np.array(fil_phases_str.split()[1:], dtype=float)
                    fil_phases = util.box(fil_phases, 2*np.pi)
                    for i in range(int(self.pars['NFIL'])):
                        fil_references_sphpolar[i] = util.cartesian_to_spherical(self.fil_references[3*i: 3*i+3])
                        
                    ax.scatter(fil_references_sphpolar[:,1], fil_references_sphpolar[:,2], c=fil_phases, cmap=colormap)
                
            plt.savefig(f'fig/fil_phase_{nfil}fil.png', bbox_inches = 'tight', format='png')
            plt.savefig(f'fig/fil_phase_{nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
            plt.show()

    def plot_ciliate(self):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        start = time.time()
        seg_states_f = open(simName + '_seg_states.dat', "r")
        body_states_f = open(simName + '_body_states.dat', "r")
        print("open files time = ",(time.time()-start))

        AR, torsion = myIo.get_ciliate_data_from_name(simName)
        nfil = int(self.pars['NFIL'])
        # Set the radius of the sphere
        radius = 0.5*AR*2.2*int(self.pars['NSEG'])
        num_points = 300

        ax.set_proj_type('ortho')
        # ax.set_proj_type('persp', 0.05)  # FOV = 157.4 deg
        # ax.view_init(elev=5., azim=45)
        # ax.dist=20
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        # ax.axis('off')
        # ax.grid(False)

        def animation_func(t):
            print(t)
            ax.cla()
            ax.set_xlim(-100, 100)
            ax.set_ylim(-100, 100)
            ax.set_zlim(-100, 100)

            body_states_str = body_states_f.readline()
            seg_states_str = seg_states_f.readline()

            body_states = np.array(body_states_str.split()[1:], dtype=float)
            seg_states = np.array(seg_states_str.split()[1:], dtype=float)
            
            for swim in range(int(self.pars['NSWIM'])):
                # blob_data = np.zeros((int(self.pars['NBLOB']), 3))
                body_pos = body_states[7*swim : 7*swim+3]
                # R = util.rot_mat(body_states[7*swim+3 : 7*swim+7])
                # # To find blob position
                # for blob in range(int(self.pars['NBLOB'])):
                #     blob_x, blob_y, blob_z = util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[3*blob:3*blob+3])
                #     ax.scatter(blob_x, blob_y, blob_z)

                # Create the sphere data points
                u = np.linspace(0, 2 * np.pi, num_points)
                v = np.linspace(0, np.pi, num_points)
                x = radius * np.outer(np.cos(u), np.sin(v))
                y = radius * np.outer(np.sin(u), np.sin(v))
                z = radius * np.outer(np.ones(np.size(u)), np.cos(v))

                # Plot the sphere
                ax.plot_surface(x+body_pos[0], y+body_pos[1], z+body_pos[2], color='grey', alpha=0.5)

                # Robot arm to find segment position (Ignored plane rotation!)
                for fil in range(int(self.pars['NFIL'])):
                    fil_data = np.zeros((int(self.pars['NSEG']), 3))
                    fil_i = int(3*fil*self.pars['NSEG'])
                    fil_data[0] = seg_states[fil_i : fil_i+3]

                    for seg in range(1, int(self.pars['NSEG'])):
                        seg_pos = seg_states[fil_i+3*(seg-1) : fil_i+3*seg]
                        fil_data[seg] = seg_pos
                    ax.plot(fil_data[:,0], fil_data[:,1], fil_data[:,2], c='black', zorder = 100)

        if(self.ciliate_video):
            plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
            ani = animation.FuncAnimation(fig, animation_func, frames=500, interval=10, repeat=False)
            plt.show()
            # FFwriter = animation.FFMpegWriter(fps=10)
            # ani.save(f'fig/ciliate_{nfil}fil_anim.mp4', writer=FFwriter)
    
        else:
            for i in range(self.plot_end_frame):
                print(" frame ", i, "/", self.frames, "          ", end="\r")
                body_states_str = body_states_f.readline()
                seg_states_str = seg_states_f.readline()
                if(i==self.plot_end_frame-1):
                    # if(i%self.plot_interval==0 and i>=self.plot_start_frame):
                    body_states = np.array(body_states_str.split()[1:], dtype=float)
                    seg_states = np.array(seg_states_str.split()[1:], dtype=float)
                    
                    for swim in range(int(self.pars['NSWIM'])):
                        blob_data = np.zeros((int(self.pars['NBLOB']), 3))
                        body_pos = body_states[7*swim : 7*swim+3]
                        R = util.rot_mat(body_states[7*swim+3 : 7*swim+7])
                        # To find blob position
                        for blob in range(int(self.pars['NBLOB'])):
                            blob_x, blob_y, blob_z = util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[3*blob:3*blob+3])
                            # ax.scatter(blob_x, blob_y, blob_z)

                        # Create the sphere data points
                        u = np.linspace(0, 2 * np.pi, num_points)
                        v = np.linspace(0, np.pi, num_points)
                        x = radius * np.outer(np.cos(u), np.sin(v))
                        y = radius * np.outer(np.sin(u), np.sin(v))
                        z = radius * np.outer(np.ones(np.size(u)), np.cos(v))

                        # Plot the sphere
                        ax.plot_surface(x, y, z, color='grey', alpha=0.5)

                        # Robot arm to find segment position (Ignored plane rotation!)
                        for fil in range(int(self.pars['NFIL'])):
                            fil_data = np.zeros((int(self.pars['NSEG']), 3))
                            fil_i = int(3*fil*self.pars['NSEG'])
                            fil_data[0] = seg_states[fil_i : fil_i+3]

                            for seg in range(1, int(self.pars['NSEG'])):
                                seg_pos = seg_states[fil_i+3*(seg-1) : fil_i+3*seg]
                                fil_data[seg] = seg_pos
                            ax.plot(fil_data[:,0], fil_data[:,1], fil_data[:,2], c='black', zorder = 100)

            plt.savefig(f'fig/ciliate_{nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
            plt.show()

    def plot_multi_ciliate(self):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        start = time.time()
        seg_states_f = open(simName + '_seg_states.dat', "r")
        body_states_f = open(simName + '_body_states.dat', "r")
        print("open files time = ",(time.time()-start))

        AR, torsion = myIo.get_ciliate_data_from_name(simName)
        nfil = int(self.pars['NFIL'])
        # Set the radius of the sphere
        radius = 0.5*AR*2.2*int(self.pars['NSEG'])
        num_points = 300

        ax.set_proj_type('ortho')
        # ax.set_proj_type('persp', 0.05)  # FOV = 157.4 deg
        # ax.view_init(elev=5., azim=45)
        # ax.dist=20
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
        # ax.axis('off')
        # ax.grid(False)

        def animation_func(t):
            print(t)
            ax.cla()
            ax.set_xlim(-100, 100)
            ax.set_ylim(-100, 100)
            ax.set_zlim(-100, 100)

            body_states_str = body_states_f.readline()
            seg_states_str = seg_states_f.readline()

            body_states = np.array(body_states_str.split()[1:], dtype=float)
            seg_states = np.array(seg_states_str.split()[1:], dtype=float)
            
            for swim in range(int(self.pars['NSWIM'])):
                # blob_data = np.zeros((int(self.pars['NBLOB']), 3))
                body_pos = body_states[7*swim : 7*swim+3]
                # R = util.rot_mat(body_states[7*swim+3 : 7*swim+7])
                # # To find blob position
                # for blob in range(int(self.pars['NBLOB'])):
                #     blob_x, blob_y, blob_z = util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[3*blob:3*blob+3])
                #     ax.scatter(blob_x, blob_y, blob_z)

                # Create the sphere data points
                u = np.linspace(0, 2 * np.pi, num_points)
                v = np.linspace(0, np.pi, num_points)
                x = radius * np.outer(np.cos(u), np.sin(v))
                y = radius * np.outer(np.sin(u), np.sin(v))
                z = radius * np.outer(np.ones(np.size(u)), np.cos(v))

                # Plot the sphere
                ax.plot_surface(x+body_pos[0], y+body_pos[1], z+body_pos[2], color='grey', alpha=0.5)

                # Robot arm to find segment position (Ignored plane rotation!)
                for fil in range(int(self.pars['NFIL'])):
                    fil_data = np.zeros((int(self.pars['NSEG']), 3))
                    fil_i = int(3*fil*self.pars['NSEG'])
                    fil_data[0] = seg_states[fil_i : fil_i+3]

                    for seg in range(1, int(self.pars['NSEG'])):
                        seg_pos = seg_states[fil_i+3*(seg-1) : fil_i+3*seg]
                        fil_data[seg] = seg_pos
                    ax.plot(fil_data[:,0], fil_data[:,1], fil_data[:,2], c='black', zorder = 100)

        if(self.ciliate_video):
            plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
            ani = animation.FuncAnimation(fig, animation_func, frames=500, interval=10, repeat=False)
            plt.show()
            # FFwriter = animation.FFMpegWriter(fps=10)
            # ani.save(f'fig/ciliate_{nfil}fil_anim.mp4', writer=FFwriter)
    
        else:
            for i in range(self.plot_end_frame):
                print(" frame ", i, "/", self.frames, "          ", end="\r")
                body_states_str = body_states_f.readline()
                seg_states_str = seg_states_f.readline()
                if(i==self.plot_end_frame-1):
                    # if(i%self.plot_interval==0 and i>=self.plot_start_frame):
                    body_states = np.array(body_states_str.split()[1:], dtype=float)
                    seg_states = np.array(seg_states_str.split()[1:], dtype=float)
                    
                    for swim in range(int(self.pars['NSWIM'])):
                        blob_data = np.zeros((int(self.pars['NBLOB']), 3))
                        body_pos = body_states[7*swim : 7*swim+3]
                        R = util.rot_mat(body_states[7*swim+3 : 7*swim+7])
                        # To find blob position
                        for blob in range(int(self.pars['NBLOB'])):
                            blob_x, blob_y, blob_z = util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[3*blob:3*blob+3])
                            # ax.scatter(blob_x, blob_y, blob_z)

                        # Create the sphere data points
                        u = np.linspace(0, 2 * np.pi, num_points)
                        v = np.linspace(0, np.pi, num_points)
                        x = radius * np.outer(np.cos(u), np.sin(v))
                        y = radius * np.outer(np.sin(u), np.sin(v))
                        z = radius * np.outer(np.ones(np.size(u)), np.cos(v))

                        # Plot the sphere
                        ax.plot_surface(x, y, z, color='grey', alpha=0.5)

                        # Robot arm to find segment position (Ignored plane rotation!)
                        for fil in range(int(self.pars['NFIL'])):
                            fil_data = np.zeros((int(self.pars['NSEG']), 3))
                            fil_i = int(3*fil*self.pars['NSEG'])
                            fil_data[0] = seg_states[fil_i : fil_i+3]

                            for seg in range(1, int(self.pars['NSEG'])):
                                seg_pos = seg_states[fil_i+3*(seg-1) : fil_i+3*seg]
                                fil_data[seg] = seg_pos
                            ax.plot(fil_data[:,0], fil_data[:,1], fil_data[:,2], c='black', zorder = 100)

            plt.savefig(f'fig/ciliate_{nfil}fil.pdf', bbox_inches = 'tight', format='pdf')
            plt.show()
    
    
    # Rods
    def compute_rod_vel(self):
        body_states_f = open(simName + '_body_states.dat', "r")

        compute_start = self.plot_start_frame
        compute_start = 0

        body_states = np.zeros(7*int(self.pars['NSWIM']))
        average_vel = np.zeros((self.plot_end_frame - compute_start, 3))

        for i in range(self.plot_end_frame):
            print(" frame ", i, "/", self.frames, "          ", end="\r")
            body_states_str = body_states_f.readline()

            if(i == compute_start-1 or i == 0):
                body_states = np.array(body_states_str.split()[1:], dtype=float)

            if(i >= compute_start):
                body_states2 = body_states
                body_states = np.array(body_states_str.split()[1:], dtype=float)
                body_disp = body_states - body_states2
                
                for swim in range(int(self.pars['NSWIM'])):
                    body_vel = body_disp[7*swim : 7*swim+3]/self.dt
                    average_vel[i-compute_start] += body_vel

        time_array = np.arange(compute_start, self.plot_end_frame )
        average_vel /= float(self.pars['NSWIM'])
        ax = plt.figure().add_subplot(1,1,1)
        ax.set_ylabel(r"<$V_x$>/L")
        ax.set_xlabel(r"time")
        ax.plot(time_array[1:], average_vel[1:,0]/self.L)
        plt.savefig('fig/rod_velocity.eps', format='eps')

        ax2 = plt.figure().add_subplot(1,1,1)
        ax2.set_ylabel(r"<$|V_x|$>/L")
        ax2.set_xlabel(r"time")
        ax2.plot(time_array[1:], average_vel[1:,0]/self.L)
        plt.savefig('fig/rod_velocity_abs.eps', format='eps')
        plt.show()

        plt.show()

    def multi_rod_vel(self):
        ac = 1 - 1./4.*np.pi
        T0 = 14.14**2/22.
        ls = ['solid', 'dashed', 'dotted']
        markers = ["^", "s", "d"]

        simDir = 'data/2304rod_sims/'
        simNames = ['test_rod_960_960_60', 'test_rod_1920_1920_60', 'test_rod_3840_3840_60']
        box = [(960,960,960), (1920,1920,960), (3840,3840,3840)]
        area_fraction = 2304*(48-28*ac) / np.array([(960**2), (1920**2), (3840**2)])
        V0_list = np.array([3.9666615301e-01, 4.1704212613e-01, 4.3713151925e-01])

        # simDir = 'data/8100rod_sims/'
        # simNames = ['test_rod_1800_1800_75', 'test_rod_3600_3600_75', 'test_rod_7200_7200_75']
        # box = [(1800,1800,1800), (3600,3600,3600), (7200,7200,7200)]
        # area_fraction = 8100*(48-28*ac) / np.array([(1800**2), (3600**2), (7200**2)])
        # V0_list = np.array([4.0106966224e-01, 4.1725594718e-01, 4.3356909706e-01])
    
        fig1 = plt.figure()
        ax = fig1.add_subplot(1,1,1)
        for ni, name in enumerate(simNames):
            Lx, Ly, Lz = box[ni]
            simName = simDir + name
        
            body_states_f = open(simName + '_body_states.dat', "r")

            compute_start = self.plot_start_frame
            compute_start = 0

            body_states = np.zeros(7*int(self.pars['NSWIM']))
            average_vel = np.zeros((self.plot_end_frame - compute_start, 3))
            # average_abs_vel = np.zeros((self.plot_end_frame - compute_start, 3))

            for i in range(self.plot_end_frame):
                print(name, " frame ", i, "/", self.frames, "          ", end="\r")
                body_states_str = body_states_f.readline()

                if(i == compute_start-1 or i == 0):
                    body_states = np.array(body_states_str.split()[1:], dtype=float)

                if(i >= compute_start):
                    body_states2 = body_states
                    body_states = np.array(body_states_str.split()[1:], dtype=float)
                    body_disp = body_states - body_states2
                    
                    for swim in range(int(self.pars['NSWIM'])):
                        body_vel = body_disp[7*swim : 7*swim+3]
                        average_vel[i-compute_start] += body_vel
                        # average_abs_vel[i-compute_start] += np.abs(body_vel)

            time_array = np.arange(compute_start, self.plot_end_frame ) / T0
            average_vel =  average_vel / float(self.pars['NSWIM']) / self.dt
            
            ax.plot(time_array[1:], average_vel[1:,0]/ V0_list[ni], c='black', linestyle=ls[ni],label='Area fraction='+'{:.2f}'.format(area_fraction[ni]*100)+'%')
            if(ni>=1):
                for ti, t in enumerate(self.plot_multi_rod_frames):
                    if(ni==1):
                        ax.scatter(time_array[t], average_vel[t, 0]/ V0_list[ni], s=60, marker=markers[ti], facecolors='none', edgecolors='black')
                    if(ni==2):
                        ax.scatter(time_array[t], average_vel[t, 0]/ V0_list[ni], s=60, marker=markers[ti], facecolors='black', edgecolors='black')
        rect = patches.Rectangle((1556, 1.4), 138, 0.2, linewidth=1, edgecolor='black', facecolor='none', ls='dotted')
        ax.add_patch(rect)
        ax.legend()
        ax.set_ylabel(r"$<V_x>/W$")
        ax.set_xlabel(r"$t/T$")
        ax.set_xlim(left=0)
        nswim = int(self.pars['NSWIM'])
        plt.savefig(f'fig/multi_velocity_{nswim}rods.eps', format='eps')
        plt.savefig(f'fig/multi_velocity_{nswim}rods.png', format='png')
        plt.show()

    def body_vel(self):
        self.frames = 31
        self.plot_end_frame = self.frames
        simDir = 'data/phase_model/ishikawa/'
        simNames = ['test_bab_160fil_6000blob_6R_2torsion_k0v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k0.5v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k1v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k1.5v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k2v0']
        ls = ['solid', 'dashed', 'dotted']
        markers = ["^", "s", "d"]
        labels = ["k=0","k=0.5","k=1","k=1.5","k=2",]
        colors = ["black","red","green","blue","purple"]
        
        ax = plt.figure().add_subplot(1,1,1)
        ax.set_xlim(0, 1)
        ax.set_ylabel(r"$V_zT/L$")
        ax.set_xlabel(r"t/T")

        for ni, name in enumerate(simNames):
            simName = simDir + name
            print(simName)

            body_states_f = open(simName + '_body_states.dat', "r")
            body_vels_f = open(simName + '_body_vels.dat', "r")

            compute_start = self.plot_start_frame
            compute_start = 0

            L = 2.2*(20-1)
            

            body_states = np.zeros(7*int(self.pars['NSWIM']))
            body_vel_array = np.zeros((self.plot_end_frame - compute_start, 3))
            body_disp_array = np.zeros((self.plot_end_frame - compute_start, 7))

            for i in range(self.plot_end_frame):
                print(" frame ", i, "/", self.frames, "          ", end="\r")
                body_states_str = body_states_f.readline()
                body_vels_str = body_vels_f.readline()

                if(i >= compute_start):
                    body_states = np.array(body_states_str.split()[1:], dtype=float)
                    body_vels = np.array(body_vels_str.split(), dtype=float)
                    body_vel_array[i-compute_start] = body_vels[:3]
                    body_disp_array[i-compute_start] = body_states

            time_array = np.arange(compute_start, self.plot_end_frame )
            ax.plot(time_array/30., body_vel_array[:,2]/L, label=labels[ni], c=colors[ni])
        ax.legend()
        plt.savefig('fig/swimmer_velocity.pdf', format='pdf')
        plt.show()
    
    def avg_body_vel(self):
        self.frames = 31
        self.plot_end_frame = self.frames
        simDir = 'data/phase_model/ishikawa/'
        simNames = ['test_bab_160fil_6000blob_6R_2torsion_k-2v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k-1.5v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k-1.1v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k-0.9v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k-0.4v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k0v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k0.3v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k0.5v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k0.7v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k0.8v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k1v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k1.1v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k1.5v0',
                    'test_bab_160fil_6000blob_6R_2torsion_k2v0']
        ks = [-2, -1.5, -1.1, -0.9, -.4, 0, .3, 0.5, .7, .8, 1, 1.1, 1.5, 2]
        avg_vel_array = np.zeros(len(ks))
        
        ax = plt.figure().add_subplot(1,1,1)
        ax.set_ylabel(r"$V_zT/L$")
        ax.set_xlabel(r"k")

        for ni, name in enumerate(simNames):
            simName = simDir + name
            print(simName)

            body_states_f = open(simName + '_body_states.dat', "r")
            body_vels_f = open(simName + '_body_vels.dat', "r")

            compute_start = self.plot_start_frame
            compute_start = 0

            L = 2.2*(20-1)
            

            body_states = np.zeros(7*int(self.pars['NSWIM']))
            body_vel_array = np.zeros((self.plot_end_frame - compute_start, 3))
            body_disp_array = np.zeros((self.plot_end_frame - compute_start, 7))

            for i in range(self.plot_end_frame):
                print(" frame ", i, "/", self.frames, "          ", end="\r")
                body_states_str = body_states_f.readline()
                body_vels_str = body_vels_f.readline()

                if(i >= compute_start):
                    body_states = np.array(body_states_str.split()[1:], dtype=float)
                    body_vels = np.array(body_vels_str.split(), dtype=float)
                    body_vel_array[i-compute_start] = body_vels[:3]
                    body_disp_array[i-compute_start] = body_states
                
                avg_vel_array[ni] += body_vels[2]/L

        ax.plot(ks, avg_vel_array/30)
        plt.savefig('fig/velocity_vs_k.pdf', format='pdf')
        plt.show()

    def plot_hist(self):
        body_states_f = open(simName + '_body_states.dat', "r")
        body_states = np.zeros(7*int(self.pars['NSWIM']))
        body_vel_x = np.zeros(int(self.pars['NSWIM']))

        ax2 = plt.figure().add_subplot(1,1,1)
        ax2.set_ylabel(r"frequency")
        ax2.set_xlabel(r"$V_x$/L")
        bins = np.linspace(-8/self.L, 10/self.L, 20)

        for i in range(max(self.plot_hist_frame) + 1):
            print(" frame ", i, "/", self.frames, "          ", end="\r")
            if(i>self.frames-1):
                raise ValueError("Exceeding available frame number")
            body_states_str = body_states_f.readline()

            if(i in self.plot_hist_frame-1):
                body_states = np.array(body_states_str.split()[1:], dtype=float)

            if(i in self.plot_hist_frame):
                body_states2 = body_states
                body_states = np.array(body_states_str.split()[1:], dtype=float)
                body_disp = body_states - body_states2
                for swim in range(int(self.pars['NSWIM'])):
                    body_vel_x[swim] = body_disp[7*swim : 7*swim+1]/self.dt
                ax2.hist(body_vel_x/self.L, bins=bins, label=f'time={i}')
        ax2.legend()
        plt.savefig('fig/rod_velocity_hist.eps', format='eps')
        plt.show()

    def plot_rod_vel(self):

        body_states_f = open(simName + '_body_states.dat', "r")

        end_pos = np.zeros((int(self.pars['NSWIM']), 3))
        end_vel = np.zeros((int(self.pars['NSWIM']), 3))

        fig1 = plt.figure(figsize=(6, 6), dpi=300)
        fig2 = plt.figure()
        ax = fig1.add_subplot(1,1,1)
        ax2 = fig2.add_subplot(1,1,1)
        bins = np.linspace(0, Lx, 15)

        for i in range(self.plot_end_frame):
            print(" frame ", i, "/", self.frames, "          ", end="\r")
            body_states_str = body_states_f.readline()

            if(i in self.plot_rod_frame-1):
                body_states = np.array(body_states_str.split()[1:], dtype=float)

            if(i in self.plot_rod_frame):
                body_states2 = body_states
                body_states = np.array(body_states_str.split()[1:], dtype=float)
                body_disp = body_states - body_states2

                x_ = body_states[::7]
                
                for swim in range(int(self.pars['NSWIM'])):
                    end_pos[swim] = body_states[7*swim : 7*swim+3]
                    end_vel[swim] = body_disp[7*swim : 7*swim+3]/self.dt

                    rod_head = np.array(util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[:3]))
                    rod_tail = np.array(util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[-3:]))

                    pdiff = rod_tail - rod_head
                    two_points = np.array([util.box(rod_head, np.array([Lx, Ly, Lz])), util.box(rod_head, np.array([Lx, Ly, Lz])) + pdiff])

                    ax.plot(two_points[:,0], two_points[:,1], c= 'r', zorder=0, lw=.3
                    )

                if self.enable_superpunto:
                    end_pos = util.box(end_pos, np.array([Lx, Ly, Lz]))

                ax2.hist(end_pos[:,1], bins=bins, label=f'time={i}')
        
        # ax.scatter(end_pos[:,0], end_pos[:,1])
        ax.quiver(end_pos[:,0], end_pos[:,1], end_vel[:,0], end_vel[:,1], zorder=1)
        ax.set_ylabel(r"y")
        ax.set_xlabel(r"x")
        ax.set_xlim(0, Lx)
        ax.set_ylim(0, Ly)
        ax.set_aspect('equal')
        # ax.axis('off')
        nswim = int(self.pars['NSWIM'])
        fig1.savefig(f'fig/rod_displacement_{nswim}rods_{int(Lx)}_{int(Ly)}_{int(Lz)}_{int(self.plot_end_frame)}.eps', bbox_inches = 'tight',  format='eps')
        fig1.savefig(f'fig/rod_displacement_{nswim}rods_{int(Lx)}_{int(Ly)}_{int(Lz)}_{int(self.plot_end_frame)}.pdf', bbox_inches = 'tight',  format='pdf')

        ax2.set_ylabel(r"frequency")
        ax2.set_xlabel(r"$y$")
        ax2.set_xlim(0, Lx)
        fig2.savefig(f'fig/pos_distribution_{nswim}rods_{int(Lx)}_{int(Ly)}_{int(Lz)}_{int(self.plot_end_frame)}.eps', bbox_inches = 'tight',  format='eps')

        plt.show()

    def plot_rod_vel_multi(self):

        simDir = 'data/2304rod_sims/'
        simNames = ['test_rod_1920_1920_60', 'test_rod_3840_3840_60']
        box = [(1920,1920,960), (3840,3840,3840)]
        markers = ["^", "s", "d"]
        marker_string = [r'$\Delta$', r'$square$', r'$\bigcirc$']
        facecolors = ['none', 'black']
        
        fig, axs = plt.subplots(3, 2, figsize=(8, 9))
        for ni, name in enumerate(simNames):
            Lx, Ly, Lz = box[ni]
            simName = simDir + name
            body_states_f = open(simName + '_body_states.dat', "r")

            end_pos = np.zeros((int(self.pars['NSWIM']), 3))
            end_vel = np.zeros((int(self.pars['NSWIM']), 3))

            current_fig = 0
            for i in range(self.plot_end_frame):
                
                print(name, " frame ", i, "/", self.frames, "          ", end="\r")
                body_states_str = body_states_f.readline()

                if(i in self.plot_multi_rod_frames-1):
                    body_states = np.array(body_states_str.split()[1:], dtype=float)

                if(i in self.plot_multi_rod_frames):
                    body_states2 = body_states
                    body_states = np.array(body_states_str.split()[1:], dtype=float)
                    body_disp = body_states - body_states2

                    x_ = body_states[::7]
                    
                    for swim in range(int(self.pars['NSWIM'])):
                        end_pos[swim] = body_states[7*swim : 7*swim+3]
                        end_vel[swim] = body_disp[7*swim : 7*swim+3]/self.dt

                        rod_head = np.array(util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[:3]))
                        rod_tail = np.array(util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[-3:]))

                        pdiff = rod_tail - rod_head
                        two_points = np.array([util.box(rod_head, np.array([Lx, Ly, Lz])), util.box(rod_head, np.array([Lx, Ly, Lz])) + pdiff])

                        axs[current_fig, ni].plot(two_points[:,0], two_points[:,1], c= 'r', zorder=0, lw=.3)

                    if self.enable_superpunto:
                        end_pos = util.box(end_pos, np.array([Lx, Ly, Lz]))

                    axs[current_fig, ni].scatter(-1, -1, marker=markers[current_fig], facecolors=facecolors[ni], edgecolors='black', label=' ')
                    axs[current_fig, ni].quiver(end_pos[:,0], end_pos[:,1], end_vel[:,0], end_vel[:,1], zorder=1)
                    axs[current_fig, ni].set_xlim(0, Lx)
                    axs[current_fig, ni].set_ylim(0, Ly)
                    axs[current_fig, ni].set_aspect('equal')
                    axs[current_fig, ni].legend(loc='upper left', frameon=False)
                    current_fig +=1
            nswim = int(self.pars['NSWIM'])
        for axi, ax in enumerate(axs.flat):
            ax.set(ylabel=r"y", xlabel=r"x")
        
        T = 14.14**2/22.
        times = [fr"t={int(self.plot_multi_rod_frames[i]/T)} T" for i in range(3)]
        pad = 5 # in points
        for ax, time in zip(axs[:,0], times):
            ax.annotate(time, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                        xycoords=ax.yaxis.label, textcoords='offset points',
                        size='large', ha='right', va='center')
        fig.tight_layout()
        fig.subplots_adjust(left=0.21, top=0.95)
        # Hide x labels and tick labels for top plots and y ticks for right plots.
        # for ax in axs.flat:
        #     ax.label_outer()
        
        # fig.supxlabel(r"x")
        # fig.supylabel(r"y")

        fig.savefig(f'fig/rod_displacement_{nswim}rods_multi.eps', bbox_inches = 'tight',  format='eps')
        fig.savefig(f'fig/rod_displacement_{nswim}rods_multi.pdf', bbox_inches = 'tight',  format='pdf')

        plt.show()