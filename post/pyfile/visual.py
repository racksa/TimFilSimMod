import numpy as np
import matplotlib.pyplot as plt
import myIo
import util
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
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

simDir = 'data/64fils_800_800_200/'
simDir = 'data/rod_sims/'
# simDir = ''
simName = simDir + 'test_rod'
superpuntoDatafileName = '../../' + simName + '_superpunto.dat'
fcmPosfileName = '../../' + simName + '_flow_pos.dat'
fcmForcefileName = '../../' + simName + '_flow_force.dat'
fcmTorquefileName = '../../' + simName + '_flow_torque.dat'

Lx = 800.
Ly = 800.
Lz = 800.*2

# rod_16384
# Lx = 2048
# Ly = 2048
# Lz = 32

# rod_7744
# Lx = 1760.
# Ly = 1760.
# Lz = 10.3125*10
# 1024*1024*6

# fil_4096
# Lx = 3840.
# Ly = 3840.
# Lz = 240.
# 256*256*16

#fil_4096_cube
# Lx = 6400.
# Ly = 6400.
# Lz = 6400.

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

enable_periodic = True

class VISUAL:

    def __init__(self):
        self.pars = myIo.read_pars('../../' + simName + '.par')
        if(self.pars['NBLOB']>0):
            self.blob_references = myIo.read_blob_references('../../' + simName + '_blob_references.dat')
        if(self.pars['NFIL']>0):
            self.fil_references = myIo.read_fil_references('../../' + simName + '_fil_references.dat')
        self.dt = self.pars['DT']*self.pars['PLOT_FREQUENCY_IN_STEPS']
        self.L = 14.14*self.pars['NBLOB']/22.
        self.frames = sum(1 for line in open('../../' + simName + '_body_states.dat'))

        self.plot_start_frame = max(0, self.frames-100)
        self.plot_end_frame = self.frames
        self.plot_interval = 5

        self.plot_hist_frame = np.array([self.frames-1])
        self.plot_seg_frame = self.plot_end_frame-1
        self.plot_rod_frame = self.plot_end_frame-1
        self.fcm_frame = self.plot_end_frame-1

        self.output_to_superpunto = False
        self.output_to_fcm = False

    def set_plot_dim(self, dim):
        self.plot_dim = dim

    def enable_superpunto(self):
        self.output_to_superpunto = True
        myIo.clean_file(superpuntoDatafileName)

    def enable_fcm(self):
        self.output_to_fcm = True
        myIo.clean_file(fcmPosfileName)
        myIo.clean_file(fcmForcefileName)
        myIo.clean_file(fcmTorquefileName)

    def write_data(self, x, r, filename, box=True, center=True, superpunto=True):
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
                            str(0),
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
        seg_states_f = open('../../' + simName + '_seg_states.dat', "r")
        body_states_f = open('../../' + simName + '_body_states.dat', "r")
        if(self.pars['NBLOB']>0):
            blob_forces_f = open('../../' + simName + '_blob_forces.dat', "r")
        seg_forces_f = open('../../' + simName + '_seg_forces.dat', "r")
        print("open files time = ",(time.time()-start))

        for i in range(self.plot_end_frame):
            print("frame ", i, "/", self.frames, end="\r")
            body_states_str = body_states_f.readline()
            if(self.pars['NFIL']>0):
                seg_states_str = seg_states_f.readline()
                if(self.output_to_fcm):
                    seg_forces_str = seg_forces_f.readline()
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
                        # To find blob position
                        for blob in range(int(self.pars['NBLOB'])):
                            blob_x, blob_y, blob_z = util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[3*blob:3*blob+3])
                            if(self.output_to_superpunto):
                                self.write_data([blob_x, blob_y, blob_z], float(self.pars['RBLOB']), superpuntoDatafileName, enable_periodic)
                            if(self.output_to_fcm and i == self.fcm_frame):
                                blob_fx, blob_fy, blob_fz = blob_forces[3*blob : 3*blob+3]
                                self.write_data([blob_x, blob_y, blob_z], 0, fcmPosfileName, box=True, center=False, superpunto=False)
                                self.write_data([blob_fx, blob_fy, blob_fz], 0, fcmForcefileName, box=False, center=False, superpunto=False)
                                self.write_data([0, 0, 0], 0, fcmTorquefileName, box=False, center=False, superpunto=False)

                        # Robot arm to find segment position (Ignored plane rotation!)
                        for fil in range(int(self.pars['NFIL'])):
                            
                            fil_i = int(4*fil*self.pars['NSEG'])
                            fil_base_x, fil_base_y, fil_base_z = body_pos + np.matmul(R, self.fil_references[3*fil : 3*fil+3])
                            old_seg_pos = np.array([fil_base_x, fil_base_y, fil_base_z])
                            if(self.output_to_superpunto):
                                self.write_data(old_seg_pos, float(self.pars['RSEG']), superpuntoDatafileName, enable_periodic, True, True)
                            if(self.output_to_fcm):
                                seg_fx, seg_fy, seg_fz, seg_tx, seg_ty, seg_tz = seg_forces[0 : 6]
                                self.write_data(old_seg_pos, 0, fcmPosfileName, box=True, center=False, superpunto=False)
                                self.write_data([seg_fx, seg_fy, seg_fz], 0, fcmForcefileName, box=False, center=False, superpunto=False)
                                self.write_data([seg_tx, seg_ty, seg_tz], 0, fcmTorquefileName, box=False, center=False, superpunto=False)

                            for seg in range(1, int(self.pars['NSEG'])):
                                q1 = seg_states[fil_i+4*(seg-1) : fil_i+4*seg]
                                q2 = seg_states[fil_i+4*seg : fil_i+4*seg+4]
                                
                                t1 = util.find_t(q1)
                                t2 = util.find_t(q2)
                                
                                seg_pos = old_seg_pos + 0.5*self.pars['DL']*(t1 + t2)
                                old_seg_pos = seg_pos
                                
                                if(self.output_to_superpunto):
                                    self.write_data(seg_pos, float(self.pars['RSEG']), superpuntoDatafileName, enable_periodic, True, True)
                                if(self.output_to_fcm):
                                    seg_fx, seg_fy, seg_fz, seg_tx, seg_ty, seg_tz = seg_forces[6*seg : 6*seg+6]
                                    self.write_data(seg_pos, 0, fcmPosfileName, box=True, center=False, superpunto=False)
                                    self.write_data([seg_fx, seg_fy, seg_fz], 0, fcmForcefileName, box=False, center=False, superpunto=False)
                                    self.write_data([seg_tx, seg_ty, seg_tz], 0, fcmTorquefileName, box=False, center=False, superpunto=False)


        if(self.output_to_fcm):
            fcm_directory = "../../../CUFCM/data/flow_data/"
            subprocess.call("cp " + fcmPosfileName + " " + fcm_directory, shell=True)
            subprocess.call("cp " + fcmForcefileName + " " + fcm_directory, shell=True)
            subprocess.call("cp " + fcmTorquefileName + " " + fcm_directory, shell=True)

    def compute_rod_vel(self):
        body_states_f = open('../../' + simName + '_body_states.dat', "r")

        body_states = np.zeros(7*int(self.pars['NSWIM']))
        average_vel = np.zeros((self.plot_end_frame-self.plot_start_frame, 3))

        for i in range(self.plot_end_frame):
            print("frame ", i, "/", self.frames, end="\r")
            body_states_str = body_states_f.readline()

            if(i == self.plot_start_frame-1 or i == 0):
                body_states = np.array(body_states_str.split()[1:], dtype=float)

            if(i >= self.plot_start_frame):
                body_states2 = body_states
                body_states = np.array(body_states_str.split()[1:], dtype=float)
                body_disp = body_states - body_states2
                
                for swim in range(int(self.pars['NSWIM'])):
                    body_vel = body_disp[7*swim : 7*swim+3]/self.dt
                    average_vel[i-self.plot_start_frame] += body_vel

        time_array = np.arange(self.plot_start_frame, self.plot_end_frame )
        average_vel /= float(self.pars['NSWIM'])
        ax = plt.figure().add_subplot(1,1,1)
        ax.set_ylabel(r"<$V_x$>/L")
        ax.set_xlabel(r"time")
        ax.plot(time_array[1:], average_vel[1:,0]/self.L)
        plt.savefig('fig/rod_velocity.eps', format='eps')
        plt.show()

    def plot_hist(self):
        body_states_f = open('../../' + simName + '_body_states.dat', "r")
        body_states = np.zeros(7*int(self.pars['NSWIM']))
        body_vel_x = np.zeros(int(self.pars['NSWIM']))

        ax2 = plt.figure().add_subplot(1,1,1)
        ax2.set_ylabel(r"frequency")
        ax2.set_xlabel(r"$V_x$/L")
        bins = np.linspace(-8/self.L, 10/self.L, 20)

        for i in range(max(self.plot_hist_frame) + 1):
            print("frame ", i, "/", self.frames, end="\r")
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

    def plot_seg_vel(self):
        seg_states_f = open('../../' + simName + '_seg_states.dat', "r")
        seg_vels_f = open('../../' + simName + '_seg_vels.dat', "r")
        body_states_f = open('../../' + simName + '_body_states.dat', "r")

        ax = plt.figure().add_subplot(1,1,1)
        end_pos = np.zeros((int(self.pars['NFIL']), 3))
        end_vel = np.zeros((int(self.pars['NFIL']), 3))

        for i in range(self.plot_end_frame):
            print("frame ", i, "/", self.frames, end="\r")
            body_states_str = body_states_f.readline()
            if(self.pars['NFIL']>0):
                seg_states_str = seg_states_f.readline()
                seg_vels_str = seg_vels_f.readline()

            if(i == self.plot_seg_frame):
                seg_states = np.array(seg_states_str.split()[1:], dtype=float)
                seg_vels = np.array(seg_vels_str.split()[1:], dtype=float)
                body_states = np.array(body_states_str.split()[1:], dtype=float)

                R = util.rot_mat(body_states[3 : 7])

                # Robot arm to find segment position (Ignored plane rotation!)
                for fil in range(int(self.pars['NFIL'])):
                    fil_i = int(4*fil*self.pars['NSEG'])
                    fil_base_x, fil_base_y, fil_base_z = np.matmul(R, self.fil_references[3*fil : 3*fil+3])
                    old_seg_pos = np.array([fil_base_x, fil_base_y, fil_base_z])

                    for seg in range(1, int(self.pars['NSEG'])):
                        q1 = seg_states[fil_i+4*(seg-1) : fil_i+4*seg]
                        q2 = seg_states[fil_i+4*seg : fil_i+4*seg+4]
                        
                        t1 = util.find_t(q1)
                        t2 = util.find_t(q2)
                        
                        seg_pos = old_seg_pos + 0.5*self.pars['DL']*(t1 + t2)
                        old_seg_pos = seg_pos
                    end_pos[fil] = seg_pos
                    end_vel[fil] = seg_vels[6*int(self.pars['NSEG'])*fil + 6*int(self.pars['NSEG']) - 6 : 6*int(self.pars['NSEG'])*fil + 6*int(self.pars['NSEG']) - 3]
            
        ax.scatter(end_pos[:,0], end_pos[:,1])
        ax.quiver(end_pos[:,0], end_pos[:,1], end_vel[:,0], end_vel[:,1])
        plt.savefig('fig/fil_vel.eps', format='eps')
        plt.show()

    def plot_rod_vel(self):
        body_states_f = open('../../' + simName + '_body_states.dat', "r")

        end_pos = np.zeros((int(self.pars['NSWIM']), 3))
        end_vel = np.zeros((int(self.pars['NSWIM']), 3))

        ax = plt.figure().add_subplot(1,1,1)
        ax2 = plt.figure().add_subplot(1,1,1)
        bins = np.linspace(0, Lx, 20)
        
        for i in range(self.plot_end_frame):
            print("frame ", i, "/", self.frames, end="\r")
            body_states_str = body_states_f.readline()

            if(i == self.plot_rod_frame-1):
                body_states = np.array(body_states_str.split()[1:], dtype=float)

            if(i == self.plot_rod_frame):
                body_states2 = body_states
                body_states = np.array(body_states_str.split()[1:], dtype=float)
                body_disp = body_states - body_states2
                
                for swim in range(int(self.pars['NSWIM'])):
                    end_pos[swim] = body_states[7*swim : 7*swim+3]
                    end_vel[swim] = body_disp[7*swim : 7*swim+3]/self.dt

                    rod_head = np.array(util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[:3]))
                    rod_tail = np.array(util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[-3:]))
                    pdiff = rod_tail - rod_head
                    two_points = np.array([util.box(rod_head, np.array([Lx, Ly, Lz])), util.box(rod_head, np.array([Lx, Ly, Lz])) + pdiff])

                    ax.plot(two_points[:,0], two_points[:,1], c= 'b', zorder=0)

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
        plt.savefig('fig/rod_vel.eps', format='eps')

        ax2.set_ylabel(r"frequency")
        ax2.set_xlabel(r"$y$")
        ax2.set_xlim(0, Lx)
        plt.savefig('fig/pos_distribution.eps', format='eps')

        plt.show()