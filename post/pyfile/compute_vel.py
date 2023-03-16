import numpy as np
import matplotlib.pyplot as plt
import myIo
import util
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import fileinput
import random
import time

overlap = {name for name in mcolors.CSS4_COLORS
           if f'xkcd:{name}' in mcolors.XKCD_COLORS}
color_list = list(overlap)
color_list.sort()
color_list.pop(2)
color_list.pop(2)

# color_list = ['#00ffff', '#faebD7', '#838bbb', '#0000ff', '	#8a2be2', '#ff4040', '#7fff00', '#ff6103', '#9932cc', '#ff1493', '#030303']

simName = 'test_fil'
# superpuntoDatafileName = '../../' + simName + '_superpunto.dat'

enable_periodic = True

class COMPUTEVEL:

    def __init__(self):
        self.pars = myIo.read_pars('../../' + simName + '.par')
        if(self.pars['NBLOB']>0):
            self.blob_references = myIo.read_blob_references('../../' + simName + '_blob_references.dat')
        if(self.pars['NFIL']>0):
            self.fil_references = myIo.read_fil_references('../../' + simName + '_fil_references.dat')
        
        self.frames = sum(1 for line in open('../../' + simName + '_body_states.dat'))
        self.dt = self.pars['DT']*self.pars['PLOT_FREQUENCY_IN_STEPS']
        self.L = 14.14

        self.plot_start_frame = 0
        self.plot_end_frame = self.frames

        self.plot_hist_frame = [10, 50, 100, 200, 580]
        self.plot_seg_frame = self.plot_end_frame-1

    def compute(self):
        body_states_f = open('../../' + simName + '_body_states.dat', "r")

        body_states = np.zeros(7*int(self.pars['NSWIM']))
        average_vel = np.zeros((self.plot_end_frame-self.plot_start_frame, 3))
        
        for i in range(self.plot_start_frame):
            body_states_f.readline()

        for i in range(self.plot_start_frame, self.plot_end_frame):
            print("frame ", i, "/", self.frames, flush=True)
            body_states2 = body_states
            body_states = np.array(body_states_f.readline().split()[1:], dtype=float)
            if(i==0):
                body_states2 = body_states
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
            body_states2 = body_states
            body_states = np.array(body_states_f.readline().split()[1:], dtype=float)
            if(i in self.plot_hist_frame):
                print(i)
                body_disp = body_states - body_states2
                for swim in range(int(self.pars['NSWIM'])):
                    body_vel_x[swim] = body_disp[7*swim : 7*swim+1]/self.dt
                ax2.hist(body_vel_x/self.L, bins=bins, label=f'time={i}')
        ax2.legend()
        plt.savefig('fig/rod_velocity_hist.eps', format='eps')
        plt.show()


    def plot_seg_height(self):
        seg_states_f = open('../../' + simName + '_seg_states.dat', "r")
        seg_vels_f = open('../../' + simName + '_seg_vels.dat', "r")
        body_states_f = open('../../' + simName + '_body_states.dat', "r")

        ax = plt.figure().add_subplot(1,1,1)
        end_pos = np.zeros((int(self.pars['NFIL']), 3))
        end_vel = np.zeros((int(self.pars['NFIL']), 3))

        for i in range(self.plot_end_frame):
            print("frame ", i, "/", self.frames, flush=True)
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

        # time_array = np.arange(self.plot_start_frame, self.plot_end_frame )
        # average_vel /= float(self.pars['NSWIM'])
        # ax = plt.figure().add_subplot(1,1,1)
        # ax.set_ylabel(r"<$V_x$>/L")
        # ax.set_xlabel(r"time")
        # ax.plot(time_array[1:], average_vel[1:,0]/self.L)
        plt.savefig('fig/fil_height.eps', format='eps')
        plt.show()


















#