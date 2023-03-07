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

simName = 'test_rod_1024'
# superpuntoDatafileName = '../../' + simName + '_superpunto.dat'

Lx = 2560.
Ly = 2560.
Lz = 2560.

# rod_16384
# Lx = 2048
# Ly = 2048
# Lz = 32

# rod_8100
# Lx = 1440.
# Ly = 1440.
# Lz = 22.5

# rod_7744
# Lx = 1760.
# Ly = 1760.
# Lz = 10.3125*10
# 1024*1024*6

# rod_4096 (actually 1024 fils)
# Lx = 2560.
# Ly = 2560.
# Lz = 120.

# rod_1024
Lx = 640.
Ly = 640.
Lz = 10.*10
# 512*512*8

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

        self.plot_start_frame = 0
        self.plot_end_frame = 464

        self.plot_hist_frame = [10, 50, 100]

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
        ax.set_ylabel(r"<$V_x$>")
        ax.set_xlabel(r"time")
        ax.plot(time_array[1:], average_vel[1:,0])
        plt.show()

    def plot_hist(self):
        body_states_f = open('../../' + simName + '_body_states.dat', "r")
        body_states = np.zeros(7*int(self.pars['NSWIM']))
        body_vel_x = np.zeros(int(self.pars['NSWIM']))

        ax2 = plt.figure().add_subplot(1,1,1)
        ax2.set_ylabel(r"frequency")
        ax2.set_xlabel(r"$V_x$")
        bins = np.linspace(-5, 6, 100)

        for i in range(max(self.plot_hist_frame) + 1):
            body_states2 = body_states
            body_states = np.array(body_states_f.readline().split()[1:], dtype=float)
            if(i in self.plot_hist_frame):
                print(i)
                body_disp = body_states - body_states2
                for swim in range(int(self.pars['NSWIM'])):
                    body_vel_x[swim] = body_disp[7*swim : 7*swim+1]/self.dt
                ax2.hist(body_vel_x, bins=bins, label=f'time={i}')
        ax2.legend()
        plt.show()






















#