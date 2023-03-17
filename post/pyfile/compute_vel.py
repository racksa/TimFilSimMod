# import numpy as np
# import matplotlib.pyplot as plt
# import myIo
# import util
# import pandas as pd
# import matplotlib.colors as mcolors
# from matplotlib.lines import Line2D
# import fileinput
# import random
# import time

# simDir = 'data/rod_sims/'
# # simDir = ''
# simName = simDir + 'test_rod'
# # superpuntoDatafileName = '../../' + simName + '_superpunto.dat'

# enable_periodic = True

# class COMPUTEVEL:

#     def __init__(self):
#         self.pars = myIo.read_pars('../../' + simName + '.par')
#         if(self.pars['NBLOB']>0):
#             self.blob_references = myIo.read_blob_references('../../' + simName + '_blob_references.dat')
#         if(self.pars['NFIL']>0):
#             self.fil_references = myIo.read_fil_references('../../' + simName + '_fil_references.dat')
        
#         self.frames = sum(1 for line in open('../../' + simName + '_body_states.dat'))
#         self.dt = self.pars['DT']*self.pars['PLOT_FREQUENCY_IN_STEPS']
#         self.L = 14.14

#         self.plot_start_frame = 0
#         self.plot_end_frame = self.frames

#         self.plot_hist_frame = [10, 50, 100, 200, 580]
#         self.plot_seg_frame = self.plot_end_frame-1

#     def compute(self):
#         body_states_f = open('../../' + simName + '_body_states.dat', "r")

#         body_states = np.zeros(7*int(self.pars['NSWIM']))
#         average_vel = np.zeros((self.plot_end_frame-self.plot_start_frame, 3))
        
#         for i in range(self.plot_start_frame):
#             body_states_f.readline()

#         for i in range(self.plot_start_frame, self.plot_end_frame):
#             print("frame ", i, "/", self.frames, flush=True)
#             body_states2 = body_states
#             body_states = np.array(body_states_f.readline().split()[1:], dtype=float)
#             if(i==0):
#                 body_states2 = body_states
#             body_disp = body_states - body_states2
            
#             for swim in range(int(self.pars['NSWIM'])):
#                 body_vel = body_disp[7*swim : 7*swim+3]/self.dt
#                 average_vel[i-self.plot_start_frame] += body_vel

#         time_array = np.arange(self.plot_start_frame, self.plot_end_frame )
#         average_vel /= float(self.pars['NSWIM'])
#         ax = plt.figure().add_subplot(1,1,1)
#         ax.set_ylabel(r"<$V_x$>/L")
#         ax.set_xlabel(r"time")
#         ax.plot(time_array[1:], average_vel[1:,0]/self.L)
#         plt.savefig('fig/rod_velocity.eps', format='eps')
#         plt.show()

#     def plot_hist(self):
#         body_states_f = open('../../' + simName + '_body_states.dat', "r")
#         body_states = np.zeros(7*int(self.pars['NSWIM']))
#         body_vel_x = np.zeros(int(self.pars['NSWIM']))

#         ax2 = plt.figure().add_subplot(1,1,1)
#         ax2.set_ylabel(r"frequency")
#         ax2.set_xlabel(r"$V_x$/L")
#         bins = np.linspace(-8/self.L, 10/self.L, 20)

#         for i in range(max(self.plot_hist_frame) + 1):
#             body_states2 = body_states
#             body_states = np.array(body_states_f.readline().split()[1:], dtype=float)
#             if(i in self.plot_hist_frame):
#                 print(i)
#                 body_disp = body_states - body_states2
#                 for swim in range(int(self.pars['NSWIM'])):
#                     body_vel_x[swim] = body_disp[7*swim : 7*swim+1]/self.dt
#                 ax2.hist(body_vel_x/self.L, bins=bins, label=f'time={i}')
#         ax2.legend()
#         plt.savefig('fig/rod_velocity_hist.eps', format='eps')
#         plt.show()





















#