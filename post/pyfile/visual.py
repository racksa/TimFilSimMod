import numpy as np
import matplotlib.pyplot as plt
import myIo
import util
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import subprocess
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

simName = 'test_fil_4096'
superpuntoDatafileName = '../../' + simName + '_superpunto.dat'
fcmPosfileName = '../../' + simName + '_flow_pos.dat'
fcmForcefileName = '../../' + simName + '_flow_force.dat'
fcmTorquefileName = '../../' + simName + '_flow_torque.dat'

Lx = 240.
Ly = 240.
Lz = 120.

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
Lx = 3840.
Ly = 3840.
Lz = 240.

# rod_1024
# Lx = 640.
# Ly = 640.
# Lz = 10.*10
# 512*512*8

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
        
        self.frames = sum(1 for line in open('../../' + simName + '_body_states.dat'))

        self.plot_start_frame = 0
        self.plot_end_frame = self.frames
        self.plot_interval = 1

        self.fcm_frame = self.frames-1

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
        plot_x = np.zeros(np.shape(x))
        if(box):
            plot_x[0] = util.box(x[0], Lx)
            plot_x[1] = util.box(x[1], Ly)
            plot_x[2] = util.box(x[2], Lz)
            if(center):
                plot_x[0] -= 0.5*Lx
                plot_x[1] -= 0.5*Ly
        if superpunto:
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

        for i in range(0, self.plot_end_frame):
            body_states = np.array(body_states_f.readline().split()[1:], dtype=float)
            seg_states = np.array(seg_states_f.readline().split()[1:], dtype=float)
            if(self.pars['NBLOB']>0):
                blob_forces = np.array(blob_forces_f.readline().split()[1:], dtype=float)
            seg_forces = np.array(seg_forces_f.readline().split()[1:], dtype=float)

            if(i%self.plot_interval==0):
                print("frame ", i, "/", self.frames, flush=True)
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
                            self.write_data(old_seg_pos, float(self.pars['RSEG']), superpuntoDatafileName, enable_periodic, False, True)

                            for seg in range(1, int(self.pars['NSEG'])):
                                q1 = seg_states[fil_i+4*(seg-1) : fil_i+4*seg]
                                q2 = seg_states[fil_i+4*seg : fil_i+4*seg+4]
                                
                                t1 = util.find_t(q1)
                                t2 = util.find_t(q2)
                                
                                seg_pos = old_seg_pos + 0.5*self.pars['DL']*(t1 + t2)
                                old_seg_pos = seg_pos
                                
                                if(self.output_to_superpunto):
                                    self.write_data(seg_pos, float(self.pars['RSEG']), superpuntoDatafileName, enable_periodic, False, True)
                                if(self.output_to_fcm):
                                    pass
        if(self.output_to_fcm):
            fcm_directory = "../../../CUFCM/data/flow_data/"
            subprocess.call("cp " + fcmPosfileName + " " + fcm_directory, shell=True)
            subprocess.call("cp " + fcmForcefileName + " " + fcm_directory, shell=True)
            subprocess.call("cp " + fcmTorquefileName + " " + fcm_directory, shell=True)

            #         else:
            #             x0, y0, z0 = util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[:3])
            #             x1, y1, z1 = util.blob_point_from_data(body_states[7*swim : 7*swim+7], self.blob_references[-3:])
            #             x_diff, y_diff, z_diff = x1-x0, y1-y0, z1-z0
            #             two_points_x = [util.box(x0, Lx), util.box(x0, Lx) + x_diff]
            #             two_points_y = [util.box(y0, Ly), util.box(y0, Ly) + y_diff]
            #             two_points_z = [util.box(z0, Lz), util.box(z0, Lz) + z_diff]

            #             # color = (1-0.25*swim/self.pars['NSWIM'], 0.2*swim/self.pars['NSWIM'], 0.55+0.45*swim/self.pars['NSWIM'], (i+1-self.plot_start_frame)/(self.plot_end_frame-self.plot_start_frame) )
            #             random.seed(3*swim)
            #             cr = random.random()
            #             random.seed(3*swim+1)
            #             cg = random.random()
            #             random.seed(3*swim+2)
            #             cb = random.random()
            #             color = (cr, cg, cb, (i+1-self.plot_start_frame)/(self.plot_end_frame-self.plot_start_frame) )

            #             if (self.plot_dim == 2):
            #                 ax.scatter(two_points_x, two_points_y)
            #                 ax.plot(two_points_x, two_points_y, c=color)
            #             if (self.plot_dim == 3):
            #                 ax.plot(two_points_x, two_points_y, two_points_z, c=color )
                    
            # if (self.plot_dim == 3):
            #     util.set_axes_equal(ax)
            #     ax.set_xlabel("x")
            #     ax.set_ylabel("y")
            #     ax.set_zlabel("z")
            # if(not self.output_to_superpunto):
            #     plt.show()
