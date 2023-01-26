import numpy as np
import matplotlib.pyplot as plt
import myIo
import util
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import random

overlap = {name for name in mcolors.CSS4_COLORS
           if f'xkcd:{name}' in mcolors.XKCD_COLORS}
color_list = list(overlap)
color_list.sort()
color_list.pop(2)
color_list.pop(2)

# color_list = ['#00ffff', '#faebD7', '#838bbb', '#0000ff', '	#8a2be2', '#ff4040', '#7fff00', '#ff6103', '#9932cc', '#ff1493', '#030303']

simName = 'test_rod'
superpuntoDatafileName = '../../' + simName + '_superpunto.dat'

boxsize = 128


class VISUAL:

    def __init__(self):
        self.pars = myIo.read_pars('../../' + simName + '.par')
        if(self.pars['NBLOB']>0):
            self.blob_references = myIo.read_blob_references('../../' + simName + '_blob_references.dat')
        self.fil_references = myIo.read_fil_references('../../' + simName + '_fil_references.dat')
        self.body_states = myIo.read_body_states('../../' + simName + '_body_states.dat')
        self.seg_states = myIo.read_seg_states('../../' + simName + '_seg_states.dat')
        
        self.frames = len(self.body_states)

        self.plot_start_frame = 0
        self.plot_end_frame = 1000
        self.plot_interval = 1

        self.output_to_superpunto = False

    def set_plot_dim(self, dim):
        self.plot_dim = dim

    def enable_superpunto(self):
        self.output_to_superpunto = True
        myIo.clean_file(superpuntoDatafileName)

    def write_data(self, x, r, filename, box=True):
        if(box):
            x[0] = util.box(x[0], boxsize)
            x[1] = util.box(x[1], boxsize)
            x[2] = util.box(x[2], boxsize)
        myIo.write_line(str(x[0]) + ' ' +\
                        str(x[1]) + ' ' +\
                        str(x[2]) + ' ' +\
                        str(r) + ' ' +\
                        str(0),
                        filename)

    def plot(self):
        if (self.plot_dim == 2):
            ax = plt.figure().add_subplot(1,1,1)
            ax.axis('equal')
        if (self.plot_dim == 3):
            ax = plt.figure().add_subplot(projection='3d')

        for i in range(self.plot_start_frame, self.plot_end_frame, self.plot_interval):
            print("frame ", i, "/", self.frames)
            if(self.output_to_superpunto):
                myIo.write_line('#', superpuntoDatafileName)
            for swim in range(int(self.pars['NSWIM'])):
                body_pos = self.body_states[i][7*swim : 7*swim+3]
                R = util.rot_mat(self.body_states[i][7*swim+3 : 7*swim+7])
                
                if(self.output_to_superpunto):
                    # To find segment position
                    for blob in range(int(self.pars['NBLOB'])):
                        blob_x, blob_y, blob_z = util.blob_point_from_data(self.body_states[i][7*swim : 7*swim+7], self.blob_references[3*blob:3*blob+3])
                        self.write_data([blob_x, blob_y, blob_z], float(self.pars['RBLOB']), superpuntoDatafileName, False)

                    # Robot arm to find segment position (Ignored plane rotation!)
                    for fil in range(int(self.pars['NFIL'])):
                        fil_i = int(4*fil*self.pars['NSEG'])
                        fil_base_x, fil_base_y, fil_base_z = body_pos + np.matmul(R, self.fil_references[3*fil : 3*fil+3])
                        old_seg_pos = np.array([fil_base_x, fil_base_y, fil_base_z]) + 0.5*self.pars['DL']*util.find_t(self.seg_states[i][fil_i+0 : fil_i+4])     
                        for seg in range(1, int(self.pars['NSEG'])):
                            q1 = self.seg_states[i][fil_i+4*(seg-1) : fil_i+4*seg]
                            q2 = self.seg_states[i][fil_i+4*seg : fil_i+4*seg+4]
                            
                            t1 = util.find_t(q1)
                            t2 = util.find_t(q2)
                            
                            seg_pos = old_seg_pos + 0.5*self.pars['DL']*(t1 + t2)
                            old_seg_pos = seg_pos

                            # ax.scatter(seg_pos[0], seg_pos[1], seg_pos[2])
                            self.write_data(seg_pos, float(self.pars['RSEG']), superpuntoDatafileName, False)

                else:
                    x0, y0, z0 = util.blob_point_from_data(self.body_states[i][7*swim : 7*swim+7], self.blob_references[:3])
                    x1, y1, z1 = util.blob_point_from_data(self.body_states[i][7*swim : 7*swim+7], self.blob_references[-3:])
                    x_diff, y_diff, z_diff = x1-x0, y1-y0, z1-z0
                    two_points_x = [util.box(x0, boxsize), util.box(x0, boxsize) + x_diff]
                    two_points_y = [util.box(y0, boxsize), util.box(y0, boxsize) + y_diff]
                    two_points_z = [util.box(z0, boxsize), util.box(z0, boxsize) + z_diff]

                    # two_points_x = [x0, x1]
                    # two_points_y = [y0, y1]
                    # two_points_z = [z0, z1]

                    # color = (1-0.25*swim/self.pars['NSWIM'], 0.2*swim/self.pars['NSWIM'], 0.55+0.45*swim/self.pars['NSWIM'], (i+1-self.plot_start_frame)/(self.plot_end_frame-self.plot_start_frame) )
                    random.seed(3*swim)
                    cr = random.random()
                    random.seed(3*swim+1)
                    cg = random.random()
                    random.seed(3*swim+2)
                    cb = random.random()
                    color = (cr, cg, cb, (i+1-self.plot_start_frame)/(self.plot_end_frame-self.plot_start_frame) )

                    if (self.plot_dim == 2):
                        ax.scatter(two_points_x, two_points_y)
                        ax.plot(two_points_x, two_points_y, c=color)
                    if (self.plot_dim == 3):
                        ax.plot(two_points_x, two_points_y, two_points_z, c=color )
                 
        if (self.plot_dim == 3):
            util.set_axes_equal(ax)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")
        # if(not self.output_to_superpunto):
        # plt.show()
























#