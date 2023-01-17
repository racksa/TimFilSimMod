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

        self.blob_references = myIo.read_blob_references('../../' + simName + '_blob_references.dat')
        self.body_states = myIo.read_body_states('../../' + simName + '_body_states.dat')
        self.pars = myIo.read_pars('../../' + simName + '.par')
        self.frames = len(self.body_states)

        self.plot_start_frame = 0
        self.plot_end_frame = 60
        self.plot_interval = 1

        self.output_to_superpunto = False

    def set_plot_dim(self, dim):
        self.plot_dim = dim

    def enable_superpunto(self):
        self.output_to_superpunto = True
        myIo.clean_file(superpuntoDatafileName)

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

                if(self.output_to_superpunto):
                    for blob in range(int(self.pars['NBLOB'])):
                        blob_x, blob_y, blob_z = util.blob_point_from_data(self.body_states[i][7*swim : 7*swim+7], self.blob_references[3*blob:3*blob+3])
                        
                        myIo.write_line(str(util.box(blob_x, boxsize)) + ' ' +\
                                        str(util.box(blob_y, boxsize)) + ' ' +\
                                        str(util.box(blob_z, boxsize)) + ' ' +\
                                        str(float(self.pars['RBLOB'])) + ' ' +\
                                        str(0),
                                        superpuntoDatafileName)
                
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

        plt.show()
























#