import numpy as np
import matplotlib.pyplot as plt
import myIo
import util
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D   

overlap = {name for name in mcolors.CSS4_COLORS
           if f'xkcd:{name}' in mcolors.XKCD_COLORS}
color_list = list(overlap)
color_list.sort()
color_list.pop(2)
color_list.pop(2)

# color_list = ['#00ffff', '#faebD7', '#838bbb', '#0000ff', '	#8a2be2', '#ff4040', '#7fff00', '#ff6103', '#9932cc', '#ff1493', '#030303']

simName = 'test_rod'


class VISUAL:

    def __init__(self):

        self.blob_references = myIo.read_blob_references('../../' + simName + '_blob_references.dat')
        self.body_states = myIo.read_body_states('../../' + simName + '_body_states.dat')
        self.pars = myIo.read_pars('../../' + simName + '.par')
        self.frames = len(self.body_states)

    def set_plot_dim(self, dim):
        self.plot_dim = dim


    def plot(self):
        xmin, xmax, ymin, ymax, zmin, zmax = 0, 0, 0, 0, 0, 0
        def expand_plot_limit(new_x, new_y, new_z):
            xmin = min(xmin, new_x)
            xmax = max(xmax, new_x)
            ymin = min(ymin, new_y)
            ymax = max(ymax, new_y)
            zmin = min(zmin, new_z)
            zmax = max(zmax, new_z)

        if (self.plot_dim == 2):
            ax = plt.figure().add_subplot(1,1,1)
            # ax.axis('equal')
        if (self.plot_dim == 3):
            ax = plt.figure().add_subplot(projection='3d')

        for i in range(0, 100, 5):
            print("frame ", i, "/", self.frames)
            for swim in range(int(self.pars['NSWIM'])):
                line_x = list()
                line_y = list()
                line_z = list()
                for blob in range(int(self.pars['NBLOB'])):
                    x, y, z = util.blob_point_from_data(self.body_states[i][7*swim : 7*swim+7], self.blob_references[3*blob : 3*blob+3])
                    if (blob == 0) or (blob == int(self.pars['NBLOB'])-1):
                        line_x.append(x)
                        line_y.append(y)
                        line_z.append(z)

                if (self.plot_dim == 2):
                    # ax.scatter(line_x, line_z, c=color_list[swim], s=1, alpha=0.2+0.8*(i/self.frames))
                    ax.plot(line_x, line_z, color=color_list[swim])
                if (self.plot_dim == 3):
                    # ax.scatter(line_x, line_y, line_z, c=color_list[swim])
                    ax.plot(line_x, line_y, line_z, c=color_list[swim])

                
                
        # if (self.plot_dim == 2):
        #     ax.set_xlim((-400, 500))
        if (self.plot_dim == 3):
            ax.set_ylim(ax.get_xlim())

        plt.show()

























#