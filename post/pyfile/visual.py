import numpy as np
import matplotlib.pyplot as plt
import myIo
import pandas as pd
# from io import *

color_list = ['r', 'g', 'b', 'y', 'black', 'cyan']

simName = 'test_rod'


class VISUAL:

    def __init__(self):

        self.blob_references = myIo.read_blob_references('../../' + simName + '_blob_references.dat')
        self.body_states = myIo.read_body_states('../../' + simName + '_body_states.dat')
        self.pars = myIo.read_pars('../../' + simName + '.par')
        self.frames = len(self.body_states)


    def plot2d(self):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        for i in range(0, self.frames, self.frames-1):
            for swim in range(int(self.pars['NSWIM'])):
                for blob in range(int(self.pars['NBLOB'])):
                    x = self.body_states[i][7*swim + 0] + self.blob_references[3*blob + 0]
                    y = self.body_states[i][7*swim + 1] + self.blob_references[3*blob + 1]
                    z = self.body_states[i][7*swim + 2] + self.blob_references[3*blob + 2]
                    ax.scatter(x, z, c=color_list[swim%len(color_list)], alpha=0.5+0.5*(i/self.frames))
        
        ax.axis('equal')

        plt.show()


    def plot3d(self):
        ax = plt.figure().add_subplot(projection='3d')

        for i in range(0, self.frames, self.frames-1):
            for swim in range(int(self.pars['NSWIM'])):
                for blob in range(int(self.pars['NBLOB'])):
                    x = self.body_states[i][7*swim + 0] + self.blob_references[3*blob + 0]
                    y = self.body_states[i][7*swim + 1] + self.blob_references[3*blob + 1]
                    z = self.body_states[i][7*swim + 2] + self.blob_references[3*blob + 2]
                    ax.scatter(x, y, z, s=self.pars['RBLOB']*100, c=color_list[i%len(color_list) ])

        plt.show()






























#