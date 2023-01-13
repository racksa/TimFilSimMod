import os
import visual
import sys

visualiser = visual.VISUAL()

if(sys.argv[1] == 'plot2d'):
    visualiser.set_plot_dim(2)
    visualiser.plot()

if(sys.argv[1] == 'plot3d'):
    if(len(sys.argv) > 2):
        if(sys.argv[2] == 'superpunto'):
            visualiser.enable_superpunto()
    visualiser.set_plot_dim(3)
    visualiser.plot()