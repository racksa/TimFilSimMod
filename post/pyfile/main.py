import os
import visual
import compute_vel
import sys

visualiser = visual.VISUAL()
computer = compute_vel.COMPUTEVEL()

if(sys.argv[1] == 'plot2d'):
    visualiser.set_plot_dim(2)
    visualiser.plot()

if(sys.argv[1] == 'plot3d'):
    if(len(sys.argv) > 2):
        if(sys.argv[2] == 'superpunto'):
            visualiser.enable_superpunto()
        if(sys.argv[2] == 'fcm'):
            visualiser.enable_fcm()
    visualiser.set_plot_dim(3)
    visualiser.plot()

if(sys.argv[1] == 'compute_vel'):
    computer.compute()

if(sys.argv[1] == 'plot_hist'):
    computer.plot_hist()

if(sys.argv[1] == 'plot_seg'):
    computer.plot_seg_height()
