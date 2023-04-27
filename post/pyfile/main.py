import os
import visual
import sys

visualiser = visual.VISUAL()
# computer = compute_vel.COMPUTEVEL()

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

if(sys.argv[1] == 'plot_seg'):
    visualiser.plot_seg_vel()

if(sys.argv[1] == 'plot_seg_force'):
    visualiser.plot_seg_force()

if(sys.argv[1] == 'plot_pattern'):
    visualiser.enable_superpunto()
    visualiser.plot_pattern()

if(sys.argv[1] == 'plot_fil3d'):
    visualiser.plot_fil3d()

if(sys.argv[1] == 'plot_phase'):
    visualiser.plot_phase_heatmap()

if(sys.argv[1] == 'compute_vel'):
    visualiser.compute_rod_vel()

if(sys.argv[1] == 'multi_vel'):
    visualiser.multi_rod_vel()

if(sys.argv[1] == 'plot_hist'):
    visualiser.plot_hist()

if(sys.argv[1] == 'plot_rod'):
    visualiser.plot_rod_vel()

if(sys.argv[1] == 'plot_rod_multi'):
    visualiser.plot_rod_vel_multi()
