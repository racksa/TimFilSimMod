import os
import visual
import sys

visualiser = visual.VISUAL()
visualiser.index = 0

if(sys.argv[1] == 'plot'):
    visualiser.read_rules()
    visualiser.plot()

if(sys.argv[1] == 'plot_phase'):
    visualiser.read_rules()
    if(len(sys.argv) > 2):
        if(sys.argv[2] == 'video'):
            visualiser.video = True
    visualiser.plot_phase()

if(sys.argv[1] == 'plot_ciliate'):
    visualiser.read_rules()
    if(len(sys.argv) > 2):
        if(sys.argv[2] == 'video'):
            visualiser.video = True
    visualiser.plot_ciliate()

if(sys.argv[1] == 'ciliate_speed'):
    visualiser.read_rules()
    visualiser.ciliate_speed()

if(sys.argv[1] == 'ciliate_traj'):
    visualiser.read_rules()
    visualiser.ciliate_traj()

if(sys.argv[1] == 'timing'):
    visualiser.read_rules()
    visualiser.timing()



if(sys.argv[1] == 'plot_fil'):
    visualiser.read_rules()
    if(len(sys.argv) > 2):
        if(sys.argv[2] == 'video'):
            visualiser.video = True
    visualiser.plot_fil()



if(sys.argv[1] == 'plot_multi_phase'):
    visualiser.read_rules()
    visualiser.plot_multi_phase()

if(sys.argv[1] == 'plot_multi_ciliate'):
    visualiser.read_rules()
    visualiser.plot_multi_ciliate()

if(sys.argv[1] == 'ciliate_multi_traj'):
    visualiser.read_rules()
    visualiser.ciliate_multi_traj()