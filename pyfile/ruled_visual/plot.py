import os
import visual
import sys

visualiser = visual.VISUAL()
visualiser.index = 2

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
    # if(len(sys.argv) > 2):
    #     if(sys.argv[2] == 'video'):
    #         visualiser.video = True
    visualiser.ciliate_speed()

if(sys.argv[1] == 'ciliate_traj'):
    visualiser.read_rules()
    visualiser.ciliate_traj()

if(sys.argv[1] == 'plot_fil'):
    visualiser.read_rules()
    if(len(sys.argv) > 2):
        if(sys.argv[2] == 'video'):
            visualiser.video = True
    visualiser.plot_fil()