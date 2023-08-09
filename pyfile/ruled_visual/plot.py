import os
import visual
import sys

visualiser = visual.VISUAL()
visualiser.index = 14


if(sys.argv[1] == 'plot'):
    visualiser.read_rules()
    visualiser.plot()

## Filaments
if(sys.argv[1] == 'plot_fil'):
    visualiser.read_rules()
    if(len(sys.argv) > 2):
        if(sys.argv[2] == 'video'):
            visualiser.video = True
    visualiser.plot_fil()


## Ciliates
# Single sim
if(sys.argv[1] == 'phase'):
    visualiser.read_rules()
    if(len(sys.argv) > 2):
        if(sys.argv[2] == 'video'):
            visualiser.video = True
    visualiser.phase()

if(sys.argv[1] == 'ciliate'):
    visualiser.read_rules()
    if(len(sys.argv) > 2):
        if(sys.argv[2] == 'video'):
            visualiser.video = True
    visualiser.ciliate()

if(sys.argv[1] == 'ciliate_speed'):
    visualiser.read_rules()
    visualiser.ciliate_speed()

if(sys.argv[1] == 'ciliate_traj'):
    visualiser.read_rules()
    visualiser.ciliate_traj()

if(sys.argv[1] == 'timing'):
    visualiser.read_rules()
    visualiser.timing()

if(sys.argv[1] == 'ciliate_forcing'):
    visualiser.read_rules()
    visualiser.ciliate_forcing()

if(sys.argv[1] == 'ciliate_dissipation'):
    visualiser.read_rules()
    visualiser.ciliate_dissipation()    


# Multi sims
if(sys.argv[1] == 'multi_phase'):
    visualiser.read_rules()
    visualiser.multi_phase()

if(sys.argv[1] == 'multi_ciliate'):
    visualiser.read_rules()
    visualiser.multi_ciliate()

if(sys.argv[1] == 'multi_ciliate_traj'):
    visualiser.read_rules()
    visualiser.multi_ciliate_traj()

if(sys.argv[1] == 'multi_ciliate_speed'):
    visualiser.read_rules()
    visualiser.multi_ciliate_speed()

if(sys.argv[1] == 'multi_timing'):
    visualiser.read_rules()
    visualiser.multi_timing()


# Summary plot
if(sys.argv[1] == 'summary_ciliate_speed'):
    visualiser.read_rules()
    visualiser.summary_ciliate_speed()

if(sys.argv[1] == 'summary_timing'):
    visualiser.read_rules()
    visualiser.summary_timing()

if(sys.argv[1] == 'summary_ciliate_dissipation'):
    visualiser.read_rules()
    visualiser.summary_ciliate_dissipation()




#