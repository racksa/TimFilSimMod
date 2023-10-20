import os
import visual
import sys

# initialise the visuliser
visualiser = visual.VISUAL()

# lists of available options
list_of_single_options = ['plot', 'plot_fil', 'phase', 'eckert', 'ciliate', 'ciliate_speed', 'ciliate_traj', 
                   'timing', 'ciliate_forcing', 'ciliate_dissipation',
                   'ciliate_svd', 'ciliate_dmd']
list_of_multi_options = ['multi_phase', 'multi_ciliate', 'multi_ciliate_traj',
                         'multi_ciliate_speed', 'multi_timing', 'multi_ciliate_svd']
list_of_summary_options = ['summary_ciliate_speed', 'summary_timing',
                           'summary_ciliate_dissipation']
list_of_all_options = list_of_single_options + list_of_multi_options + list_of_summary_options


# execute the plotting function
if(sys.argv[1] in list_of_all_options):
    visualiser.read_rules()
    if('interpolate' in sys.argv):
            visualiser.interpolate = True

    if(sys.argv[1] in list_of_single_options):
        if(len(sys.argv) > 2):
            if(sys.argv[2].isdigit()):
                visualiser.index = int(sys.argv[2])
        if('video' in sys.argv):
            visualiser.video = True
        

    if hasattr(visualiser, sys.argv[1]):
        method_to_call = getattr(visualiser, sys.argv[1])
        if callable(method_to_call):
            method_to_call()
        else:
            print(f"'{sys.argv[1]}' is not a callable method.")











#