import os
import visual
import sys

visualiser = visual.VISUAL()

if(sys.argv[1] == 'plot2d'):
    visualiser.plot2d()

if(sys.argv[1] == 'plot3d'):
    visualiser.plot3d()