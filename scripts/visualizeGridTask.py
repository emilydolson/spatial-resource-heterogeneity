#!/usr/bin/python

import sys
from avidaSpatialTools import *

if len(sys.argv) < 2:
    print "USAGE: ./visualizeGridTask.py <grid_task file>"

make_species_grid(sys.argv[1])
#plt.imshow(make_n_tasks_grid(sys.argv[1]), interpolation="none")
plt.show()
