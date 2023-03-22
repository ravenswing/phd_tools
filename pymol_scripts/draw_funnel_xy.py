#!/usr/bin/python

# Script for drawing a funnel in PyMol.
#
# USAGE:
#
# draw_funnel selection [, dimension [, s_cent [, beta_cent [, wall_width,
#     [, wall_buffer [, lower_wall [, upper_wall [, dim_step [, angle_sample ]]]]]]]]]
#
# ARGUMENTS:
#
# selection = string, a single atom or pseudoatom that will be the origin of the funnel
# dimension = integer, a dimension in which the funnel will be drawn; x=0, y=1, z=2
# s_cent = float, inflection point of the sigmoid function
# beta_cent = float, steepness of the sigmoid function
# wall_width = float, the radius of the funnel at the widest point (excluding the buffer)
# wall_buffer = float, the radius of the funnel at the narrowest point
# lower_wall = float, the starting point of the funnel from the origin along the selected dimension
# upper_wall = float, the ending point of the funnel from the origin along the selected dimension
# dim_step = float, determines how often the funnel points are drawn along the selected dimension
# angle_sample = float, sets the angle between adjacent points at each funnel ring

# by Antonija (based on Giulio's script)

import numpy as np
from pymol import cmd

def draw_funnel(selection, dimension=2, s_cent=25, beta_cent=0.3, wall_width=18.5, wall_buffer=1.5,
                lower_wall=0, upper_wall=45, dim_step=2.5, angle_sample=18):
    dimension = int(dimension)
    s_cent = float(s_cent)
    beta_cent = float(beta_cent)
    wall_width = float(wall_width)
    wall_buffer = float(wall_buffer)
    lower_wall = float(lower_wall)
    upper_wall = float(upper_wall)
    dim_step = float(dim_step)
    angle_sample = float(angle_sample)
    # Get coords of the origin
    origin = cmd.get_coords(selection, 1)[0]
    print 'Origin:', origin


    dimList = [0,1,2]
    del dimList[dimension]

    # Iterate in the selected dimension
    for dim in np.arange(origin[dimension] + lower_wall, origin[dimension] + upper_wall, dim_step):
        # Iterate around a circle, with its radius defined by the sigmoid function
        radius = (wall_width / (1 + np.exp(beta_cent * (dim - (origin[dimension] + s_cent))))) + wall_buffer
        print dim, radius
        for angle in np.arange(-np.pi, np.pi, 2 * np.pi / angle_sample):
            # Generate pseudoatom. Assign b-factor = dim for coloring purposes
            pos = [0,0,0]
            pos[dimension] = dim
            pos[dimList[0]] = origin[dimList[0]] + radius * np.sin(angle)
            pos[dimList[1]] = origin[dimList[1]] + radius * np.cos(angle)

            cmd.pseudoatom('funnel', pos=pos, b=dim)

    cmd.color('orange', selection='funnel')
    cmd.show_as('nonbonded', 'funnel')
    cmd.show('spheres', 'not funnel and ' + selection)

cmd.extend("draw_funnel", draw_funnel)