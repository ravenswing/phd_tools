"""
===============================================================================
                                GRAPHICS AND PLOTTING
===============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt


def hills_plot(pdb, funnel_side):
    """ Plot a simple line graph of HILLS file """

    hills = [[], []]
    filename = './monitor/{m}_{f}/{m}.hills'.format(f=funnel_side, m=pdb)
    with open(filename) as f:
        lines = f.readlines()
        data = [l for l in lines if l[0] not in ("@", "#")]
        data = [[float(val) for val in line.split()] for line in data]
        hills[0] = ([l[0] for l in data])
        hills[1] = ([l[5] for l in data])
    plt.figure()
    plt.plot([x/1000 for x in hills[0]], hills[1], label=pdb)
    plt.legend()
    plt.title('{m}  |  {}  |  Hill Heights'.format(funnel_side, m=pdb))
    plt.savefig('monitor/{m}_{f}/{m}-{f}_Heights.png'\
            .format(f=funnel_side, m=pdb), bbox_inches='tight', dpi=300)


def diffusion_plots(pdb, funnel_side, num_cvs):
    """ plots the diffusion plots for X cvs."""

    if num_cvs == 2:
        colvar_data = [[], [], []]
        filename = './monitor/{m}_{f}/{m}.colvar'.format(f=funnel_side, m=pdb)
        with open(filename) as f:
            lines = f.readlines()
            data = [l for l in lines if l[0] not in ("@", "#")]
            data = [[float(val) for val in line.split()] for line in data]
        [x_name, y_name] = lines[0].split()[2:4]
        colvar_data = [[]]*len(data[0])
        for i in np.arange(len(data[0])):
            colvar_data[i] = [l[i] for l in data]
        plt.figure()
        plt.plot([x/1000 for x in colvar_data[0]], colvar_data[1], label=pdb)
        plt.axhline(0.0, color='k', linestyle='--')
        plt.axhline(4.5, color='k', linestyle='--')
        plt.legend()
        plt.xlabel('Simulation Time / ns')
        plt.ylabel(x_name+' / nm')
        plt.ylim(-0.2, 5.0)
        plt.title('{m}  |  {}  |  Projection on Z - pp.proj '\
                .format(funnel_side, m=pdb))
        plt.savefig('monitor/{m}_{f}/{m}-{f}_proj.png'\
                .format(f=funnel_side, m=pdb), bbox_inches='tight', dpi=300)
        plt.figure()
        plt.plot([x/1000 for x in colvar_data[0]], colvar_data[2], label=pdb)
        plt.axhline(0.0, color='k', linestyle='--')
        plt.axhline(2.0, color='k', linestyle='--')
        plt.legend()
        plt.xlabel('Simulation Time / ns')
        plt.ylabel(y_name+' / nm')
        plt.ylim(-0.1, 2.1)
        plt.title('{m}  |  {}  |  Distance from Z - pp.ext'\
                .format(funnel_side, m=pdb))
        plt.savefig('monitor/{m}_{f}/{m}-{f}_ext.png'\
                .format(f=funnel_side, m=pdb), bbox_inches='tight', dpi=300)


def two_cv_contour(pdb, funnel_side, vmax):
    """ Plot a contour plot for 2 CVs"""
    # read in and remove blanks from .fes file
    fes = [[], [], []]
    filename = './monitor/{m}_{f}/{m}.fes'.format(f=funnel_side, m=pdb)
    with open(filename) as f:
        lines = f.readlines()
        data = [l for l in lines if l[0] not in ("@", "#")]
    data = [[float(val) for val in line.split()] for line in data]
    breaks = [i for i, e in enumerate(data) if e == []]
    # remove blank lines
    for index in sorted(breaks, reverse=True):
        del data[index]
    # get number of bins and CV names from header
    nbins = int([l.split()[-1] for l in lines if "nbins_pp.proj" in l][0])
    [x_name, y_name] = lines[0].split()[2:4]
    # organise data into plotable arrays
    split_data = [data[i:i + nbins] for i in range(0, len(data), nbins)]
    z = []
    for block in split_data:
        z.append([l[2] for l in block])
        fes[1].append(block[0][1])
    fes[0] = [l[0] for l in split_data[0]]
    fes[2] = np.asarray(z)
    with open('./monitor/vmax.dat', 'a+') as f:
        f.write('{} {} {}\r\n'.format(pdb, funnel_side, np.amax(fes[2])))

    x, y = np.meshgrid(fes[0], fes[1])
    old_method = False
    if old_method:
        fes[2] = fes[2]-np.amax(fes[2])
        conts = np.arange(-vmax, 0.0, 1.0)
    else:
        #conts = np.append(np.arange(0.0,20.0,2.0), vmax)
        conts = np.arange(0.0, vmax, 1.0)
        print(conts)
    f_x = np.linspace(0.0, 4.5, 1000) # funnel lower & upper walls
    sc = 2.5
    b = 1.5     # funnel beta-cent
    f = 0.15    # funnel wall buffer
    h = 1.7     # funnel wall width
    f_y = h*(1./(1.+np.exp(b*(f_x-sc))))+f

    plt.figure()
    plt.contourf(x, y, fes[2]/4.184, conts, cmap='CMRmap', )
    plt.colorbar(label='Free Energy Surface / kcal/pdb')
    plt.plot(f_x, f_y, 'k')
    plt.xlim(-0.2, 5.0)
    plt.ylim(-0.1, 2.0)
    plt.title('{m}  |  {}  |  Free Energy Surface'.format(funnel_side, m=pdb))
    plt.xlabel(x_name+' / nm')
    plt.ylabel(y_name+' / nm')
    plt.savefig('monitor/{m}_{f}/{m}-{f}_FES.png'\
            .format(f=funnel_side, m=pdb), bbox_inches='tight', dpi=300)
