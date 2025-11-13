# Description: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-05-16


# -------------- Required Packages
import os
import sys
import yaml
import numpy as np


# -------------- Import UniformBox
from vatpy import const
from vatpy import UniformBox


# -------------- Parameters
print('\n Welcome to this IC generator script (powered by Vatpy):')

# Parameter file & group:
paramfile = sys.argv[1]
paramgroup = sys.argv[2]

print(f'  * Parameters: {paramgroup} in {paramfile}')

with open(paramfile, 'r') as file:
    f = yaml.safe_load(file)
    param = f[paramgroup]

    # Please provide all parameter values in cgs units
    ulength = 100 * const['pc']
    umass = const['Msol']
    uvel = 1e5

    filename = 'uniform_box'
    savepath = os.getcwd()
    boxsize = param['boxsize'] * const['pc']

    N = param['N']
    ambientdensity = float(param['ambientdensity'])
    T = float(param['gastemperature'])

    pos_star = np.array(param['pos_star']) * const['pc']
    vel_star = np.array([0, 0, 0])
    mass_star = param['mass_star'] * const['Msol']

    pos_sink = np.array(param['pos_sink']) * const['pc']
    vel_sink = np.array([0, 0, 0])
    mass_sink = float(param['mass_sink']) * const['Msol']

    relax = param['relax']
    NumRelax = param['NumRelax']
    meshrelaxdir = '/uniform_box/'

# -------------- Run script
UB = UniformBox(ulength=ulength, umass=umass, uvel=uvel, boxsize=boxsize)
UB.generate(N=N, ambientdensity=ambientdensity, T=T, pos_star=pos_star,
            vel_star=vel_star, mass_star=mass_star, pos_sink=pos_sink,
            vel_sink=vel_sink, mass_sink=mass_sink, filename=filename,
            savepath=savepath, relax=relax, NumRelax=NumRelax,
            meshrelaxdir=meshrelaxdir)

print('  * Done!\n')

# -------------- End of file
