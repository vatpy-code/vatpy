'''
Description:
Authour(s): Jonathan Petersson
Last updated: 2025-05-16
'''


# -------------- Required Packages
import os
import numpy as np


# -------------- Import UniformBox
from vatpy import const
from vatpy import UniformBox


# -------------- Parameters
# Please provide all parameter values in cgs units

ulength = 100 * const['pc']
umass = const['Msol']
uvel = 1e5

filename = 'uniform_box'
savepath = os.getcwd()
boxsize = 1000 * const['pc']

N = 100000
ambientdensity = 1e-25

pos_star = np.array([250, 250, 500]) * const['pc']
vel_star = np.array([0, 0, 0])
mass_star = 40 * const['Msol']

pos_sink = np.array([750, 750, 500]) * const['pc']
vel_sink = np.array([0, 0, 0])
mass_sink = 1e5 * const['Msol']

relax = True
NumRelax = 1
meshrelaxdir = '/uniform_box/'

# -------------- Run script
print('\n Welcome to this IC generator script (powered by Vatpy):')

UB = UniformBox(ulength=ulength, umass=umass, uvel=uvel, boxsize=boxsize)
UB.generate(N=N, ambientdensity=ambientdensity, pos_star=pos_star,
            vel_star=vel_star, mass_star=mass_star, pos_sink=pos_sink,
            vel_sink=vel_sink, mass_sink=mass_sink, filename=filename,
            savepath=savepath, relax=relax, NumRelax=NumRelax,
            meshrelaxdir=meshrelaxdir)

print('  * Done!\n')

# -------------- End of file
