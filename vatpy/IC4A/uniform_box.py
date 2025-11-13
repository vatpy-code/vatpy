'''
Description:
Authour(s): Jonathan Petersson
Last updated: 2025-05-16
'''


# -------------- Required packages
import os
import numpy as np
import h5py

from .write import write_ic_file
from .check import check_ic_file
from .relax import do_meshrelax

from vatpy import const


# -------------- Function to modify IC file after mesh relax
def modify_after_meshrelax(ic):
    '''
    Description: TODO
    '''

    return None


# -------------- UniformBox
class UniformBox:
    '''
    UniformBox: Class to generate initial conditions of an uniform box, with
                and without the possibility to add a star particle and/or a
                black hole to the simulation domain.
    '''

    def __init__(self, ulength, umass, uvel, boxsize):
        # Simulation set-up:
        self.ulength = ulength
        self.umass = umass
        self.uvel = uvel
        self.utime = ulength/uvel
        self.udens = umass/(ulength**3)
        self.ucoldens = umass/(ulength**2)
        self.uinterg = uvel**2
        self.uangmom = ulength * uvel
        self.boxsize = boxsize

    ##########################################################################
    ##########################################################################
    def uniform_background(self, N, ambientdensity, T=1e4):
        '''
        Description: Generate a background of uniform density.
        '''
        # Coordinates:
        pos = self.boxsize * np.random.rand(N, 3)

        # Mass:
        mass = np.full(N, ambientdensity)

        # Velocities:
        vel = np.zeros((N, 3))

        # Internal Energy:
        mu = 1 + 4*0.1
        interg = np.full(N, (3 * const['kb'] * T) / (2 * mu * const['mp']))

        ic_uniform_background = {
            'Coordinates': pos / self.ulength,
            'Velocities': vel / self.uvel,
            'Masses': mass / self.udens,
            'InternalEnergy': interg / self.uinterg,
        }

        return ic_uniform_background

    def star_particle(self, pos, vel, mass):
        '''
        Description: Generate a star particle at a given location.
        '''
        ic_star = {
            'Coordinates': np.array([[pos[0], pos[1], pos[2]]]) / self.ulength,
            'Velocities': np.array([[vel[0], vel[1], vel[2]]]) / self.uvel,
            'Masses': np.array([mass]) / self.umass
         }

        return ic_star

    def sink_particle(self, pos, vel, mass):
        '''
        Description: Generate a sink particle at a given location.
        '''
        ic_sink = {
            'Coordinates': np.array([[pos[0], pos[1], pos[2]]]) / self.ulength,
            'Velocities': np.array([[vel[0], vel[1], vel[2]]]) / self.uvel,
            'Masses': np.array([mass]) / self.umass
         }

        return ic_sink

    ##########################################################################
    ##########################################################################
    def generate(self, N, ambientdensity, T, pos_star, vel_star, mass_star,
                 pos_sink, vel_sink, mass_sink, filename, savepath,
                 check=False, relax=False, NumRelax=1, meshrelaxdir=None,
                 wait='manual'):
        '''
        Description: Function to generate a box of uniform background density,
                     with a star particle and sink particle on each side of the
                     simulation domain.
        '''

        # Generate ICs for PartType0:
        ic_gas = self.uniform_background(N=N, ambientdensity=ambientdensity,
                                         T=T)

        # Generate ICs for PartType4:
        ic_star = self.star_particle(pos=pos_star, vel=vel_star,
                                     mass=mass_star)

        # Generate ICs for PartType5:
        ic_sink = self.sink_particle(pos=pos_sink, vel=vel_sink,
                                     mass=mass_sink)

        # Generate ParticleIDs:
        ic_gas['ParticleIDs'] = np.arange(0, N) + 1
        ic_star['ParticleIDs'] = np.array([ic_gas['ParticleIDs'][-1] + 1])
        ic_sink['ParticleIDs'] = np.array([ic_star['ParticleIDs'][-1] + 1])

        # Write IC to file:
        parttypes = {'PartType0': ic_gas}
        if mass_star > 0:
            parttypes['PartType4'] = ic_star
        if mass_sink > 0:
            parttypes['PartType5'] = ic_sink

        write_ic_file(filename=filename, savepath=savepath,
                      boxsize=self.boxsize / self.ulength, parttypes=parttypes)

        # Check IC-file:
        if check is True:
            check_ic_file(f'{savepath}/{filename}.hdf5')

        # Mesh relaxation:
        if relax is True:
            jobdir = os.getcwd() + '/meshrelax/' + meshrelaxdir
            modify_args = {}
            do_meshrelax(ic=f'{savepath}/{filename}.hdf5', runs=NumRelax,
                         jobdir=jobdir, modify=modify_after_meshrelax,
                         modify_args=modify_args, filename=filename,
                         savepath=savepath, wait=wait)

        # Final edit(s):
        if mass_star > 0:
            with h5py.File(f'{savepath}/{filename}.hdf5', 'r+') as f:
                # Units:
                umass = f['Parameters'].attrs['UnitMass_in_g']

                # Modify the IC file:
                stellar_masses = np.zeros(50)
                stellar_masses[0] = mass_star / umass

                f['PartType4']['StellarFormationTime'] = np.array([0.]),
                f['PartType4']['StellarMasses'] = np.array([stellar_masses]),
                f['PartType4']['NumberOfSupernovae'] = np.array([1])

        return None

# -------------- End of file
