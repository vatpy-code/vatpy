'''
Description: TODO
Authour(s): Jonathan Petersson
Last updated: 2025-06-04
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
def modify_after_meshrelax_ucloud(ic, centre, R, rho):
    '''
    Description: TODO
    '''

    with h5py.File(ic, 'r+') as f:
        # Units:
        ulength = f['Header'].attrs['UnitLength_in_cm']
        umass = f['Header'].attrs['UnitMass_in_g']
        udens = umass/(ulength**3)

        # Modify the ICs:
        pos = f['PartType0']['Coordinates'] * ulength
        mass = f['PartType0']['Masses'] * umass
        mask = np.linalg.norm(pos - centre, axis=1) < R
        mass[mask] = rho
        mass[~mask] = 1e-25
        f['PartType0']['Masses'][:] = mass / udens

    return None


def modify_after_meshrelax_r2cloud(ic, centre, R, rho):
    '''
    Description: TODO
    '''

    with h5py.File(ic, 'r+') as f:
        # Units:
        ulength = f['Header'].attrs['UnitLength_in_cm']
        umass = f['Header'].attrs['UnitMass_in_g']
        udens = umass/(ulength**3)

        # Modify the ICs:
        pos = f['PartType0']['Coordinates'] * ulength
        mass = f['PartType0']['Masses'] * umass
        radius = np.linalg.norm(pos - centre, axis=1)
        mask1 = radius < R
        mask2 = radius < 1.0 * const['pc']
        mass[mask1] = rho * (1.0 * const['pc'] / radius[mask1])**2
        mass[mask2] = rho
        mass[~mask1] = 1e-28
        f['PartType0']['Masses'][:] = mass / udens

    return 0


class CollapsingCloud:
    '''
    CloudCollapse: initial conditions for collapsing cloud(s)
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
    def ucloud(self, N, centre, kick, R, rho, T, P, rot):
        '''
        Description: TODO
        '''

        # Coordinates:
        pos = self.boxSize * np.random.rand(N, 3)

        # Mass:
        mass = np.full(N, 1e-25)
        mask = np.linalg.norm(pos - centre, axis=1) < R
        mass[mask] = rho

        # Velocities:
        vel = self.velocity(pos - centre, R, kick, P, rot)

        # Internal Energy:
        mu = 1 + 4*0.1
        interg = np.full(N, (3 * self.kb * 1e4) / (2 * mu * self.mp))
        interg[mask] = (3 * self.kb * T) / (2 * mu * self.mp)

        ic_cloud = {
            'Coordinates': pos / self.ulength,
            'Velocities': vel / self.uvel,
            'Masses': mass / self.udens,
            'InternalEnergy': interg / self.uinterg,
        }

        return ic_cloud

    def velocity(self, pos, R, kick, P, rot):
        '''
        Description: TODO
        '''

        # Velocity:
        if P > 0:
            # Rotation axis:
            rot = rot / np.linalg.norm(rot)

            # Coordinate projection onto the plane of the rotation axis:
            pos1 = np.dot(pos, rot)[:, np.newaxis] * rot
            pos2 = pos - pos1

            # Velocity vectors based on the angular velocity:
            vec = np.cross(rot, pos2)
            vec = vec / np.linalg.norm(vec)
            w = 2 * np.pi / P
            v = w * np.linalg.norm(pos2, axis=1)
            vel = v[:, np.newaxis] * vec
            mask = np.linalg.norm(pos, axis=1) < R
            vel[~mask] = 0
        else:
            vel = np.zeros(np.shape(pos))

        vel += kick

        return vel

    ##########################################################################
    ##########################################################################
    def r2cloud(self, N, centre, kick, R, rho, T):
        '''
        Description: TODO
        '''

        # Coordinates:
        N_cloud = int(0.9 * N)
        N_medium = int(0.1 * N)
        pos_cloud = R * (2 * np.random.rand(N_cloud, 3) - 1) + centre
        pos_medium = self.boxSize * np.random.rand(N_medium, 3)
        pos = np.append(pos_cloud, pos_medium, axis=0)

        # Mass:
        dens = np.full(N, 1e-28)
        radius = np.linalg.norm(pos - centre, axis=1)
        mask1 = radius < R
        mask2 = radius < 1.0 * self.pc
        dens[mask1] = rho * (1.0 * self.pc / radius[mask1])**2
        dens[mask2] = rho

        # Velocities:
        vel = np.zeros(np.shape(pos))
        vel += kick

        # Internal Energy:
        mu = 1 + 4*0.1
        interg = np.full(N, (3 * self.kb * 1e4) / (2 * mu * self.mp))
        interg[mask1] = (3 * self.kb * T) / (2 * mu * self.mp)

        ic_cloud = {
            'Coordinates': pos / self.ulength,
            'Velocities': vel / self.uvel,
            'Masses': dens / self.udens,
            'InternalEnergy': interg / self.uinterg,
        }

        return ic_cloud

    ##########################################################################
    ##########################################################################
    def u2clouds(self, N, centre1, centre2, kick1, kick2, R1, R2, rho1, rho2,
                 T1, T2):
        '''
        Description: TODO
        '''

        # Coordinates:
        pos = self.boxSize * np.random.rand(N, 3)

        # Mass:
        mass = np.full(N, 1e-24)
        mask1 = np.linalg.norm(pos - centre1, axis=1) < R1
        mask2 = np.linalg.norm(pos - centre2, axis=1) < R2
        mass[mask1] = rho1
        mass[mask2] = rho2

        # Velocities:
        vel = np.zeros((N, 3))
        vel[mask1] += kick1
        vel[mask2] += kick2

        # Internal Energy:
        mu = 1 + 4*0.1
        interg = np.full(N, (3 * self.kb * 1e4) / (2 * mu * self.mp))
        interg[mask1] = (3 * self.kb * T1) / (2 * mu * self.mp)
        interg[mask2] = (3 * self.kb * T2) / (2 * mu * self.mp)

        ic_clouds = {
            'Coordinates': pos / self.ulength,
            'Velocities': vel / self.uvel,
            'Masses': mass / self.udens,
            'InternalEnergy': interg / self.uinterg,
        }

        return ic_clouds

    ##########################################################################
    ##########################################################################
    def sink_particle(self, x, y, z, vx, vy, vz, m):
        '''
        Description: TODO
        '''

        ic_sink = {
            'Coordinates': np.array([[x, y, z]]) / self.ulength,
            'Velocities': np.array([[vx, vy, vz]]) / self.uvel,
            'Masses': np.array([m]) / self.umass,
         }

        return ic_sink

    ##########################################################################
    ##########################################################################
    def generate_ucloud(self, N, pos_cloud, vel_cloud, R, rho, T, P, rot,
                        pos_sink, vel_sink, mass_sink, filename, savepath,
                        check=False, relax=False, NumRelax=1, wait='manual'):
        '''
        Description: TODO
        '''

        # Generate ICs for PartType0:
        ic_gas = self.ucloud(N=N, centre=pos_cloud, kick=vel_cloud, R=R,
                             rho=rho, T=T, P=P, rot=rot)

        # Generate ICs for PartType5:
        ic_sink = self.sink_particle(x=pos_sink[0], y=pos_sink[1],
                                     z=pos_sink[2], vx=vel_sink[0],
                                     vy=vel_sink[1], vz=vel_sink[2],
                                     m=mass_sink)

        # Generate ParticleIDs:
        ic_gas['ParticleIDs'] = np.arange(0, N) + 1
        ic_sink['ParticleIDs'] = np.array([ic_gas['ParticleIDs'][-1] + 1])

        # Write IC-file:
        write_ic_file(filename=filename, savepath=savepath,
                      boxsize=self.boxsize / self.ulength,
                      partTypes={'PartType0': ic_gas, 'PartType5': ic_sink})

        # Check IC-file:
        if check is True:
            check_ic_file(f'{savepath}/{filename}.hdf5')

        # Mesh relaxation:
        if relax is True:
            jobdir = os.getcwd() + '/meshrelax/cloudcollapse/'
            modify = modify_after_meshrelax_ucloud
            modify_args = {'centre': pos_cloud, 'R': R, 'rho': rho}
            do_meshrelax(ic=f'{savepath}/{filename}.hdf5', runs=NumRelax,
                         jobdir=jobdir, modify=modify, modify_args=modify_args,
                         filename=filename, savepath=savepath, wait=wait)

            # Final update of the velocity:
            with h5py.File(f'{savepath}/{filename}.hdf5', 'r+') as f:
                ulength = f['Header'].attrs['UnitLength_in_cm']
                uvel = f['Header'].attrs['UnitVelocity_in_cm_per_s']
                pos = f['PartType0']['Coordinates'] * ulength
                pos = pos - pos_cloud
                vel = self.velocity(pos=pos, R=R, kick=vel_cloud, P=P, rot=rot)
                f['PartType0']['Velocities'][:] = vel / uvel

        return None

    ##########################################################################
    ##########################################################################
    def generate_r2cloud(self, N, pos_cloud, vel_cloud, R, rho, T, pos_sink,
                         vel_sink, mass_sink, filename, savepath, check=False,
                         relax=False, NumRelax=1, wait='manual'):
        '''
        Description: TODO
        '''

        # Generate ICs for PartType0:
        ic_gas = self.r2cloud(N=N, centre=pos_cloud, kick=vel_cloud, R=R,
                              rho=rho, T=T)

        # Generate ICs for PartType5:
        ic_sink = self.sink_particle(x=pos_sink[0], y=pos_sink[1],
                                     z=pos_sink[2], vx=vel_sink[0],
                                     vy=vel_sink[1], vz=vel_sink[2],
                                     m=mass_sink)

        # Generate ParticleIDs:
        ic_gas['ParticleIDs'] = np.arange(0, N) + 1
        ic_sink['ParticleIDs'] = np.array([ic_gas['ParticleIDs'][-1] + 1])

        # Write IC-file:
        write_ic_file(filename=filename, savepath=savepath,
                      boxsize=self.boxsize / self.ulength,
                      partTypes={'PartType0': ic_gas, 'PartType5': ic_sink})

        # Check IC-file:
        if check is True:
            check_ic_file(f'{savepath}/{filename}.hdf5')

        # Mesh relaxation:
        if relax is True:
            jobdir = os.getcwd() + '/meshrelax/cloudcollapse/'
            modify = modify_after_meshrelax_r2cloud
            modify_args = {'centre': pos_cloud, 'R': R, 'rho': rho}
            do_meshrelax(ic=f'{savepath}/{filename}.hdf5', runs=NumRelax,
                         jobdir=jobdir, modify=modify, modify_args=modify_args,
                         filename=filename, savepath=savepath, wait=wait)

            # Final update of the velocity:
            with h5py.File(f'{savepath}/{filename}.hdf5', 'r+') as f:
                vel = f['PartType0']['Velocities'][:]
                vel = np.zeros(np.shape(vel))
                f['PartType0']['Velocities'][:] = vel

        return None

    ##########################################################################
    ##########################################################################
    def generate_u2clouds(self, N, pos_cloud1, pos_cloud2, vel_cloud1,
                          vel_cloud2, R1, R2, rho1, rho2, T1, T2, pos_sink,
                          vel_sink, mass_sink, filename, savepath, check=False,
                          relax=False, NumRelax=1, wait='manual'):
        '''
        Description: TODO
        '''

        # Generate ICs for PartType0:
        ic_gas = self.twoclouds(N=N, centre1=pos_cloud1, centre2=pos_cloud2,
                                kick1=vel_cloud1, kick2=vel_cloud2, R1=R1,
                                R2=R2, rho1=rho1, rho2=rho2, T1=T1, T2=T2)

        # Generate ICs for PartType5:
        ic_sink = self.sink_particle(x=pos_sink[0], y=pos_sink[1],
                                     z=pos_sink[2], vx=vel_sink[0],
                                     vy=vel_sink[1], vz=vel_sink[2],
                                     m=mass_sink)

        # Generate ParticleIDs:
        ic_gas['ParticleIDs'] = np.arange(0, N) + 1
        ic_sink['ParticleIDs'] = np.array([ic_gas['ParticleIDs'][-1] + 1])

        # Write IC-file:
        write_ic_file(filename=filename, savepath=savepath,
                      boxsize=self.boxsize / self.ulength,
                      partTypes={'PartType0': ic_gas, 'PartType5': ic_sink})

        # Check IC-file:
        if check is True:
            check_ic_file(f'{savepath}/{filename}.hdf5')

        # Mesh relaxation:
        if relax is True:
            jobdir = os.getcwd() + '/meshrelax/cloudcollapse/'
            modify = None
            modify_args = None
            do_meshrelax(ic=f'{savepath}/{filename}.hdf5', runs=NumRelax,
                         jobdir=jobdir, modify=modify, modify_args=modify_args,
                         filename=filename, savepath=savepath, wait=wait)

        return None

# -------------- End of file
