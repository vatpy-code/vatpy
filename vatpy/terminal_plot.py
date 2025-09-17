# Description: File containing the TerminalPlot class.
# Authour(s): Jonathan Petersson
# Last updated: 2025-08-27


# -------------- Required packages
import os
import numpy as np
import cmasher as cmr
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.spatial.transform import Rotation
from scipy.ndimage import gaussian_filter
from labellines import labelLines

from .read import read_hdf5
from .constants import const
from .get_gas_property import number_density, temperature
from .interpolation import interpolate_to_2d, interpolate_to_2d_kdtree
from .get_data import get_black_hole_data

import configv


# -------------- TerminalPlot
class TerminalPlot:
    '''Class to make simple, but informative, visual plots of Arepo snapshots,
    directly in the terminal (or notebook). Most importantly, it contains
    funtions to generate column density maps of the gas surface density, as
    well column density maps of various chemical species, such as HI, HII, and
    H2. The class also contains functions to inspect the surface density of
    dark matter and stellar components. For more details, see each function's
    individual description.

    Args:
        file (str): File to analyse
        style (str): Matplotlib style option
        path (str): Path to save file(s) at
        name (str): Name to save file(s) as
        format (str): Format to save file(s) as (e.g. png/jpg/pdf/tiff)
        vmin (float): Global vmin value
        vmax (float): Global vmax value
        xlim (tuple): Global x-axis min and max values
        ylim (tuple): Global y-axis min and max values
        ulengthselect (str): Unit length (see configv.py for more details)
        show (bool): If True, try to display the generated figure(s)
    '''
    def __init__(self, file, style=configv.mplstyle, path='./vplots',
                 name=None, format='png', vmin=None, vmax=None, xlim=None,
                 ylim=None, ulengthselect=configv.unit_for_length, show=True):

        print(f'  * Setting up an environment to analyse {file}')

        # Variables:
        self.file = file
        self.path = path
        self.name = name
        self.format = format
        self.vmin, self.vmax = vmin, vmax
        self.xlim, self.ylim = xlim, ylim
        self.show = show

        # Mpl style:
        self.style = style
        plt.style.use(f'{configv.homedir}/vatpy/mpl/{self.style}.mplstyle')

        # Unit selection:
        self.ulengthselect = ulengthselect
        if self.ulengthselect == 'kpc':
            self.ulength = const['kpc']
        elif self.ulengthselect == 'pc':
            self.ulength = const['pc']
        else:
            self.ulength = 1

    ##########################################################################
    ##########################################################################
    def get_ranges(self, boxsize, box, xrange, yrange, zrange, bhfocus):
        '''
        Description: Function to get x, y, and z ranges.
        '''
        # Determine the coordinate ranges:
        if not box:
            if not xrange:
                if bhfocus:
                    xrange = (-boxsize/2, boxsize/2)
                else:
                    xrange = (0, boxsize)
            if not yrange:
                if bhfocus:
                    yrange = (-boxsize/2, boxsize/2)
                else:
                    yrange = (0, boxsize)
            if not zrange:
                if bhfocus:
                    zrange = (-boxsize/2, boxsize/2)
                else:
                    zrange = (0, boxsize)
        else:
            xrange = (box[0], box[1])
            yrange = (box[0], box[1])
            zrange = (box[0], box[1])

        return xrange, yrange, zrange

    def do_rotation(self, boxsize, axis, rotate, pos, bhfocus):
        '''
        Description: Function to rotate the position of particles.
        '''
        # If applicable, rotate the position of particles:
        if rotate != 0:
            rotation = Rotation.from_euler(axis, rotate, degrees=True)
            if bhfocus:
                pos = rotation.apply(pos)
            else:
                pos = rotation.apply(pos - boxsize/2)
                pos += boxsize/2

        return pos

    def save(self, fig, funcname):
        '''
        Description: Function to save figures.
        '''
        # Name selection:
        if self.name:
            savename = self.name
        else:
            savename = funcname

        # Figure name:
        filestr = self.file.split('.')
        figname = f'{savename}_{filestr[0]}'
        if self.name:
            figname = self.name + '_' + filestr[0][-3:]
        figname += f'.{self.format}'

        # Check if Vatpy plot directory already exists:
        if os.path.isdir(self.path):
            print('  * Path to save figure at detected')
        else:
            print('  * Path to save figure at NOT detected')
            print('  * Creating a \'vplots\' directory')
            os.mkdir(f'{os.getcwd()}/vplots/')

        # Save figure:
        fig.savefig(f'{self.path}/{figname}')
        print('  * Figure saved as')
        print(f'    - name: {figname}')
        print(f'    - at: {self.path}')

        return None

    def display(self):
        '''
        Description: TODO
        '''

        if self.show is not True:
            print('  * Display of figure is turned OFF')
            plt.close()
        else:
            print('  * Display of figure is now running')
            plt.show()

        return None

    ##########################################################################
    ##########################################################################
    def info(self):
        '''Function to provide some general information about the given
        snapshot, such as the physical time, size of the simulation domain,
        number of particles, etc.
        '''
        # Read the data:
        h, iu = read_hdf5(file=self.file)

        time = h['Header'].attrs['Time'] * iu['utime'] / const['Myr']
        boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / const['kpc']
        numpart = h['Header'].attrs['NumPart_ThisFile']

        print('  * Snapshot information')
        print(f'  | Time    : {round(time, 2)} Myr')
        print(f'  | BoxSize : {round(boxsize, 2)} kpc')
        print('  |')
        print('  | Number of particles')
        print(f'  | PartType0 (gas)      : {numpart[0]}')
        print(f'  | PartType1 (halo)     : {numpart[1]}')
        print(f'  | PartType2 (disk)     : {numpart[2]}')
        print(f'  | PartType3 (bulge)    : {numpart[3]}')
        print(f'  | PartType4 (stars)    : {numpart[4]}')
        print(f'  | PartType5 (bndry)    : {numpart[5]}')
        print(f'  | Total Particle Count : {np.sum(numpart)}')
        print('  | ')

        mt = h['Header'].attrs['MassTable']
        print('  | Table of Particle Masses [i.u.]')
        print(f'  | {mt[0]}, {mt[1]}, {mt[2]}, {mt[3]}, {mt[4]}, {mt[5]}')

        if numpart[5] == 1:
            print('  |')
            print('  * A central BH detected')
            Pbh = h['PartType5']['Coordinates'][0]
            print(f'  | Coordinates [i.u.] : ({Pbh[0]}, {Pbh[1]}, {Pbh[2]})')
            Vbh = h['PartType5']['Velocities'][0]
            print(f'  | Velocities [i.u.]  : ({Vbh[0]}, {Vbh[1]}, {Vbh[2]})')
            Mbh = h['PartType5']['Masses'][0] * iu['umass'] / const['Msol']
            print(f'  | Mass [i.u.]        : {Mbh}')
            IDbh = h['PartType5']['ParticleIDs'][0]
            print(f'  | Particle ID        : {IDbh}')

        print('  |')

        return None

    ##########################################################################
    ##########################################################################
    def density(self, axis='z', rotate=0, quantity='mass', bins=100,
                interpolation='kdtree', bhfocus=False, funcname='dens',
                xrange=None, yrange=None, zrange=None, box=None, cut=None):
        '''Function to generate a column density map of the gas surface
        density, with the possibility to also show the column density of
        various chemical species, such as HI, HII, and H2. This is achieved by
        first interpolating the selected gas quantity onto a grid, and later
        doing a sum along the line-of-sight.

        Args:
            axis (str): Rotation axis
            rotate (float): Amount of rotation
            quantity (str): Gas quantity to visualise
                            (options: mass/n/HI/HII/H2/CO/He/e)
        '''
        # Read the data:
        print(f'  * Reading data of {self.file}')
        h, iu = read_hdf5(file=self.file)
        pos = h['PartType0']['Coordinates'] * iu['ulength'] / self.ulength
        dens = h['PartType0']['Density'] * iu['udens']
        time = h['Header'].attrs['Time'] * iu['utime'] / const['Myr']
        boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / self.ulength

        print('  * Generating a gas surface density map')

        # Centre the data on the black hole:
        if bhfocus:
            bh = (h['PartType5']['Coordinates'][0] * iu['ulength']
                  / self.ulength)
            pos -= bh

        # Coordinate ranges:
        xrange, yrange, zrange = self.get_ranges(boxsize=boxsize, box=box,
                                                 xrange=xrange, yrange=yrange,
                                                 zrange=zrange,
                                                 bhfocus=bhfocus)

        # Selection of gas quantity:
        if (quantity != 'mass'):
            num = number_density(h, iu)
            dens = num[quantity]

        # Rotation of particle positions:
        pos = self.do_rotation(boxsize=boxsize, axis=axis, rotate=rotate,
                               pos=pos, bhfocus=bhfocus)

        # Interpolation:
        if interpolation == 'kdtree':
            interpDens = interpolate_to_2d_kdtree(pos=pos, unit=self.ulength,
                                                  values=dens, bins=bins,
                                                  xrange=xrange, yrange=yrange,
                                                  zrange=zrange, cut=cut)
        else:
            interpDens = interpolate_to_2d(pos=pos, unit=self.ulength,
                                           values=dens, bins=bins,
                                           xrange=xrange, yrange=yrange,
                                           zrange=zrange, cut=cut)

        # Plot:
        fig, ax = plt.subplots(figsize=(8, 6.4))
        fig.subplots_adjust(left=0.18, right=0.82, bottom=0.14, top=0.94,
                            wspace=0, hspace=0)

        im = ax.imshow(np.log10(interpDens), vmin=self.vmin, vmax=self.vmax,
                       extent=(xrange[0], xrange[1], yrange[0], yrange[1]),
                       origin='lower', cmap=configv.cmap['gas'])
        if bhfocus:
            ax.scatter(0, 0, s=40, c='k')
        ax.text(0.95, 0.05, f'{time:.2f} Myr', color='k',
                ha='right', va='bottom', transform=ax.transAxes,
                bbox={'facecolor': 'white', 'edgecolor': 'none',
                      'boxstyle': 'round', 'alpha': 0.5})
        ax.set_aspect('equal')
        ax.set_xlabel(f'$x$ [{self.ulengthselect}]')
        ax.set_ylabel(f'$y$ [{self.ulengthselect}]')
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)

        # Colorbar:
        if not cut:
            dim = 2
        else:
            dim = 3
        cbar_label = {
            'HII': (r'$\log_{10}(\Sigma_{\mathrm{HII}}$'
                    r' $[\mathrm{cm}^{-%d}])$' % dim),
            'H2': (r'$\log_{10}(\Sigma_{\mathrm{H}_2}$'
                   r' $[\mathrm{cm}^{-%d}])$' % dim),
            'HI': (r'$\log_{10}(\Sigma_{\mathrm{HI}}$'
                   r' $[\mathrm{cm}^{-%d}])$' % dim),
            'CO': (r'$\log_{10}(\Sigma_{\mathrm{CO}}$'
                   r' $[\mathrm{cm}^{-%d}])$' % dim),
            'He': (r'$\log_{10}(\Sigma_{\mathrm{He}}$'
                   r' $[\mathrm{cm}^{-%d}])$' % dim),
            'e': (r'$\log_{10}(\Sigma_{e^{-}}$'
                  r' $[\mathrm{cm}^{-%d}])$' % dim),
            'n': (r'$\log_{10}(\Sigma_\mathrm{Gas}$'
                  r' $[\mathrm{cm}^{-%d}])$' % dim),
            'mass': (r'$\log_{10}(\Sigma_\mathrm{Gas}$'
                     r' $[\mathrm{g} \ \mathrm{cm}^{-%d}])$' % dim)
        }

        div = make_axes_locatable(ax)
        cax = div.append_axes('right', size='5%', pad=0)
        fig.colorbar(im, cax=cax, label=cbar_label[quantity])

        # Save:
        print('  * Figure generated successfully!')
        self.save(fig=fig, funcname=funcname)

        # Display figure:
        self.display()

        return None

    ##########################################################################
    ##########################################################################
    def temperature(self, axis='z', rotate=0, bins=100, interpolation='kdtree',
                    bhfocus=False, funcname='temp', xrange=None, yrange=None,
                    zrange=None, box=None, cut=None):
        '''
        Description: TODO
        '''
        # Read the data:
        print(f'  * Reading data of {self.file}')
        h, iu = read_hdf5(file=self.file)
        pos = h['PartType0']['Coordinates'] * iu['ulength'] / self.ulength
        dens = h['PartType0']['Density'] * iu['udens']
        time = h['Header'].attrs['Time'] * iu['utime'] / const['Myr']
        boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / self.ulength
        temp = temperature(h, iu)

        print('  * Generating a gas temperature (density-weighted) map')

        # Centre the data on the black hole:
        if bhfocus:
            bh = (h['PartType5']['Coordinates'][0] * iu['ulength']
                  / self.ulength)
            pos -= bh

        # Coordinate ranges:
        xrange, yrange, zrange = self.get_ranges(boxsize=boxsize, box=box,
                                                 xrange=xrange, yrange=yrange,
                                                 zrange=zrange,
                                                 bhfocus=bhfocus)

        # Rotation of particle positions:
        pos = self.do_rotation(boxsize=boxsize, axis=axis, rotate=rotate,
                               pos=pos, bhfocus=bhfocus)

        # Interpolation:
        if interpolation == 'kdtree':
            interpTemp = interpolate_to_2d_kdtree(pos=pos, unit=self.ulength,
                                                  values=temp, bins=bins,
                                                  xrange=xrange, yrange=yrange,
                                                  zrange=zrange, cut=cut,
                                                  weights=dens)
        else:
            interpTemp = interpolate_to_2d(pos=pos, unit=self.ulength,
                                           values=temp, bins=bins,
                                           xrange=xrange, yrange=yrange,
                                           zrange=zrange, cut=cut,
                                           weights=dens)

        # Plot:
        fig, ax = plt.subplots(figsize=(8, 6.4))
        fig.subplots_adjust(left=0.18, right=0.82, bottom=0.14, top=0.94,
                            wspace=0, hspace=0)

        im = ax.imshow(np.log10(interpTemp), vmin=self.vmin, vmax=self.vmax,
                       extent=(xrange[0], xrange[1], yrange[0], yrange[1]),
                       origin='lower', cmap='afmhot')
        ax.text(0.95, 0.05, f'{time:.2f} Myr', color='k',
                ha='right', va='bottom', transform=ax.transAxes,
                bbox={'facecolor': 'white', 'edgecolor': 'none',
                      'boxstyle': 'round', 'alpha': 0.5})
        ax.set_aspect('equal')
        ax.set_xlabel(f'$x$ [{self.ulengthselect}]')
        ax.set_ylabel(f'$y$ [{self.ulengthselect}]')
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)

        div = make_axes_locatable(ax)
        cax = div.append_axes('right', size='5%', pad=0)
        fig.colorbar(im, cax=cax, label=r'$\log_{10}(T \ [\mathrm{K}])$')

        # Save:
        print('  * Figure generated successfully!')
        self.save(fig=fig, funcname=funcname)

        # Display figure:
        if self.show is not True:
            print('  * Interactive display of figure is NOT allowed')
            plt.close()
        else:
            print('  * Interactive display of figure is now running')
            plt.show()

        return None

    ##########################################################################
    ##########################################################################
    def resolution(self, bins=100, levels=5, smooth=0, funcname='resol'):
        '''
        Description: TODO
        '''
        # Read the data:
        h, iu = read_hdf5(file=self.file)
        mass = h['PartType0']['Masses'] * iu['umass']
        dens = h['PartType0']['Density'] * iu['udens']
        radius = ((3*mass) / (4*np.pi*dens))**(1/3)

        # 2D Histograms:
        H0, xedges0, yedges0 = np.histogram2d(np.log10(dens),
                                              np.log10(radius / const['pc']),
                                              bins=bins)
        H1, xedges1, yedges1 = np.histogram2d(np.log10(dens),
                                              np.log10(mass / const['Msol']),
                                              bins=bins)

        # Gaussian filter:
        if smooth > 0:
            H0 = gaussian_filter(H0, sigma=smooth)
            H1 = gaussian_filter(H1, sigma=smooth)

        # Log scale:
        with np.errstate(divide='ignore'):
            H0 = np.log10(H0.T)
            H1 = np.log10(H1.T)

        # Figure:
        fig, ax = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
        fig.subplots_adjust(left=0.15, right=0.85, bottom=0.15, top=0.95,
                            wspace=0, hspace=0)

        # Density vs Radius:
        cf = ax[0].contourf(H0, levels=levels, extent=(xedges0[0], xedges0[-1],
                            yedges0[0], yedges0[-1]), origin='lower',
                            cmap=configv.cmap['default'])
        ax[0].set_ylabel(r'$\log_{10}(r_\mathrm{cell} \ [\mathrm{pc}])$')
        ax[0].grid()
        ax[0].set_xlim(self.xlim)

        cax = ax[0].inset_axes([1, -1, 0.05, 2])
        fig.colorbar(cf, cax=cax, label=r'$\log_{10}(N_\mathrm{cells})$')

        # Density vs Mass:
        ax[1].contourf(H1, levels=levels, extent=(xedges1[0], xedges1[-1],
                       yedges1[0], yedges1[-1]), origin='lower',
                       cmap=configv.cmap['default'])
        ax[1].set_ylabel(r'$\log_{10}(M_\mathrm{cell} \ [\mathrm{M}_\odot])$')
        ax[1].set_xlabel(r'$\log_{10}(\rho_\mathrm{cell}$' +
                         r' $[\mathrm{g} \ \mathrm{cm}^{-3}])$')
        ax[1].grid()

        # Jeans length:
        ax[0].autoscale(False)
        xlim = ax[0].get_xlim()
        xlim = np.array(xlim)

        mu = (1 + 4 * 0.1)
        jeans10 = np.sqrt((15 * self.kb * 10) / (4 * np.pi * self.G * mu *
                                                 self.mp * 10**(xlim)))
        jeans100 = np.sqrt((15 * self.kb * 100) / (4 * np.pi * self.G * mu *
                                                   self.mp * 10**(xlim)))

        ax[0].plot(xlim, np.log10(jeans10 / const['pc']), c='k', ls='--',
                   lw=1, alpha=0.8, label=r'$\Lambda_J$($T=10$ K)')
        ax[0].plot(xlim, np.log10(jeans100 / const['pc']), c='k', ls='--',
                   lw=1, alpha=0.8, label=r'$\Lambda_J$($T=100$ K)')
        labelLines(ax[0].get_lines(), xvals=[-24, -24], ha='center',
                   fontsize=14)

        # Save:
        print('  * Figure generated successfully!')
        self.save(fig=fig, funcname=funcname)

        # Display figure:
        if self.show is not True:
            print('  * Interactive display of figure is NOT allowed')
            plt.close()
        else:
            print('  * Interactive display of figure is now running')
            plt.show()

        return None

    ##########################################################################
    ##########################################################################
    def stellar(self, axis='z', rotate=0, bins=100, bhfocus=False,
                funcname='stellar', xrange=None, yrange=None, zrange=None,
                box=None):
        '''
        Description: TODO
        '''
        # Read the data:
        print(f'  * Reading data of {self.file}')
        h, iu = read_hdf5(file=self.file)
        boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / const['kpc']
        time = h['Header'].attrs['Time'] * iu['utime'] / const['Myr']

        # Stellar component:
        pos_disk = h['PartType2']['Coordinates'] * iu['ulength'] / const['kpc']
        mass_disk = h['PartType2']['Masses'] * iu['umass'] / const['Msol']

        pos = pos_disk
        mass = mass_disk

        print('  * Generating a stellar surface density map')

        # Centre the data on the black hole:
        if bhfocus:
            bh = (h['PartType5']['Coordinates'][0] * iu['ulength']
                  / self.ulength)
            pos -= bh

        # Coordinate ranges:
        xrange, yrange, zrange = self.get_ranges(boxsize=boxsize, box=box,
                                                 xrange=xrange, yrange=yrange,
                                                 zrange=zrange,
                                                 bhfocus=bhfocus)

        # Bins in x and y:
        xbins, dx = np.linspace(xrange[0], xrange[1], bins, retstep=True)
        ybins, dy = np.linspace(yrange[0], yrange[1], bins, retstep=True)

        # Remove particles outside the zrange:
        mask = (pos[:, 2] > np.min(zrange)) * (pos[:, 2] < np.max(zrange))
        pos, mass = pos[mask], mass[mask]

        # Rotation of particle positions:
        pos = self.do_rotation(boxsize=boxsize, axis=axis, rotate=rotate,
                               pos=pos, bhfocus=bhfocus)

        # Histogram 2D:
        H, xedges, yedges = np.histogram2d(pos[:, 0], pos[:, 1],
                                           bins=(xbins, ybins), weights=mass)
        with np.errstate(divide='ignore'):
            H = np.log10(H.T / (dx * dy))

        # Plot:
        fig, ax = plt.subplots(figsize=(8, 6.4))
        fig.subplots_adjust(left=0.18, right=0.82, bottom=0.14, top=0.94,
                            wspace=0, hspace=0)

        im = ax.imshow(H, origin='lower',
                       extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                       vmin=self.vmin, vmax=self.vmax, cmap='bone')
        ax.text(0.95, 0.05, f'{round(time, 2)} Myr', color='k',
                ha='right', va='bottom', transform=ax.transAxes,
                bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.5,
                      'boxstyle': 'round'})
        bone = mpl.colormaps['bone']
        ax.set_facecolor(bone(0))
        ax.set_aspect('equal')
        ax.set_xlabel('$x$ [kpc]')
        ax.set_ylabel('$y$ [kpc]')
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)
        div = make_axes_locatable(ax)
        cax = div.append_axes('right', size='5%', pad=0)
        fig.colorbar(im, cax=cax, label=r'$\log_{10}$($\Sigma_\star$'
                     + r' [M$_\odot$ kpc$^{-2}$])')

        # Save:
        print('  * Figure generated successfully!')
        self.save(fig=fig, funcname=funcname)

        # Display figure:
        if self.show is not True:
            print('  * Interactive display of figure is NOT allowed')
            plt.close()
        else:
            print('  * Interactive display of figure is now running')
            plt.show()

        return None

    ##########################################################################
    ##########################################################################
    def darkmatter(self, axis='z', rotate=0, bins=100, bhfocus=False,
                   funcname='dm', xrange=None, yrange=None, zrange=None,
                   box=None):
        '''
        Description: TODO
        '''
        # Read the data:
        print(f'  * Reading data of {self.file}')
        h, iu = read_hdf5(file=self.file)
        boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / const['kpc']
        time = h['Header'].attrs['Time'] * iu['utime'] / const['Myr']

        # DM component:
        pos = h['PartType1']['Coordinates'] * iu['ulength'] / const['kpc']
        if 'Masses' in h['PartType1']:
            mass = h['PartType1']['Masses'] * iu['umass'] / const['Msol']
        else:
            mass = np.full(len(pos), h['Header'].attrs['MassTable'][1]
                           * iu['umass'] / const['Msol'])

        print('  * Generating a dark matter surface density map')

        # Centre the data on the black hole:
        if bhfocus:
            bh = (h['PartType5']['Coordinates'][0] * iu['ulength']
                  / self.ulength)
            pos -= bh

        # Coordinate ranges:
        xrange, yrange, zrange = self.get_ranges(boxsize=boxsize, box=box,
                                                 xrange=xrange, yrange=yrange,
                                                 zrange=zrange,
                                                 bhfocus=bhfocus)

        # Bins in x and y:
        xbins, dx = np.linspace(xrange[0], xrange[1], bins, retstep=True)
        ybins, dy = np.linspace(yrange[0], yrange[1], bins, retstep=True)

        # Remove particles outside the zrange:
        mask = (pos[:, 2] > np.min(zrange)) * (pos[:, 2] < np.max(zrange))
        pos, mass = pos[mask], mass[mask]

        # Rotation of particle positions:
        pos = self.do_rotation(boxsize=boxsize, axis=axis, rotate=rotate,
                               pos=pos, bhfocus=bhfocus)

        # Histogram 2D:
        H, xedges, yedges = np.histogram2d(pos[:, 0], pos[:, 1],
                                           bins=(xbins, ybins), weights=mass)
        with np.errstate(divide='ignore'):
            H = np.log10(H.T / (dx * dy))

        # Plot:
        fig, ax = plt.subplots(figsize=(8, 4), layout='constrained')
        im = ax.imshow(H, origin='lower',
                       extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                       vmin=self.vmin, vmax=self.vmax, cmap='magma')
        ax.text(0.95, 0.05, f'{round(time, 2)} Myr', color='k',
                ha='right', va='bottom', transform=ax.transAxes,
                bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.5,
                      'boxstyle': 'round'})
        bone = mpl.colormaps['magma']
        ax.set_facecolor(bone(0))
        ax.set_aspect('equal')
        ax.set_xlabel('$x$ [kpc]')
        ax.set_ylabel('$y$ [kpc]')
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)
        div = make_axes_locatable(ax)
        cax = div.append_axes('right', size='5%', pad=0)
        fig.colorbar(im, cax=cax, label=r'$\log_{10}$($\Sigma_\star$'
                     + r' [M$_\odot$ kpc$^{-2}$])')

        # Save:
        print('  * Figure generated successfully!')
        self.save(fig=fig, funcname=funcname)

        # Display figure:
        if self.show is not True:
            print('  * Interactive display of figure is NOT allowed')
            plt.close()
        else:
            print('  * Interactive display of figure is now running')
            plt.show()

        return None

    ##########################################################################
    ##########################################################################
    def star_formation(self, axis='z', rotate=0, bins=100, sfb=100,
                       interpolation='kdtree', bhfocus=False, funcname='sf',
                       xrange=None, yrange=None, zrange=None, box=None,
                       cut=None):
        '''
        Description: TODO
        '''
        # Read the data:
        print(f'  * Reading data of {self.file}')
        h, iu = read_hdf5(file=self.file)
        boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / self.ulength
        time = h['Header'].attrs['Time'] * iu['utime'] / const['Myr']
        pos_gas = h['PartType0']['Coordinates'] * iu['ulength'] / const['kpc']
        dens_gas = h['PartType0']['Density'] * iu['udens']
        pos_stars = (h['PartType4']['Coordinates'] * iu['ulength']
                     / const['kpc'])
        mass_stars = h['PartType4']['Masses'] * iu['umass'] / const['Msol']
        time_stars = (h['PartType4']['StellarFormationTime'] * iu['utime']
                      / const['Myr'])

        print('  * Generating a star formation rate surface density map')

        # Centre the data on the black hole:
        if bhfocus:
            bh = (h['PartType5']['Coordinates'][0] * iu['ulength']
                  / self.ulength)
            pos_gas -= bh
            pos_stars -= bh

        # Coordinate ranges:
        xrange, yrange, zrange = self.get_ranges(boxsize=boxsize, box=box,
                                                 xrange=xrange, yrange=yrange,
                                                 zrange=zrange,
                                                 bhfocus=bhfocus)

        # Rotation of particle positions:
        pos_gas = self.do_rotation(boxsize=boxsize, axis=axis, rotate=rotate,
                                   pos=pos_gas, bhfocus=bhfocus)
        pos_stars = self.do_rotation(boxsize=boxsize, axis=axis, rotate=rotate,
                                     pos=pos_stars, bhfocus=bhfocus)

        # Interpolation:
        if interpolation == 'kdtree':
            interpDens = interpolate_to_2d_kdtree(pos=pos_gas,
                                                  unit=self.ulength,
                                                  values=dens_gas, bins=bins,
                                                  xrange=xrange, yrange=yrange,
                                                  zrange=zrange, cut=cut)
        else:
            interpDens = interpolate_to_2d(pos=pos_gas, unit=self.ulength,
                                           values=dens_gas, bins=bins,
                                           xrange=xrange, yrange=yrange,
                                           zrange=zrange, cut=cut)

        # Star formation surface density:
        mask_sf = (time - time_stars < 10)
        xbins, dx = np.linspace(np.min(xrange), np.max(xrange), sfb,
                                retstep=True)
        ybins, dy = np.linspace(np.min(yrange), np.max(yrange), sfb,
                                retstep=True)
        H = np.histogram2d(pos_stars[:, 0][mask_sf], pos_stars[:, 1][mask_sf],
                           bins=[xbins, ybins], weights=mass_stars[mask_sf])[0]
        with np.errstate(divide='ignore'):
            H = np.log10(H.T/(dx * dy) / (10 * 1e6))

        # Plot:
        fig, ax = plt.subplots(figsize=(8, 6.4))
        fig.subplots_adjust(left=0.18, right=0.82, bottom=0.14, top=0.94,
                            wspace=0, hspace=0)

        im_gas = ax.imshow(np.log10(interpDens), origin='lower',
                           extent=(xrange[0], xrange[1], yrange[0], yrange[1]),
                           vmin=self.vmin, vmax=self.vmax, cmap='Greys')
        ax.text(0.95, 0.05, f'{round(time, 2)} Myr', color='k',
                ha='right', va='bottom', transform=ax.transAxes,
                bbox={'facecolor': 'white', 'edgecolor': 'none',
                      'alpha': 0.5, 'boxstyle': 'round'})
        ax.set_aspect('equal')
        ax.set_xlabel('$x$ [kpc]')
        ax.set_ylabel('$y$ [kpc]')
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)

        div = make_axes_locatable(ax)
        cax = div.append_axes('right', size='5%', pad=0)
        fig.colorbar(im_gas, cax=cax, label=r'$\log_{10}(\Sigma_\mathrm{Gas}$'
                     + r' $[\mathrm{g} \ \mathrm{cm}^{-2}])$')

        im_sf = ax.imshow(H, origin='lower', cmap='winter',
                          extent=(xrange[0], xrange[1], yrange[0], yrange[1]))
        cax = ax.inset_axes([0.02, 0.94, 0.7, 0.04])
        cb = fig.colorbar(im_sf, cax=cax, orientation='horizontal',
                          location='bottom')
        cb.set_label(label=r'$\log_{10}(\Sigma_\mathrm{SFR}$'
                     + r' $[\mathrm{M}_\odot$ yr$^{-1}$'
                     + f'{self.ulengthselect}' + '$^{-2}$' + '])', size=12)
        cb.ax.tick_params(labelsize=12)

        # Save:
        print('  * Figure generated successfully!')
        self.save(fig=fig, funcname=funcname)

        # Display figure:
        if self.show is not True:
            print('  * Interactive display of figure is NOT allowed')
            plt.close()
        else:
            print('  * Interactive display of figure is now running')
            plt.show()

        return None

    ##########################################################################
    ##########################################################################
    def stellar_age(self, axis='z', rotate=0, bins=100, age=100,
                    interpolation='kdtree', bhfocus=False, funcname='sa',
                    xrange=None, yrange=None, zrange=None, box=None, cut=None):
        '''
        Description: TODO
        '''
        # Read the data:
        print(f'  * Reading data of {self.file}')
        h, iu = read_hdf5(file=self.file)
        boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / self.ulength
        time = h['Header'].attrs['Time'] * iu['utime'] / const['Myr']
        pos_gas = h['PartType0']['Coordinates'] * iu['ulength'] / self.ulength
        dens_gas = h['PartType0']['Density'] * iu['udens']
        pos_stars = (h['PartType4']['Coordinates'] * iu['ulength']
                     / self.ulength)
        time_stars = (h['PartType4']['StellarFormationTime'] * iu['utime']
                      / const['Myr'])

        print('  * Generating a map highlighting the stellar age of newly' +
              ' formed star particles')

        # Centre the data on the black hole:
        if bhfocus:
            bh = (h['PartType5']['Coordinates'][0] * iu['ulength']
                  / self.ulength)
            pos_gas -= bh
            pos_stars -= bh

        # Coordinate ranges:
        xrange, yrange, zrange = self.get_ranges(boxsize=boxsize, box=box,
                                                 xrange=xrange, yrange=yrange,
                                                 zrange=zrange,
                                                 bhfocus=bhfocus)

        # Rotation of particle positions:
        pos_gas = self.do_rotation(boxsize=boxsize, axis=axis, rotate=rotate,
                                   pos=pos_gas, bhfocus=bhfocus)
        pos_stars = self.do_rotation(boxsize=boxsize, axis=axis, rotate=rotate,
                                     pos=pos_stars, bhfocus=bhfocus)

        # Interpolation:
        if interpolation == 'kdtree':
            interpDens = interpolate_to_2d_kdtree(pos=pos_gas,
                                                  unit=self.ulength,
                                                  values=dens_gas, bins=bins,
                                                  xrange=xrange, yrange=yrange,
                                                  zrange=zrange, cut=cut)
        else:
            interpDens = interpolate_to_2d(pos=pos_gas, unit=self.ulength,
                                           values=dens_gas, bins=bins,
                                           xrange=xrange, yrange=yrange,
                                           zrange=zrange, cut=cut)

        # Stars:
        mask = ((pos_stars[:, 0] > xrange[0]) * (pos_stars[:, 0] < xrange[1])
                * (pos_stars[:, 1] > yrange[0]) * (pos_stars[:, 1] < yrange[1])
                * (pos_stars[:, 2] > zrange[0]) * (pos_stars[:, 2] < zrange[1])
                * (np.abs(time - time_stars) < age))

        # Plot:
        fig, ax = plt.subplots(figsize=(8, 6.4))
        fig.subplots_adjust(left=0.18, right=0.82, bottom=0.14, top=0.94,
                            wspace=0, hspace=0)

        im_gas = ax.imshow(np.log10(interpDens), origin='lower', cmap='Greys',
                           extent=(xrange[0], xrange[1], yrange[0], yrange[1]),
                           vmin=self.vmin, vmax=self.vmax)
        ax.text(0.95, 0.05, f'{round(time, 2)} Myr', color='k',
                ha='right', va='bottom', transform=ax.transAxes,
                bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.5,
                      'boxstyle': 'round'})
        ax.set_aspect('equal')
        ax.set_xlabel('$x$ [kpc]')
        ax.set_ylabel('$y$ [kpc]')
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)

        div = make_axes_locatable(ax)
        cax = div.append_axes('right', size='5%', pad=0)
        fig.colorbar(im_gas, cax=cax, label=r'$\log_{10}(\Sigma_\mathrm{Gas}$'
                     + r' $[\mathrm{g} \ \mathrm{cm}^{-2}])$')

        time_diff = np.abs(time - time_stars[mask])
        im_sa = ax.scatter(pos_stars[:, 0][mask], pos_stars[:, 1][mask],
                           c=time_diff, s=10, marker='.', cmap='viridis',
                           vmin=0, vmax=np.min([age, np.max(time_diff)]))
        cax = ax.inset_axes([0.02, 0.94, 0.7, 0.04])
        cb = fig.colorbar(im_sa, cax=cax, orientation='horizontal',
                          location='bottom')
        cb.set_label(label='Stellar age [Myr]', size=12)
        cb.ax.tick_params(labelsize=12)

        # Save:
        print('  * Figure generated successfully!')
        self.save(fig=fig, funcname=funcname)

        # Display figure:
        if self.show is not True:
            print('  * Interactive display of figure is NOT allowed')
            plt.close()
        else:
            print('  * Interactive display of figure is now running')
            plt.show()

        return None

    ##########################################################################
    ##########################################################################
    def black_hole_evolution(self, vcr, funcname='bhevol'):
        '''
        Description: TODO
        '''
        # Snapshot range:
        file_list = os.listdir()
        snap_list = [i for i in file_list if 'snap_' in i]
        snap_list = [i for i in snap_list if 'sink_' not in i]
        snum_list = [int(i[5:-5]) for i in snap_list]
        n = np.min(snum_list)
        N = int(self.file[len(self.file)-8:-5])

        # Obtain the black hole data:
        BHData = get_black_hole_data(output_dir=os.getcwd(), n=n, N=N,
                                     vcr=vcr)

        # Collection of plots:
        ls = '-'
        lw = 2
        c = 'tab:blue'
        fs = 14

        fig, ax = plt.subplots(2, 5, figsize=(14, 5.5), sharex=True,
                               layout='constrained')

        # Sink mass:
        ax[0, 0].plot(BHData['Time'], np.log10(BHData['MassSink']),
                      ls=ls, lw=lw, c=c)
        ax[0, 0].set_title(r'$\log_{10}(M_\mathrm{Sink}$' +
                           r' $[\mathrm{M}_\odot])$', fontsize=fs)
        ax[0, 0].set_xlim(0, np.max(BHData['Time']))
        ax[0, 0].grid()

        # BH growth:
        ax[0, 1].plot(BHData['Time'], np.log10(BHData['MassBH']),
                      ls=ls, lw=lw, c=c)
        ax[0, 1].set_title(r'$\log_{10}(M_\mathrm{BH} \ [\mathrm{M}_\odot])$',
                           fontsize=fs)
        ax[0, 1].grid()

        # Gas reservoir:
        GasReserv = np.array(BHData['MassReserv'])
        GasReserv[GasReserv <= 0] = 1e-99
        ax[0, 2].plot(BHData['Time'], np.log10(GasReserv), ls=ls, lw=lw, c=c)
        ax[0, 2].set_title(r'$\log_{10}(M_\mathrm{Reserv}$' +
                           r' $[\mathrm{M}_\odot])$', fontsize=fs)
        ax[0, 2].set_ylim(-9, 5)
        ax[0, 2].grid()

        # Gas accretion disk:
        AccDisk = np.array(BHData['MassDisk'])
        AccDisk[AccDisk <= 0] = 1e-99
        ax[0, 3].plot(BHData['Time'], np.log10(AccDisk), ls=ls, lw=lw, c=c)
        ax[0, 3].set_title(r'$\log_{10}(M_\mathrm{Disk}$' +
                           r' $[\mathrm{M}_\odot])$', fontsize=fs)
        ax[0, 3].set_ylim(-9, 5)
        ax[0, 3].grid()

        # Relative error:
        mass_diff = ((np.array(BHData['MassSink'])
                      - np.array(BHData['MassReserv'])
                      - np.array(BHData['MassDisk'])
                      - np.array(BHData['MassBH']))
                     / np.array(BHData['MassSink']))
        ax[0, 4].plot(BHData['Time'], mass_diff, lw=lw, c=c)
        ax[0, 4].set_title('Relative Error', fontsize=fs)
        ax[0, 4].set_yscale('linear')
        ax[0, 4].grid()

        # Sink accretion rate:
        FracEddSink = np.array(BHData['MdotSink'])/np.array(BHData['MdotEdd'])
        FracEddSink[FracEddSink <= 0] = 1e-99
        ax[1, 0].stairs(np.log10(FracEddSink), BHData['Time'], baseline=-99,
                        lw=lw, color=c, alpha=0.8, fill=True, rasterized=True)
        ax[1, 0].axhline(0, c='k', ls=':', lw=1, zorder=9)
        ax[1, 0].axhline(np.log10(0.02), c='k', ls='--', lw=1, zorder=9)
        ax[1, 0].set_xlabel('Time [Myr]')
        ax[1, 0].set_title(r'$\log_{10}(\dot{M}_\mathrm{Sink}/$' +
                           r'$\dot{M}_\mathrm{Edd})$', fontsize=fs)
        ax[1, 0].set_ylim(-9, 1)
        ax[1, 0].grid()

        MdotSink = np.array(BHData['MdotSink'])
        MdotSink[MdotSink <= 0] = 1e-99
        ax[1, 1].stairs(np.log10(MdotSink), BHData['Time'], baseline=-99,
                        lw=lw, color=c, alpha=0.8, fill=True, rasterized=True)
        ax[1, 1].set_xlabel('Time [Myr]')
        ax[1, 1].set_title(r'$\log_{10}(\dot{M}_\mathrm{Sink}$' +
                           r'$ [\mathrm{M}_\odot \ \mathrm{yr}^{-1}])$',
                           fontsize=fs)
        ax[1, 1].set_ylim(-9, -3)
        ax[1, 1].grid()

        # BH Accretion rate:
        FracEddBH = np.array(BHData['MdotBH']) / np.array(BHData['MdotEdd'])
        FracEddBH[FracEddBH <= 0] = 1e-99
        ax[1, 2].plot(BHData['TimeMid'], np.log10(FracEddBH), lw=lw, c=c,
                      zorder=10)
        ax[1, 2].axhline(0, c='k', ls=':', lw=1, zorder=9)
        ax[1, 2].axhline(np.log10(0.02), c='k', ls='--', lw=1, zorder=9)
        ax[1, 2].set_xlabel('Time [Myr]')
        ax[1, 2].set_title(r'$\log_{10}(\dot{M}_\mathrm{BH}/$' +
                           r'$\dot{M}_\mathrm{Edd})$', fontsize=fs)
        ax[1, 2].set_ylim(-9, 1)
        ax[1, 2].grid()

        MdotBH = np.array(BHData['MdotBH'])
        MdotBH[MdotBH <= 0] = 1e-99
        ax[1, 3].plot(BHData['TimeMid'], np.log10(MdotBH), lw=2, c=c,
                      zorder=10)
        ax[1, 3].set_xlabel('Time [Myr]')
        ax[1, 3].set_title(r'$\log_{10}(\dot{M}_\mathrm{BH}$' +
                           r' [M$_\odot$ yr$^{-1}$])', fontsize=fs)
        ax[1, 3].set_ylim(-9, -3)
        ax[1, 3].grid()

        # Circularisation radius:
        if len(BHData['CircRadius']) > 0:
            ax[1, 4].plot(BHData['Time'], BHData['CircRadius'], lw=lw, c=c)
            ax[1, 4].set_xlabel('Time [Myr]')
            ax[1, 4].set_title(r'$R_\mathrm{circ}$ [left: pc, right: ' +
                               r'$r_\mathrm{s}$]', fontsize=fs)
            ax[1, 4].set_yscale('linear')
            ax[1, 4].grid()

            ax2 = ax[1, 4].twinx()
            rs = (2 * const['G'] * BHData['MassBH'][0] * const['Msol'] /
                  np.power(const['c'], 2))
            Rcirc = np.array(BHData['CircRadius']) * const['pc'] / rs
            ax2.plot(BHData['Time'], Rcirc, lw=lw, c=c)
            ax2.set_yscale('linear')
        else:
            ax[1, 4].set_visible(False)

        # Save:
        print('  * Figure generated successfully!')
        self.save(fig=fig, funcname=funcname)

        # Display figure:
        if self.show is not True:
            print('  * Interactive display of figure is NOT allowed')
            plt.close()
        else:
            print('  * Interactive display of figure is now running')
            plt.show()

        return None

    ##########################################################################
    ##########################################################################
    def ffmpeg(self, framedir, skip):
        '''
        Description: TODO
        '''

        if os.path.isdir(f'./vframes/{framedir}'):
            fnr = 0
            if skip:
                fnr = skip
            frame = '000'[:3-len(str(fnr))] + str(fnr)
            print(f'  * Searching for {framedir} frames')
            while os.path.isfile(f'./vframes/{framedir}/' +
                                 f'{framedir}_{frame}.png'):
                fnr += 1
                frame = '000'[:3-len(str(fnr))] + str(fnr)
            print(f'  * Found {fnr - skip} frames!')

        start = 0
        file_split = self.file.split('.')
        vframes = int(file_split[0][-3:])
        if skip:
            start = skip
            vframes = vframes - skip + 1
        if not skip:
            vframes += 1
        if vframes > fnr:
            print('  * Trying to generate a movie up to a given snapshot' +
                  ' without the neccessary frames needed to do so,' +
                  ' therefore, exiting the ffmpeg function')
            return None
        fps = 15
        res_width = 1920
        res_height = 1080

        f = os.system(f'ffmpeg -r {fps} -f image2 -y' +
                      f' -s {res_width}x{res_height}' +
                      f' -start_number {start}' +
                      f' -i {os.getcwd()}/vframes/{framedir}/' +
                      f'{framedir}_%03d.png' +
                      f' -vframes {vframes}' +
                      ' -vcodec libx264' +
                      ' -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white"' +
                      ' -crf 25' +
                      ' -pix_fmt yuv420p' +
                      f' {os.getcwd()}/{framedir}.mp4' +
                      ' > ffmpeg.out 2> ffmpeg.err')
        if f == 0:
            print('  * Film generated successfully!')
        else:
            print('  * Error: Something went wrong, please check the ffmpeg' +
                  ' output/error file')

        return None

# -------------- End of file
