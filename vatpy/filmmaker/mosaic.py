# Description: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-10-10


# -------------- Required packages
import numpy as np
import matplotlib.pyplot as plt

from ..read import read_hdf5
from ..constants import const
from ..get_data.get_image_data import get_image_data

import configv


# -------------- Declare function(s)
def plot_mosaic(start, end):
    '''TODO
    '''
    plt.style.use(f'{configv.homedir}/vatpy/mpl/mosaic.mplstyle')

    # Parameters:
    bins = 200
    bins_sf = 100
    bins_stars = 200
    bins_stellar = 200
    box = (-2, 2)

    vmin_dens, vmax_dens = -5, -2
    vmin_temp, vmax_temp = 2.5, 5.5
    vmin_bfield, vmax_bfield = None, None

    vmin_HI, vmax_HI = 19, 21.5
    vmin_HII, vmax_HII = 18, 20
    vmin_H2, vmax_H2 = 12, 19

    vmin_sf, vmax_sf = -3, 0
    vmin_stars, vmax_stars = 5, 6.5
    vmin_stellar, vmax_stellar = 5, 7

    cmap_dens = 'inferno'
    cmap_temp = 'cmr.amber'
    cmap_HI = 'cmr.rainforest'
    cmap_HII = 'cmr.ember'
    cmap_H2 = 'cmr.cosmic'
    cmap_sf = 'winter'
    cmap_stars = 'autumn'
    cmap_stellar = 'bone'
    cmap_age = 'cmr.guppy_r'
    cmap_bfield = 'plasma'

    # Loop over snapshots:
    for i in range(start, end+1):
        snap = '000'[:3-len(str(i))] + str(i)
        file = 'snap_' + snap + '.hdf5'
        print(f'  * Generating frame for snapshot {i}')
        h, iu = read_hdf5(file=file)
        time = h['Header'].attrs['Time'] * iu['utime'] / const['Myr']
        boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / const['kpc']

        # Gaseous component:
        quantity = ['mass', 'temp', 'HI', 'H2', 'HII', 'bfield']
        im_gas_face = get_image_data(file=file, bins=bins, box=box,
                                     quantity=quantity, halfboxfocus=True)
        im_gas_edge = get_image_data(file=file, bins=bins, box=box,
                                     quantity=quantity, axis='x', rotate=90,
                                     halfboxfocus=True)

        # Newly formed stars:
        bin_edges_sf, bin_sizes_sf = np.linspace(box[0], box[1], bins_sf,
                                                 retstep=True)
        bin_edges_stars, bin_sizes_stars = np.linspace(box[0], box[1],
                                                       bins_stars,
                                                       retstep=True)

        if 'PartType4' in h.keys():
            pos_stars = (h['PartType4']['Coordinates'] * iu['ulength'] /
                         const['kpc'])
            pos_stars -= boxsize / 2
            mass_stars = h['PartType4']['Masses'] * iu['umass'] / const['Msol']
            if 'InitialMass' in h['PartType4'].keys():
                birthmass = (h['PartType4']['InitialMass'] * iu['umass'] /
                             const['Msol'])
            else:
                birthmass = (h['PartType4']['Masses'] * iu['umass'] /
                             const['Msol'])
            birthtime = (h['PartType4']['StellarFormationTime'] * iu['utime'] /
                         const['Myr'])
            mask_sf = (time - birthtime) < 10

            H_stars_face = np.histogram2d(pos_stars[:, 0], pos_stars[:, 1],
                                          bins=bin_edges_stars,
                                          weights=mass_stars)[0]
            H_sf_face = np.histogram2d(pos_stars[:, 0][mask_sf],
                                       pos_stars[:, 1][mask_sf],
                                       bins=bin_edges_sf,
                                       weights=birthmass[mask_sf])[0]
            H_stars_edge = np.histogram2d(pos_stars[:, 0], pos_stars[:, 2],
                                          bins=bin_edges_stars,
                                          weights=mass_stars)[0]
            H_sf_edge = np.histogram2d(pos_stars[:, 0][mask_sf],
                                       pos_stars[:, 2][mask_sf],
                                       bins=bin_edges_sf,
                                       weights=birthmass[mask_sf])[0]
            with np.errstate(divide='ignore'):
                H_stars_face = np.log10(H_stars_face.T /
                                        (bin_sizes_stars * bin_sizes_stars))
                H_sf_face = np.log10(H_sf_face.T /
                                     (bin_sizes_sf * bin_sizes_sf) / (1e7))
                H_stars_edge = np.log10(H_stars_edge.T /
                                        (bin_sizes_stars * bin_sizes_stars))
                H_sf_edge = np.log10(H_sf_edge.T /
                                     (bin_sizes_sf * bin_sizes_sf) / (1e7))
        else:
            H_stars_face = np.full(shape=(len(bin_edges_stars[:-1]),
                                          len(bin_edges_stars[:-1])),
                                   fill_value=np.nan)
            H_sf_face = np.full(shape=(len(bin_edges_sf[:-1]),
                                       len(bin_edges_sf[:-1])),
                                fill_value=np.nan)
            H_stars_edge = np.full(shape=(len(bin_edges_stars[:-1]),
                                          len(bin_edges_stars[:-1])),
                                   fill_value=np.nan)
            H_sf_edge = np.full(shape=(len(bin_edges_sf[:-1]),
                                       len(bin_edges_sf[:-1])),
                                fill_value=np.nan)

        # Stellar component:
        bin_edges_stellar, bin_sizes_stellar = np.linspace(box[0], box[1],
                                                           bins_stellar,
                                                           retstep=True)

        if 'PartType2' in h.keys():
            pos_stellar = (h['PartType2']['Coordinates'] * iu['ulength'] /
                           const['kpc'])
            pos_stellar -= boxsize / 2
            mass_stellar = (h['PartType2']['Masses'] * iu['umass'] /
                            const['Msol'])

            H_stellar_face = np.histogram2d(pos_stellar[:, 0],
                                            pos_stellar[:, 1],
                                            bins=bin_edges_stellar,
                                            weights=mass_stellar)[0]
            H_stellar_edge = np.histogram2d(pos_stellar[:, 0],
                                            pos_stellar[:, 2],
                                            bins=bin_edges_stellar,
                                            weights=mass_stellar)[0]
            with np.errstate(divide='ignore'):
                H_stellar_face = np.log10(H_stellar_face.T /
                                          (bin_sizes_stellar *
                                           bin_sizes_stellar))
                H_stellar_edge = np.log10(H_stellar_edge.T /
                                          (bin_sizes_stellar *
                                           bin_sizes_stellar))
        else:
            H_stellar_face = np.full(shape=(bin_edges_stellar,
                                            bin_edges_stellar),
                                     fill_value=np.nan)
            H_stellar_edge = np.full(shape=(bin_edges_stellar,
                                            bin_edges_stellar),
                                     fill_value=np.nan)

        # Figure:
        fig, ax = plt.subplots(4, 5, figsize=(16, 3/5 * 16),
                               gridspec_kw={'height_ratios': [1, 1/2, 1, 1/2],
                                            'width_ratios': [1, 1, 1, 1, 1]})
        fig.subplots_adjust(left=0.06, right=0.94, bottom=0.06, top=0.94,
                            wspace=0, hspace=0)

        # Density:
        im = ax[0, 0].imshow(np.log10(im_gas_face['mass']),
                             extent=im_gas_face['extent'], origin='lower',
                             vmin=vmin_dens, vmax=vmax_dens, cmap=cmap_dens)
        ax[0, 0].set_aspect('equal')
        ax[0, 0].set_xticks([])
        yticks = ax[0, 0].get_yticks()
        ax[0, 0].set_yticks(yticks[1:-1])
        ax[0, 0].set_ylabel('$y$ [kpc]')
        cax = ax[0, 0].inset_axes([0.05, 0.93, 0.6, 0.03])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(color='w', label=r'$\log_{10}(\Sigma_\text{Gas} \ $' +
                       r'$[\text{g} \ \text{cm}^{-2}])$')
        cbar.ax.tick_params(colors='w', which='both')
        cbar.outline.set_edgecolor('w')

        im = ax[1, 0].imshow(np.log10(im_gas_edge['mass']),
                             extent=im_gas_edge['extent'], origin='lower',
                             vmin=vmin_dens, vmax=vmax_dens, cmap=cmap_dens)
        ax[1, 0].set_ylim(box[0]/2, box[1]/2)
        ax[1, 0].set_aspect('equal')
        ax[1, 0].set_xticks([])
        ax[1, 0].set_ylabel('$z$ [kpc]')

        # Temperature:
        im = ax[0, 1].imshow(np.log10(im_gas_face['temp']),
                             extent=im_gas_face['extent'], origin='lower',
                             vmin=vmin_temp, vmax=vmax_temp, cmap=cmap_temp)
        ax[0, 1].set_aspect('equal')
        ax[0, 1].set_xticks([])
        ax[0, 1].set_yticks([])
        cax = ax[0, 1].inset_axes([0.05, 0.93, 0.6, 0.03])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(color='w', label=r'$\log_{10}(T \ [\mathrm{K}])$')
        cbar.ax.tick_params(colors='w', which='both')
        cbar.outline.set_edgecolor('w')

        im = ax[1, 1].imshow(np.log10(im_gas_edge['temp']),
                             extent=im_gas_edge['extent'], origin='lower',
                             vmin=vmin_temp, vmax=vmax_temp, cmap=cmap_temp)
        ax[1, 1].set_ylim(box[0]/2, box[1]/2)
        ax[1, 1].set_aspect('equal')
        ax[1, 1].set_xticks([])
        ax[1, 1].set_yticks([])

        # Star formation:
        im = ax[0, 2].imshow(np.log10(im_gas_face['mass']),
                             extent=im_gas_face['extent'], origin='lower',
                             vmin=vmin_dens, vmax=vmax_dens, cmap='Greys')
        ax[0, 2].set_aspect('equal')
        ax[0, 2].set_xticks([])
        ax[0, 2].set_yticks([])

        im = ax[1, 2].imshow(np.log10(im_gas_edge['mass']),
                             extent=im_gas_edge['extent'], origin='lower',
                             vmin=vmin_dens, vmax=vmax_dens, cmap='Greys')
        ax[1, 2].set_ylim(box[0]/2, box[1]/2)
        ax[1, 2].set_aspect('equal')
        ax[1, 2].set_xticks([])
        ax[1, 2].set_yticks([])

        im = ax[0, 2].pcolormesh(bin_edges_sf, bin_edges_sf, H_sf_face,
                                 vmin=vmin_sf, vmax=vmax_sf, cmap=cmap_sf)
        cax = ax[0, 2].inset_axes([0.05, 0.93, 0.6, 0.03])
        fig.colorbar(im, cax=cax, orientation='horizontal',
                     label=r'$\log_{10}(\Sigma_\text{SFR} \ $' +
                     r'$[\text{M}_\odot \ \text{yr}^{-1} \ \text{kpc}^{-2}])$')
        ax[1, 2].pcolormesh(bin_edges_sf, bin_edges_sf, H_sf_edge,
                            vmin=vmin_sf, vmax=vmax_sf, cmap=cmap_sf)

        # Star particles:
        im = ax[0, 3].imshow(np.log10(im_gas_face['mass']),
                             extent=im_gas_face['extent'], origin='lower',
                             vmin=vmin_dens, vmax=vmax_dens, cmap='Greys')
        ax[0, 3].set_aspect('equal')
        ax[0, 3].set_xticks([])
        ax[0, 3].set_yticks([])

        im = ax[1, 3].imshow(np.log10(im_gas_edge['mass']),
                             extent=im_gas_edge['extent'], origin='lower',
                             vmin=vmin_dens, vmax=vmax_dens, cmap='Greys')
        ax[1, 3].set_ylim(box[0]/2, box[1]/2)
        ax[1, 3].set_aspect('equal')
        ax[1, 3].set_xticks([])
        ax[1, 3].set_yticks([])

        im = ax[0, 3].pcolormesh(bin_edges_stars, bin_edges_stars,
                                 H_stars_face,
                                 vmin=vmin_stars, vmax=vmax_stars,
                                 cmap=cmap_stars)
        cax = ax[0, 3].inset_axes([0.05, 0.93, 0.6, 0.03])
        fig.colorbar(im, cax=cax, orientation='horizontal',
                     label=r'$\log_{10}(\Sigma_\text{StarP} \ $' +
                     r'$[\text{M}_\odot \ \text{kpc}^{-2}])$')
        ax[1, 3].pcolormesh(bin_edges_stars, bin_edges_stars, H_stars_edge,
                            vmin=vmin_stars, vmax=vmax_stars, cmap=cmap_stars)

        # Stellar age:
        im = ax[0, 4].imshow(np.log10(im_gas_face['mass']),
                             extent=im_gas_face['extent'], origin='lower',
                             vmin=vmin_dens, vmax=vmax_dens, cmap='Greys')
        ax[0, 4].set_xlim(box[0], box[1])
        ax[0, 4].set_ylim(box[0], box[1])
        ax[0, 4].set_aspect('equal')
        ax[0, 4].set_xticks([])
        ax[0, 4].set_yticks([])

        im = ax[1, 4].imshow(np.log10(im_gas_edge['mass']),
                             extent=im_gas_edge['extent'], origin='lower',
                             vmin=vmin_dens, vmax=vmax_dens, cmap='Greys')
        ax[1, 4].set_xlim(box[0], box[1])
        ax[1, 4].set_ylim(box[0]/2, box[1]/2)
        ax[1, 4].set_aspect('equal')
        ax[1, 4].set_xticks([])
        ax[1, 4].set_yticks([])

        if 'PartType4' in h.keys():
            mask_age = (time - birthtime) < 100
            age = (time - birthtime)[mask_age]
            x_stars, y_stars, z_stars = pos_stars[mask_age].T

            im = ax[0, 4].scatter(x_stars, y_stars, c=age, vmin=0, vmax=100,
                                  cmap=cmap_age, s=4, marker='.')
            cax = ax[0, 4].inset_axes([0.05, 0.93, 0.6, 0.03])
            fig.colorbar(im, cax=cax, orientation='horizontal',
                         label='Stellar Age [Myr]')

            im = ax[1, 4].scatter(x_stars, z_stars, c=age, vmin=0, vmax=100,
                                  cmap=cmap_age, s=4, marker='.')

        # HI surface density:
        im = ax[2, 0].imshow(np.log10(im_gas_face['HI']),
                             extent=im_gas_face['extent'], origin='lower',
                             vmin=vmin_HI, vmax=vmax_HI, cmap=cmap_HI)
        ax[2, 0].set_facecolor(plt.get_cmap(cmap_HI)(0))
        cax = ax[2, 0].inset_axes([0.05, 0.93, 0.6, 0.03])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(color='w', label=r'$\log_{10}(\Sigma_\text{HI} \ $' +
                       r'$[\text{cm}^{-2}])$')
        cbar.ax.tick_params(colors='w', which='both')
        cbar.outline.set_edgecolor('w')
        ax[2, 0].set_aspect('equal')
        ax[2, 0].set_xticks([])
        ax[2, 0].set_yticks(yticks[1:-1])
        ax[2, 0].set_ylabel(r'$y$ [kpc]')

        ax[3, 0].imshow(np.log10(im_gas_edge['HI']),
                        extent=im_gas_edge['extent'], origin='lower',
                        vmin=vmin_HI, vmax=vmax_HI, cmap=cmap_HI)
        ax[3, 0].set_facecolor(plt.get_cmap(cmap_HI)(0))
        ax[3, 0].set_ylim(box[0]/2, box[1]/2)
        ax[3, 0].set_aspect('equal')
        xticks = ax[3, 0].get_xticks()
        ax[3, 0].set_xticks(xticks[1:-1])
        ax[3, 0].set_xlabel(r'$x$ [kpc]')
        ax[3, 0].set_ylabel(r'$z$ [kpc]')

        # HII surface density:
        im = ax[2, 1].imshow(np.log10(im_gas_face['HII']),
                             extent=im_gas_face['extent'], origin='lower',
                             vmin=vmin_HII, vmax=vmax_HII, cmap=cmap_HII)
        ax[2, 1].set_facecolor(plt.get_cmap(cmap_HII)(0))
        cax = ax[2, 1].inset_axes([0.05, 0.93, 0.6, 0.03])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(color='w', label=r'$\log_{10}(\Sigma_\text{HII} \ $' +
                       r'$[\text{cm}^{-2}])$')
        cbar.ax.tick_params(colors='w', which='both')
        cbar.outline.set_edgecolor('w')
        ax[2, 1].set_aspect('equal')
        ax[2, 1].set_xticks([])
        ax[2, 1].set_yticks([])

        ax[3, 1].imshow(np.log10(im_gas_edge['HII']),
                        extent=im_gas_edge['extent'], origin='lower',
                        vmin=vmin_HII, vmax=vmax_HII, cmap=cmap_HII)
        ax[3, 1].set_facecolor(plt.get_cmap(cmap_HII)(0))
        ax[3, 1].set_ylim(box[0]/2, box[1]/2)
        ax[3, 1].set_aspect('equal')
        ax[3, 1].set_yticks([])
        ax[3, 1].set_xticks(xticks[1:-1])
        ax[3, 1].set_xlabel(r'$x$ [kpc]')

        # H2 surface density:
        im = ax[2, 2].imshow(np.log10(im_gas_face['H2']),
                             extent=im_gas_face['extent'], origin='lower',
                             vmin=vmin_H2, vmax=vmax_H2, cmap=cmap_H2)
        ax[2, 2].set_facecolor(plt.get_cmap(cmap_H2)(0))
        cax = ax[2, 2].inset_axes([0.05, 0.93, 0.6, 0.03])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(color='w', label=r'$\log_{10}(\Sigma_{\text{H}_2} \ $' +
                       r'$[\text{cm}^{-2}])$')
        cbar.ax.tick_params(colors='w', which='both')
        cbar.outline.set_edgecolor('w')
        ax[2, 2].set_aspect('equal')
        ax[2, 2].set_xticks([])
        ax[2, 2].set_yticks([])

        ax[3, 2].imshow(np.log10(im_gas_edge['H2']),
                        extent=im_gas_edge['extent'], origin='lower',
                        vmin=vmin_H2, vmax=vmax_H2, cmap=cmap_H2)
        ax[3, 2].set_facecolor(plt.get_cmap(cmap_H2)(0))
        ax[3, 2].set_ylim(box[0]/2, box[1]/2)
        ax[3, 2].set_aspect('equal')
        ax[3, 2].set_xticks(xticks[1:-1])
        ax[3, 2].set_yticks([])
        ax[3, 2].set_xlabel(r'$x$ [kpc]')

        # Stellar component:
        im = ax[2, 3].pcolormesh(bin_edges_stellar, bin_edges_stellar,
                                 H_stellar_face, vmin=vmin_stellar,
                                 vmax=vmax_stellar, cmap=cmap_stellar)
        ax[2, 3].set_facecolor(plt.get_cmap(cmap_stellar)(0))
        ax[2, 3].set_aspect('equal')
        ax[2, 3].set_xticks([])
        ax[2, 3].set_yticks([])
        cax = ax[2, 3].inset_axes([0.05, 0.93, 0.6, 0.03])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(color='w', label=r'$\log_{10}(\Sigma_\star \ $' +
                       r'$[\mathrm{M}_\odot \ \mathrm{kpc}^{-2}])$')
        cbar.ax.tick_params(colors='w', which='both')
        cbar.outline.set_edgecolor('w')

        ax[3, 3].pcolormesh(bin_edges_stellar, bin_edges_stellar,
                            H_stellar_edge, vmin=vmin_stellar,
                            vmax=vmax_stellar, cmap=cmap_stellar)
        ax[3, 3].set_facecolor(plt.get_cmap(cmap_stellar)(0))
        ax[3, 3].set_ylim(box[0]/2, box[1]/2)
        ax[3, 3].set_aspect('equal')
        ax[3, 3].set_xticks(xticks[1:-1])
        ax[3, 3].set_yticks([])
        ax[3, 3].set_xlabel('$x$ [kpc]')

        # Magnetic field strength:
        im = ax[2, 4].imshow(np.log10(im_gas_face['bfield'] / 1e-6),
                             extent=im_gas_face['extent'], origin='lower',
                             vmin=vmin_bfield, vmax=vmax_bfield,
                             cmap=cmap_bfield)
        ax[2, 4].set_facecolor(plt.get_cmap(cmap_bfield)(0))
        ax[2, 4].set_aspect('equal')
        ax[2, 4].set_xticks([])
        ax[2, 4].set_yticks([])
        cax = ax[2, 4].inset_axes([0.05, 0.93, 0.6, 0.03])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(color='w', label=r'$\log_{10}(|\vec{B}| \ [\mu G])$')
        cbar.ax.tick_params(colors='w', which='both')
        cbar.outline.set_edgecolor('w')

        im = ax[3, 4].imshow(np.log10(im_gas_edge['bfield'] / 1e-6),
                             extent=im_gas_edge['extent'], origin='lower',
                             vmin=vmin_bfield, vmax=vmax_bfield,
                             cmap=cmap_bfield)
        ax[3, 4].set_facecolor(plt.get_cmap(cmap_bfield)(0))
        ax[3, 4].set_ylim(box[0]/2, box[1]/2)
        ax[3, 4].set_aspect('equal')
        ax[3, 4].set_xticks(xticks[1:-1])
        ax[3, 4].set_yticks([])
        ax[3, 4].set_xlabel('$x$ [kpc]')

        # Time:
        ax[0, 4].text(1, 1, f'{time:.2f} Myr', fontsize=12, color='k',
                      ha='right', va='bottom', transform=ax[0, 4].transAxes)

        fig.savefig(f'./vframes/mosaic/mosaic_{snap}.png')
        plt.close()

    return None

# -------------- End of file
