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
    plt.style.use(f'{configv.homedir}/vatpy/mpl/noctua.mplstyle')

    # Parameters:
    bins = 300
    bins_sf = 200
    bins_stars = 200
    bins_stellar = 200
    box = (-2, 2)

    vmin_dens, vmax_dens = -5, -2
    vmin_temp, vmax_temp = 2, 6

    vmin_HI, vmax_HI = 19, 22
    vmin_HII, vmax_HII = 18, 20
    vmin_H2, vmax_H2 = 12.5, 20

    vmin_sf, vmax_sf = -3, 0
    vmin_stars, vmax_stars = 5, 6
    vmin_stellar, vmax_stellar = 4, 6

    cmap_dens = 'cmr.rainforest'
    cmap_temp = 'cmr.amber'
    cmap_sf = 'winter'
    cmap_stars = 'summer'
    cmap_age = 'cmr.guppy_r'

    # Loop over snapshots:
    for i in range(start, end+1):
        snap = '000'[:3-len(str(i))] + str(i)
        file = 'snap_' + snap + '.hdf5'
        print(f'  * Generating frame for snapshot {i}')
        h, iu = read_hdf5(file=file)
        time = h['Header'].attrs['Time'] * iu['utime'] / const['Myr']
        boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / const['kpc']

        # Gaseous component:
        quantity = ['mass', 'temp', 'HI', 'H2', 'HII']
        im_gas_face = get_image_data(file=file, bins=bins, box=box,
                                     quantity=quantity)
        im_gas_edge = get_image_data(file=file, bins=bins, box=box,
                                     quantity=quantity, axis='x', rot=90)

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
            H_stars_face = np.full(shape=(len(bin_sizes_stars)-1,
                                          len(bin_sizes_stars)-1),
                                   fill_value=np.nan)
            H_sf_face = np.full(shape=(len(bin_sizes_sf)-1,
                                       len(bin_sizes_sf)-1), fill_value=np.nan)
            H_stars_edge = np.full(shape=(len(bin_sizes_stars)-1,
                                          len(bin_sizes_stars)-1),
                                   fill_value=np.nan)
            H_sf_edge = np.full(shape=(len(bin_sizes_sf)-1,
                                       len(bin_sizes_sf)-1), fill_value=np.nan)

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
        fig, ax = plt.subplots(4, 4, figsize=(16, 3/4 * 16),
                               gridspec_kw={'height_ratios': [1, 1/2, 1, 1/2],
                                            'width_ratios': [1, 1, 1, 1]})
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
        cax = ax[0, 0].inset_axes([0.05, 0.95, 0.6, 0.015])
        fig.colorbar(im, cax=cax, orientation='horizontal', label=r'''
                     $\log_{10}(\Sigma_\mathrm{Gas} \
                     [\mathrm{g} \ \mathrm{cm}^{-2}])$''')

        im = ax[1, 0].imshow(np.log10(im_gas_edge['mass']),
                             extent=im_gas_edge['extent'], origin='lower',
                             vmin=vmin_dens, vmax=vmax_dens, cmap=cmap_dens)
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
        cax = ax[0, 1].inset_axes([0.05, 0.95, 0.6, 0.015])
        fig.colorbar(im, cax=cax, orientation='horizontal', label=r'''
                     $\log_{10}(T \ [\mathrm{K}])$''')

        im = ax[1, 1].imshow(np.log10(im_gas_edge['temp']),
                             extent=im_gas_edge['extent'], origin='lower',
                             vmin=vmin_temp, vmax=vmax_temp, cmap=cmap_temp)
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
        ax[1, 2].set_aspect('equal')
        ax[1, 2].set_xticks([])
        ax[1, 2].set_yticks([])
        # ax[1, 2].set_ylim(zrange[0], zrange[1])

        im = ax[0, 2].pcolormesh(bin_edges_sf, bin_edges_sf, H_sf_face,
                                 vmin=vmin_sf, vmax=vmax_sf, cmap=cmap_sf)
        cax = ax[0, 2].inset_axes([0.05, 0.95, 0.6, 0.015])
        fig.colorbar(im, cax=cax, orientation='horizontal', label=r'''
                     $\log_{10}(\Sigma_\mathrm{SFR} \ [\mathrm{M}_\odot \
                     \mathrm{yr}^{-1} \ \mathrm{kpc}^{-2}])$''')
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
        ax[1, 3].set_aspect('equal')
        ax[1, 3].set_xticks([])
        ax[1, 3].set_yticks([])
        # ax[1, 3].set_ylim(zrange[0], zrange[1])

        im = ax[0, 3].pcolormesh(bin_edges_stars, bin_edges_stars,
                                 H_stars_face,
                                 vmin=vmin_stars, vmax=vmax_stars,
                                 cmap=cmap_stars)
        cax = ax[0, 3].inset_axes([0.05, 0.95, 0.6, 0.015])
        fig.colorbar(im, cax=cax, orientation='horizontal', label=r'''
                     $\log_{10}(\Sigma_\mathrm{StarP} \ [\mathrm{M}_\odot \
                     \mathrm{kpc}^{-2}])$''')
        ax[1, 3].pcolormesh(bin_edges_stars, bin_edges_stars, H_stars_edge,
                            vmin=vmin_stars, vmax=vmax_stars, cmap=cmap_stars)

        # HI surface density:
        im = ax[2, 0].imshow(np.log10(im_gas_face['HI']),
                             extent=im_gas_face['extent'], origin='lower',
                             vmin=vmin_HI, vmax=vmax_HI, cmap=cmap_HI)
        #ax[2, 0].set_facecolor(cmap_HI(0))
        cax = ax[2, 0].inset_axes([0.05, 0.95, 0.6, 0.015])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(label='$\log_{10}(\Sigma_\mathrm{HI} \ [\mathrm{cm}^{-2}])$', color='w')
        cbar.ax.tick_params(colors='w')
        cbar.outline.set_edgecolor('w')
        ax[2,0].set_aspect('equal')
        ax[2,0].set_xticks([])
        ax[2,0].set_yticks(yticks[1:-1])
        ax[2,0].set_ylabel('$y$ [kpc]')

        ax[3,0].imshow(imHI_edge, extent=(xrange[0], xrange[1], zrange[0], zrange[1]), 
                    vmin=vmin_HI, vmax=vmax_HI, origin='lower', cmap=cmap_H)
        ax[3,0].set_facecolor(cubehelix(0))
        ax[3,0].set_aspect('equal')
        xticks = ax[3,0].get_xticks()
        ax[3,0].set_xticks(xticks[1:-1])
        ax[3,0].set_xlabel('$x$ [kpc]')
        ax[3,0].set_ylabel('$z$ [kpc]')

        # HII surface density:
        im = ax[2,1].imshow(imHII_face, extent=(xrange[0], xrange[1], yrange[0], yrange[1]), 
                            vmin=vmin_HII, vmax=vmax_HII, origin='lower', cmap=cmap_H)
        ax[2,1].set_facecolor(cubehelix(0))
        cax = ax[2,1].inset_axes([0.05, 0.95, 0.6, 0.015])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(label='$\log_{10}(\Sigma_\mathrm{HII} \ [\mathrm{cm}^{-2}])$', color='w')
        cbar.ax.tick_params(colors='w')
        cbar.outline.set_edgecolor('w')
        ax[2,1].set_aspect('equal')
        ax[2,1].set_xticks([])
        ax[2,1].set_yticks([])

        ax[3,1].imshow(imHII_edge, extent=(xrange[0], xrange[1], zrange[0], zrange[1]), 
                    vmin=vmin_HII, vmax=vmax_HII, origin='lower', cmap=cmap_H)
        ax[3,1].set_facecolor(cubehelix(0))
        ax[3,1].set_aspect('equal')
        ax[3,1].set_yticks([])
        ax[3,1].set_xticks(xticks[1:-1])
        ax[3,1].set_xlabel('$x$ [kpc]')

        # H2 surface density:
        im = ax[2,2].imshow(imH2_face, extent=(xrange[0], xrange[1], yrange[0], yrange[1]), 
                            vmin=vmin_H2, vmax=vmax_H2, origin='lower', cmap=cmap_H)
        ax[2,2].set_facecolor(cubehelix(0))
        cax = ax[2,2].inset_axes([0.05, 0.95, 0.6, 0.015])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(label='$\log_{10}(\Sigma_{\mathrm{H}_2} \ [\mathrm{cm}^{-2}])$', color='w')
        cbar.ax.tick_params(colors='w')
        cbar.outline.set_edgecolor('w')
        ax[2,2].set_aspect('equal')
        ax[2,2].set_xticks([])
        ax[2,2].set_yticks([])

        ax[3,2].imshow(imH2_edge, extent=(xrange[0], xrange[1], zrange[0], zrange[1]), 
                    vmin=vmin_H2, vmax=vmax_H2, origin='lower', cmap=cmap_H)
        ax[3,2].set_facecolor(cubehelix(0))
        ax[3,2].set_aspect('equal')
        ax[3,2].set_yticks([])
        ax[3,2].set_xticks(xticks[1:-1])
        ax[3,2].set_xlabel('$x$ [kpc]')

        # Stellar component:
        im = ax[2,3].pcolormesh(bins[0], bins[1], imStellar_face, vmin=vmin_stellar, vmax=vmax_stellar, cmap='bone')
        #im = ax[2,3].imshow(imStellar_face, extent=(xrange[0], xrange[1], yrange[0], yrange[1]), 
        #                    vmin=vmin_stellar, vmax=vmax_stellar, origin='lower', cmap='bone')
        ax[2,3].set_facecolor(bone(0))
        ax[2,3].set_aspect('equal')
        ax[2,3].set_xticks([])
        ax[2,3].set_yticks([])
        cax = ax[2,3].inset_axes([0.05, 0.95, 0.6, 0.015])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(label='$\log_{10}(\Sigma_\star \ [\mathrm{M}_\odot \ \mathrm{kpc}^{-2}])$', color='w')
        cbar.ax.tick_params(colors='w')
        cbar.outline.set_edgecolor('w')
        
        ax[3,3].pcolormesh(bins[0], bins[2], imStellar_edge, vmin=vmin_stellar, vmax=vmax_stellar, cmap='bone')
        #ax[3,3].imshow(imStellar_edge, extent=(xrange[0], xrange[1], zrange[0], zrange[1]), 
        #               vmin=vmin_stellar, vmax=vmax_stellar, origin='lower', cmap='bone')
        ax[3,3].set_facecolor(bone(0))
        ax[3,3].set_aspect('equal')
        ax[3,3].set_yticks([])
        ax[3,3].set_xticks(xticks[1:-1])
        ax[3,3].set_xlabel('$x$ [kpc]')

        # Time:
        ax[0,3].text(1, 1, f'{time:.2f} Myr', fontsize=12, color='k', ha='right', va='bottom', transform=ax[0,3].transAxes)

        fig.savefig(f'./vframes/mosaic/mosaic_{snap}.png')
        plt.close()

    return None

# -------------- End of file
