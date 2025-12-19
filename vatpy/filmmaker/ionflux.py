# Description: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-12-18


# -------------- Required packages
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from ..read import read_hdf5, read_dump
from ..constants import const
from ..get_data.get_image_data import get_image_data

import configv


# -------------- Declare function(s)
def plot_ionflux(start, end):
    '''TODO
    '''
    plt.style.use(f'{configv.homedir}/vatpy/mpl/mosaic.mplstyle')

    # Parameters:
    bins = 100
    box1 = (-2, 2)
    box2 = (-250, 250)
    box3 = (-12.5, 12.5)
    vlen = 1.5

    vmin_dens_lvl1, vmax_dens_lvl1 = -5, -2
    vmin_dens_lvl2, vmax_dens_lvl2 = -5, -2
    vmin_dens_lvl3, vmax_dens_lvl3 = -5, -2

    vmin_temp_lvl1, vmax_temp_lvl1 = 2, 6
    vmin_temp_lvl2, vmax_temp_lvl2 = 2, 6
    vmin_temp_lvl3, vmax_temp_lvl3 = 2, 6

    vmin_flux_lvl1, vmax_flux_lvl1 = -18, -12
    vmin_flux_lvl2, vmax_flux_lvl2 = -18, -12
    vmin_flux_lvl3, vmax_flux_lvl3 = -18, -12

    cmap_dens = 'inferno'
    cmap_temp = 'cmr.amber'
    cmap_flux = 'viridis'

    bh_acc_radius = 1
    label1 = '4x4 kpc'
    label2 = '500x500 pc'
    label3 = '25x25 pc'

    # Loop over snapshots:
    for i in range(start, end+1):
        snap = '000'[:3-len(str(i))] + str(i)
        file = 'snap_' + snap + '.hdf5'
        print(f'  * Generating frame for snapshot {i}')
        h, iu = read_hdf5(file=file)
        time = h['Header'].attrs['Time'] * iu['utime'] / const['Myr']
        boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / const['pc']
        bh_pos = h['PartType5']['Coordinates'][0] * iu['ulength'] / const['pc']
        bh_pos -= boxsize / 2

        # Angular momentum vector:
        dump = f'sink_snap_{snap}'
        if os.path.isfile(dump):
            sinks = read_dump(dump, spin=True, bh=True, hm=True, rcirc=True,
                              rad=True, nfreq=5)[2]
            angmom = sinks['AngularMomentum'][0] * iu['uangmom']
            dirvec = vlen * angmom / np.linalg.norm(angmom)

        # Gaseous component:
        quantity = ['mass', 'temp', 'F056', 'F112', 'F136', 'F152', 'F246']
        im_gas_face1 = get_image_data(file=file, bins=bins, box=box1,
                                      quantity=quantity, halfboxfocus=True)
        print('    - Image 1/6', end='\r')
        im_gas_edge1 = get_image_data(file=file, bins=bins, box=box1,
                                      quantity=quantity, axis='x', rotate=90,
                                      halfboxfocus=True)
        print('    - Image 2/6', end='\r')
        im_gas_face2 = get_image_data(file=file, bins=bins, box=box2,
                                      quantity=quantity, unitlength='pc',
                                      halfboxfocus=True)
        print('    - Image 3/6', end='\r')
        im_gas_edge2 = get_image_data(file=file, bins=bins, box=box2,
                                      quantity=quantity, unitlength='pc',
                                      axis='x', rotate=90, halfboxfocus=True)
        print('    - Image 4/6', end='\r')
        im_gas_face3 = get_image_data(file=file, bins=bins, box=box3,
                                      quantity=quantity, unitlength='pc',
                                      blackholefocus=True)
        print('    - Image 5/6', end='\r')
        im_gas_edge3 = get_image_data(file=file, bins=bins, box=box3,
                                      quantity=quantity, unitlength='pc',
                                      axis='x', rotate=90, blackholefocus=True)
        print('    - Image 6/6', end='\r')

        label_dens = (r'$\log_{10}(\Sigma_\mathrm{Gas} \ [\mathrm{g}$' +
                      r'$ \ \mathrm{cm}^{-2}])$')
        label_temp = r'$\log_{10}(T \ [\mathrm{K}])$'
        label_flux = r'$\log_{10}(R \ [\mathrm{s}^{-1} \ \mathrm{cm}^{-3}])$'

        def make_subplot_panel(im_gas_face, im_gas_edge, ax_face, ax_edge,
                               quantity, vmin, vmax, cmap, label,
                               add_moving_bh=False, add_fixed_bh=False,
                               add_acc_radius=False, add_dirvec=False):
            # Face-on:
            im = ax_face.imshow(np.log10(im_gas_face[quantity]),
                                vmin=vmin, vmax=vmax,
                                extent=im_gas_face['extent'],
                                origin='lower', cmap=cmap)
            cax = ax_face.inset_axes([0.1, 0.88, 0.6, 0.05])
            cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
            cbar.set_label(label=label, color='w', size='small')
            cbar.ax.tick_params(colors='w', which='both', labelsize='small')
            cbar.outline.set_edgecolor('w')
            ax_face.set_xticks([])
            ax_face.set_yticks([])
            ax_face.set_facecolor(plt.get_cmap(cmap)(0))

            # Edge-on:
            im = ax_edge.imshow(np.log10(im_gas_edge[quantity]),
                                vmin=vmin, vmax=vmax,
                                extent=im_gas_edge['extent'],
                                origin='lower', cmap=cmap)
            ax_edge.set_xticks([])
            ax_edge.set_yticks([])
            ax_edge.set_ylim(im_gas_edge['extent'][1]/2,
                             im_gas_edge['extent'][2]/2)
            ax_edge.set_facecolor(plt.get_cmap(cmap)(0))

            if add_moving_bh is True:
                ax_face.plot(bh_pos[0], bh_pos[1], c='k', marker='o', ms=2)
                ax_edge.plot(bh_pos[0], bh_pos[2], c='k', marker='o', ms=2)

            if add_fixed_bh is True:
                for a in [ax_face, ax_edge]:
                    a.plot([0], [0], c='k', marker='o', ms=5)
                    a.plot([0], [0], c='k', marker='o', ms=5,
                           markeredgecolor='w', markeredgewidth=0.5,
                           ls='none')
                    if add_acc_radius is True:
                        a.plot(bh_acc_radius * np.cos(np.linspace(0, 2*np.pi,
                                                                  360)),
                               bh_acc_radius * np.sin(np.linspace(0, 2*np.pi,
                                                                  360)),
                               ls='-', c='tab:cyan', lw=4)
                        a.plot(bh_acc_radius * np.cos(np.linspace(0, 2*np.pi,
                                                                  360)),
                               bh_acc_radius * np.sin(np.linspace(0, 2*np.pi,
                                                                  360)),
                               ls='--', c='k', lw=2)

            if os.path.isfile(dump):
                if add_dirvec is True:
                    arrow = patches.FancyArrow(x=0, y=0, dx=dirvec[0],
                                               dy=dirvec[1], width=0.04,
                                               color='k', head_width=0.3,
                                               head_length=0.3)
                    ax_face.add_patch(arrow)
                    arrow = patches.FancyArrow(x=0, y=0, dx=dirvec[0],
                                               dy=dirvec[2], width=0.03,
                                               color='k', head_width=0.3,
                                               head_length=0.3)
                    ax_edge.add_patch(arrow)

            return None

        # Figure:
        fig, ax = plt.subplots(6, 7, figsize=(15, 15*(4.5/7)),
                               gridspec_kw={'height_ratios': [1, 1/2, 1, 1/2,
                                                              1, 1/2],
                                            'width_ratios': [1, 1, 1, 1, 1,
                                                             1, 1]})
        fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95,
                            wspace=0, hspace=0)

        # Density level 1:
        make_subplot_panel(im_gas_face1, im_gas_edge1, ax[0, 0], ax[1, 0],
                           'mass', vmin_dens_lvl1, vmax_dens_lvl1, cmap_dens,
                           label_dens)

        # Density level 2:
        make_subplot_panel(im_gas_face2, im_gas_edge2, ax[2, 0], ax[3, 0],
                           'mass', vmin_dens_lvl2, vmax_dens_lvl2, cmap_dens,
                           label_dens, add_moving_bh=True)

        # Density level 3:
        make_subplot_panel(im_gas_face3, im_gas_edge3, ax[4, 0], ax[5, 0],
                           'mass', vmin_dens_lvl3, vmax_dens_lvl3, cmap_dens,
                           label_dens, add_fixed_bh=True, add_dirvec=True)

        # Temperature level 1:
        make_subplot_panel(im_gas_face1, im_gas_edge1, ax[0, 1], ax[1, 1],
                           'temp', vmin_temp_lvl1, vmax_temp_lvl1, cmap_temp,
                           label_temp)

        # Temperature level 2:
        make_subplot_panel(im_gas_face2, im_gas_edge2, ax[2, 1], ax[3, 1],
                           'temp', vmin_temp_lvl2, vmax_temp_lvl2, cmap_temp,
                           label_temp, add_moving_bh=True)

        # Temperature level 3:
        make_subplot_panel(im_gas_face3, im_gas_edge3, ax[4, 1], ax[5, 1],
                           'temp', vmin_temp_lvl3, vmax_temp_lvl3, cmap_temp,
                           label_temp, add_fixed_bh=True, add_dirvec=True)

        # F056 level 1:
        make_subplot_panel(im_gas_face1, im_gas_edge1, ax[0, 2], ax[1, 2],
                           'F056', vmin_flux_lvl1, vmax_flux_lvl1,
                           cmap_flux, label_flux)

        # F056 level 2:
        make_subplot_panel(im_gas_face2, im_gas_edge2, ax[2, 2], ax[3, 2],
                           'F056', vmin_flux_lvl2, vmax_flux_lvl2,
                           cmap_flux, label_flux, add_moving_bh=True)

        # F056 level 3:
        make_subplot_panel(im_gas_face3, im_gas_edge3, ax[4, 2], ax[5, 2],
                           'F056', vmin_flux_lvl3, vmax_flux_lvl3,
                           cmap_flux, label_flux, add_fixed_bh=True,
                           add_dirvec=True)

        # F112 level 1:
        make_subplot_panel(im_gas_face1, im_gas_edge1, ax[0, 3], ax[1, 3],
                           'F112', vmin_flux_lvl1, vmax_flux_lvl1,
                           cmap_flux, label_flux)

        # F122 level 2:
        make_subplot_panel(im_gas_face2, im_gas_edge2, ax[2, 3], ax[3, 3],
                           'F112', vmin_flux_lvl2, vmax_flux_lvl2,
                           cmap_flux, label_flux, add_moving_bh=True)

        # F112 level 3:
        make_subplot_panel(im_gas_face3, im_gas_edge3, ax[4, 3], ax[5, 3],
                           'F112', vmin_flux_lvl3, vmax_flux_lvl3,
                           cmap_flux, label_flux, add_fixed_bh=True,
                           add_dirvec=True)

        # F136 level 1:
        make_subplot_panel(im_gas_face1, im_gas_edge1, ax[0, 4], ax[1, 4],
                           'F136', vmin_flux_lvl1, vmax_flux_lvl1,
                           cmap_flux, label_flux)

        # F136 level 2:
        make_subplot_panel(im_gas_face2, im_gas_edge2, ax[2, 4], ax[3, 4],
                           'F136', vmin_flux_lvl2, vmax_flux_lvl2,
                           cmap_flux, label_flux, add_moving_bh=True)

        # F136 level 3:
        make_subplot_panel(im_gas_face3, im_gas_edge3, ax[4, 4], ax[5, 4],
                           'F136', vmin_flux_lvl3, vmax_flux_lvl3,
                           cmap_flux, label_flux, add_fixed_bh=True,
                           add_dirvec=True)

        # F152 level 1:
        make_subplot_panel(im_gas_face1, im_gas_edge1, ax[0, 5], ax[1, 5],
                           'F152', vmin_flux_lvl1, vmax_flux_lvl1,
                           cmap_flux, label_flux)

        # F152 level 2:
        make_subplot_panel(im_gas_face2, im_gas_edge2, ax[2, 5], ax[3, 5],
                           'F152', vmin_flux_lvl2, vmax_flux_lvl2,
                           cmap_flux, label_flux, add_moving_bh=True)

        # F152 level 3:
        make_subplot_panel(im_gas_face3, im_gas_edge3, ax[4, 5], ax[5, 5],
                           'F152', vmin_flux_lvl3, vmax_flux_lvl3,
                           cmap_flux, label_flux, add_fixed_bh=True,
                           add_dirvec=True)

        # F246 level 1:
        make_subplot_panel(im_gas_face1, im_gas_edge1, ax[0, 6], ax[1, 6],
                           'F246', vmin_flux_lvl1, vmax_flux_lvl1,
                           cmap_flux, label_flux)

        # F246 level 2:
        make_subplot_panel(im_gas_face2, im_gas_edge2, ax[2, 6], ax[3, 6],
                           'F246', vmin_flux_lvl2, vmax_flux_lvl2,
                           cmap_flux, label_flux, add_moving_bh=True)

        # F246 level 3:
        make_subplot_panel(im_gas_face3, im_gas_edge3, ax[4, 6], ax[5, 6],
                           'F246', vmin_flux_lvl3, vmax_flux_lvl3,
                           cmap_flux, label_flux, add_fixed_bh=True,
                           add_dirvec=True)

        # Time:
        ax[0, 6].text(1, 1, f'{time:.2f} Myr', fontsize='medium', color='k',
                      ha='right', va='bottom', transform=ax[0, 6].transAxes)

        # Info:
        ax[0, 0].set_ylabel(label1)
        ax[2, 0].set_ylabel(label2)
        ax[4, 0].set_ylabel(label3)

        # Freq. bins:
        ax[0, 2].text(0.95, 0.05, '5.6+ eV', fontsize='small', color='w',
                      ha='right', va='bottom', transform=ax[0, 2].transAxes)
        ax[0, 3].text(0.95, 0.05, '11.2+ eV', fontsize='small', color='w',
                      ha='right', va='bottom', transform=ax[0, 3].transAxes)
        ax[0, 4].text(0.95, 0.05, '13.6+ eV', fontsize='small', color='w',
                      ha='right', va='bottom', transform=ax[0, 4].transAxes)
        ax[0, 5].text(0.95, 0.05, '15.2+ eV', fontsize='small', color='w',
                      ha='right', va='bottom', transform=ax[0, 5].transAxes)
        ax[0, 6].text(0.95, 0.05, '24.6+ eV', fontsize='small', color='w',
                      ha='right', va='bottom', transform=ax[0, 6].transAxes)

        fig.savefig(f'./vframes/ionflux/ionflux_{snap}.png')
        plt.close()

    return None

# -------------- End of file
