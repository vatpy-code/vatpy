# Description: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-11-14


# -------------- Required packages
import numpy as np
import matplotlib.pyplot as plt

from ..read import read_hdf5
from ..constants import const
from ..get_data.get_image_data import get_image_data

import configv


# -------------- Declare function(s)
def plot_ionrate(start, end):
    '''TODO
    '''
    plt.style.use(f'{configv.homedir}/vatpy/mpl/mosaic.mplstyle')

    # Parameters:
    bins = 100
    box1 = (-2, 2)
    box2 = (-250, 250)
    box3 = (-12.5, 12.5)

    vmin_dens_lvl1, vmax_dens_lvl1 = -5, -2
    vmin_dens_lvl2, vmax_dens_lvl2 = -5, -2
    vmin_dens_lvl3, vmax_dens_lvl3 = -5, -2

    vmin_temp_lvl1, vmax_temp_lvl1 = 2, 6
    vmin_temp_lvl2, vmax_temp_lvl2 = 2, 6
    vmin_temp_lvl3, vmax_temp_lvl3 = 2, 6

    vmin_rate_lvl1, vmax_rate_lvl1 = None, None
    vmin_rate_lvl2, vmax_rate_lvl2 = None, None
    vmin_rate_lvl3, vmax_rate_lvl3 = None, None

    cmap_dens = 'inferno'
    cmap_temp = 'cmr.amber'
    cmap_rate = 'viridis'

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

        # Gaseous component:
        quantity = ['mass', 'temp', 'HIPI', 'H2PD']
        im_gas_face1 = get_image_data(file=file, bins=bins, box=box1,
                                      quantity=quantity, halfboxfocus=True)
        im_gas_edge1 = get_image_data(file=file, bins=bins, box=box1,
                                      quantity=quantity, axis='x', rotate=90,
                                      halfboxfocus=True)
        im_gas_face2 = get_image_data(file=file, bins=bins, box=box2,
                                      quantity=quantity, unitlength='pc',
                                      halfboxfocus=True)
        im_gas_edge2 = get_image_data(file=file, bins=bins, box=box2,
                                      quantity=quantity, unitlength='pc',
                                      axis='x', rotate=90, halfboxfocus=True)
        im_gas_face3 = get_image_data(file=file, bins=bins, box=box3,
                                      quantity=quantity, unitlength='pc',
                                      blackholefocus=True)
        im_gas_edge3 = get_image_data(file=file, bins=bins, box=box3,
                                      quantity=quantity, unitlength='pc',
                                      axis='x', rotate=90, blackholefocus=True)

        label_dens = (r'$\log_{10}(\Sigma_\mathrm{Gas} \ [\mathrm{g}$' +
                      r'$ \ \mathrm{cm}^{-2}])$')
        label_temp = r'$\log_{10}(T \ [\mathrm{K}])$'
        label_rate = r'$\log_{10}(R \ [\mathrm{s}^{-1} \ \mathrm{cm}^{-3}])$'

        def make_subplot_panel(im_gas_face, im_gas_edge, ax_face, ax_edge,
                               quantity, vmin, vmax, cmap, label,
                               add_fixed_bh=False):
            # Face-on:
            im = ax_face.imshow(np.log10(im_gas_face[quantity]),
                                vmin=vmin, vmax=vmax,
                                extent=im_gas_face['extent'],
                                origin='lower', cmap=cmap)
            cax = ax_face.inset_axes([0.05, 0.95, 0.6, 0.015])
            cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
            cbar.set_label(label=label, color='k')
            cbar.ax.tick_params(colors='k')
            cbar.outline.set_edgecolor('k')
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

            if add_fixed_bh is True:
                for a in [ax_face, ax_edge]:
                    a.plot([0], [0], c='k', marker='o', ms=10,
                           markeredgecolor='tab:cyan', markeredgewidth=2,
                           ls='none')
                    a.plot(bh_acc_radius*np.cos(np.linspace(0, 2*np.pi, 360)),
                           bh_acc_radius*np.sin(np.linspace(0, 2*np.pi, 360)),
                           ls='-', c='tab:cyan', lw=4)
                    a.plot(bh_acc_radius*np.cos(np.linspace(0, 2*np.pi, 360)),
                           bh_acc_radius*np.sin(np.linspace(0, 2*np.pi, 360)),
                           ls='--', c='k', lw=2)

            return None

        # Figure:
        fig, ax = plt.subplots(6, 4, figsize=(12, 12*(4.5/4)),
                               gridspec_kw={'height_ratios': [1, 1/2, 1, 1/2,
                                                              1, 1/2],
                                            'width_ratios': [1, 1, 1, 1]})
        fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95,
                            wspace=0, hspace=0)

        # Density level 1:
        make_subplot_panel(im_gas_face1, im_gas_edge1, ax[0, 0], ax[1, 0],
                           'mass', vmin_dens_lvl1, vmax_dens_lvl1, cmap_dens,
                           label_dens)

        # Density level 2:
        make_subplot_panel(im_gas_face2, im_gas_edge2, ax[2, 0], ax[3, 0],
                           'mass', vmin_dens_lvl2, vmax_dens_lvl2, cmap_dens,
                           label_dens)

        # Density level 3:
        make_subplot_panel(im_gas_face3, im_gas_edge3, ax[4, 0], ax[5, 0],
                           'mass', vmin_dens_lvl3, vmax_dens_lvl3, cmap_dens,
                           label_dens, add_fixed_bh=True)

        # Temperature level 1:
        make_subplot_panel(im_gas_face1, im_gas_edge1, ax[0, 1], ax[1, 1],
                           'temp', vmin_temp_lvl1, vmax_temp_lvl1, cmap_temp,
                           label_temp)

        # Temperature level 2:
        make_subplot_panel(im_gas_face2, im_gas_edge2, ax[2, 1], ax[3, 1],
                           'temp', vmin_temp_lvl2, vmax_temp_lvl2, cmap_temp,
                           label_temp)

        # Temperature level 3:
        make_subplot_panel(im_gas_face3, im_gas_edge3, ax[4, 1], ax[5, 1],
                           'temp', vmin_temp_lvl3, vmax_temp_lvl3, cmap_temp,
                           label_temp, add_fixed_bh=True)

        # HIPI level 1:
        make_subplot_panel(im_gas_face1, im_gas_edge1, ax[0, 2], ax[1, 2],
                           'HIPI', vmin_rate_lvl1, vmax_rate_lvl1,
                           cmap_rate, label_rate)

        # HIPI level 2:
        make_subplot_panel(im_gas_face2, im_gas_edge2, ax[2, 2], ax[3, 2],
                           'HIPI', vmin_rate_lvl2, vmax_rate_lvl2,
                           cmap_rate, label_rate)

        # HIPI level 3:
        make_subplot_panel(im_gas_face3, im_gas_edge3, ax[4, 2], ax[5, 2],
                           'HIPI', vmin_rate_lvl3, vmax_rate_lvl3,
                           cmap_rate, label_rate, add_fixed_bh=True)

        # H2PD level 1:
        make_subplot_panel(im_gas_face1, im_gas_edge1, ax[0, 3], ax[1, 3],
                           'H2PD', vmin_rate_lvl1, vmax_rate_lvl1,
                           cmap_rate, label_rate)

        # H2PD level 2:
        make_subplot_panel(im_gas_face2, im_gas_edge2, ax[2, 3], ax[3, 3],
                           'H2PD', vmin_rate_lvl2, vmax_rate_lvl2,
                           cmap_rate, label_rate)

        # HIPI level 3:
        make_subplot_panel(im_gas_face3, im_gas_edge3, ax[4, 3], ax[5, 3],
                           'H2PD', vmin_rate_lvl3, vmax_rate_lvl3,
                           cmap_rate, label_rate, add_fixed_bh=True)

        # Time:
        ax[0, 3].text(1, 1, f'{time:.2f} Myr', fontsize=18, color='k',
                      ha='right', va='bottom', transform=ax[0, 3].transAxes)

        # Info:
        ax[0, 0].set_ylabel(label1)
        ax[2, 0].set_ylabel(label2)
        ax[4, 0].set_ylabel(label3)

        # ax[0, 0].text(0, 1.02, label1, fontsize=18, color='k',
        #               ha='left', va='bottom', transform=ax[0, 0].transAxes)
        # ax[0, 1].text(0, 1.02, label2, fontsize=18, color='k',
        #               ha='left', va='bottom', transform=ax[0, 1].transAxes)
        # ax[0, 2].text(0, 1.02, label3, fontsize=18, color='k',
        #               ha='left', va='bottom', transform=ax[0, 2].transAxes)

        fig.savefig(f'./vframes/ionrate/ionrate_{snap}.png')
        plt.close()

    return None

# -------------- End of file
