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
def plot_zoom(start, end):
    '''TODO
    '''
    plt.style.use(f'{configv.homedir}/vatpy/mpl/mosaic.mplstyle')

    # Parameters:
    bins = 200
    box1 = (-2, 2)
    box2 = (-250, 250)
    box3 = (-12.5, 12.5)

    vmin_dens_lvl1, vmax_dens_lvl1 = -5, -2
    vmin_dens_lvl2, vmax_dens_lvl2 = -5, -2
    vmin_dens_lvl3, vmax_dens_lvl3 = -5, -2

    vmin_temp_lvl1, vmax_temp_lvl1 = 2, 6
    vmin_temp_lvl2, vmax_temp_lvl2 = 2, 6
    vmin_temp_lvl3, vmax_temp_lvl3 = 2, 6

    cmap_dens = 'inferno'
    cmap_temp = 'cmr.amber'

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
        quantity = ['mass', 'temp']
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

        # Figure:
        fig, ax = plt.subplots(4, 3, figsize=(16, 15),
                               gridspec_kw={'height_ratios': [1, 1/2, 1, 1/2],
                                            'width_ratios': [1, 1, 1]})
        fig.subplots_adjust(left=0.02, right=0.98, bottom=0.05, top=0.95,
                            wspace=0, hspace=0)

        # Density level 1 (face-on):
        im = ax[0, 0].imshow(np.log10(im_gas_face1['mass']),
                             vmin=vmin_dens_lvl1, vmax=vmax_dens_lvl1,
                             extent=im_gas_face1['extent'],
                             origin='lower', cmap=cmap_dens)
        cax = ax[0, 0].inset_axes([0.05, 0.95, 0.6, 0.015])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(label=r'$\log_{10}(\Sigma_\mathrm{Gas} \ [\mathrm{g}$' +
                       r'$ \ \mathrm{cm}^{-2}])$', color='k')
        cbar.ax.tick_params(colors='k')
        cbar.outline.set_edgecolor('k')
        ax[0, 0].set_xticks([])
        ax[0, 0].set_yticks([])

        # Density level 1 (edge-on):
        im = ax[1, 0].imshow(np.log10(im_gas_edge1['mass']),
                             vmin=vmin_dens_lvl1, vmax=vmax_dens_lvl1,
                             extent=im_gas_edge1['extent'],
                             origin='lower', cmap=cmap_dens)
        ax[1, 0].set_xticks([])
        ax[1, 0].set_yticks([])
        ax[1, 0].set_ylim(im_gas_edge1['extent'][1]/2,
                          im_gas_edge1['extent'][2]/2)

        # Density level 2 (face-on):
        im = ax[0, 1].imshow(np.log10(im_gas_face2['mass']),
                             vmin=vmin_dens_lvl2, vmax=vmax_dens_lvl2,
                             extent=im_gas_face2['extent'],
                             origin='lower', cmap=cmap_dens)
        cax = ax[0, 1].inset_axes([0.05, 0.95, 0.6, 0.015])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(label=r'$\log_{10}(\Sigma_\mathrm{Gas} \ [\mathrm{g}$' +
                       r'$ \ \mathrm{cm}^{-2}])$', color='k')
        cbar.ax.tick_params(colors='k')
        cbar.outline.set_edgecolor('k')
        ax[0, 1].set_xticks([])
        ax[0, 1].set_yticks([])

        # Density level 2 (edge-on):
        im = ax[1, 1].imshow(np.log10(im_gas_edge2['mass']),
                             vmin=vmin_dens_lvl2, vmax=vmax_dens_lvl2,
                             extent=im_gas_edge2['extent'],
                             origin='lower', cmap=cmap_dens)
        ax[1, 1].set_xticks([])
        ax[1, 1].set_yticks([])
        ax[1, 1].set_ylim(im_gas_edge2['extent'][1]/2,
                          im_gas_edge2['extent'][2]/2)

        # Density level 3 (face-on):
        im = ax[0, 2].imshow(np.log10(im_gas_face3['mass']),
                             vmin=vmin_dens_lvl3, vmax=vmax_dens_lvl3,
                             extent=im_gas_face3['extent'],
                             origin='lower', cmap=cmap_dens)
        cax = ax[0, 2].inset_axes([0.05, 0.95, 0.6, 0.015])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(label=r'$\log_{10}(\Sigma_\mathrm{Gas} \ [\mathrm{g}$' +
                       r'$ \ \mathrm{cm}^{-2}])$', color='k')
        cbar.ax.tick_params(colors='k')
        cbar.outline.set_edgecolor('k')
        ax[0, 2].set_xticks([])
        ax[0, 2].set_yticks([])

        ax[0, 2].plot([0], [0], c='k', marker='o', ms=10,
                      markeredgecolor='tab:cyan', markeredgewidth=2, ls='none')
        ax[0, 2].plot(bh_acc_radius * np.cos(np.linspace(0, 2*np.pi, 360)),
                      bh_acc_radius * np.sin(np.linspace(0, 2*np.pi, 360)),
                      ls='-', c='tab:cyan', lw=4)
        ax[0, 2].plot(bh_acc_radius * np.cos(np.linspace(0, 2*np.pi, 360)),
                      bh_acc_radius * np.sin(np.linspace(0, 2*np.pi, 360)),
                      ls='--', c='k', lw=2)

        # Density level 3 (edge-on):
        im = ax[1, 2].imshow(np.log10(im_gas_edge3['mass']),
                             vmin=vmin_dens_lvl3, vmax=vmax_dens_lvl3,
                             extent=im_gas_edge3['extent'],
                             origin='lower', cmap=cmap_dens)
        ax[1, 2].set_xticks([])
        ax[1, 2].set_yticks([])
        ax[1, 2].set_ylim(im_gas_edge3['extent'][1]/2,
                          im_gas_edge3['extent'][2]/2)

        ax[1, 2].plot([0], [0], c='k', marker='o', ms=10,
                      markeredgecolor='tab:cyan', markeredgewidth=2, ls='none')
        ax[1, 2].plot(bh_acc_radius * np.cos(np.linspace(0, 2*np.pi, 360)),
                      bh_acc_radius * np.sin(np.linspace(0, 2*np.pi, 360)),
                      ls='-', c='tab:cyan', lw=4)
        ax[1, 2].plot(bh_acc_radius * np.cos(np.linspace(0, 2*np.pi, 360)),
                      bh_acc_radius * np.sin(np.linspace(0, 2*np.pi, 360)),
                      ls='--', c='k', lw=2)

        # Temperature level 1 (face-on):
        im = ax[2, 0].imshow(np.log10(im_gas_face1['temp']),
                             vmin=vmin_temp_lvl1, vmax=vmax_temp_lvl1,
                             extent=im_gas_face1['extent'],
                             origin='lower', cmap=cmap_temp)
        cax = ax[2, 0].inset_axes([0.05, 0.95, 0.6, 0.015])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(label=r'$\log_{10}(T \ [\mathrm{K}])$', color='k')
        cbar.ax.tick_params(colors='k')
        cbar.outline.set_edgecolor('k')
        ax[2, 0].set_xticks([])
        ax[2, 0].set_yticks([])

        # Temperature level 1 (edge-on):
        im = ax[3, 0].imshow(np.log10(im_gas_edge1['temp']),
                             vmin=vmin_temp_lvl1, vmax=vmax_temp_lvl1,
                             extent=im_gas_edge1['extent'],
                             origin='lower', cmap=cmap_temp)
        ax[3, 0].set_xticks([])
        ax[3, 0].set_yticks([])
        ax[3, 0].set_ylim(im_gas_edge1['extent'][1]/2,
                          im_gas_edge1['extent'][2]/2)

        # Temperature level 2 (face-on):
        im = ax[2, 1].imshow(np.log10(im_gas_face2['temp']),
                             vmin=vmin_temp_lvl2, vmax=vmax_temp_lvl2,
                             extent=im_gas_face2['extent'],
                             origin='lower', cmap=cmap_temp)
        cax = ax[2, 1].inset_axes([0.05, 0.95, 0.6, 0.015])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(label=r'$\log_{10}(T \ [\mathrm{K}])$', color='k')
        cbar.ax.tick_params(colors='k')
        cbar.outline.set_edgecolor('k')
        ax[2, 1].set_xticks([])
        ax[2, 1].set_yticks([])

        # Temperature level 2 (edge-on):
        im = ax[3, 1].imshow(np.log10(im_gas_edge2['temp']),
                             vmin=vmin_temp_lvl2, vmax=vmax_temp_lvl2,
                             extent=im_gas_edge2['extent'],
                             origin='lower', cmap=cmap_temp)
        ax[3, 1].set_xticks([])
        ax[3, 1].set_yticks([])
        ax[3, 1].set_ylim(im_gas_edge2['extent'][1]/2,
                          im_gas_edge2['extent'][2]/2)

        # Temperature level 3 (face-on):
        im = ax[2, 2].imshow(np.log10(im_gas_face3['temp']),
                             vmin=vmin_temp_lvl3, vmax=vmax_temp_lvl3,
                             extent=im_gas_face3['extent'],
                             origin='lower', cmap=cmap_temp)
        cax = ax[2, 2].inset_axes([0.05, 0.95, 0.6, 0.015])
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_label(label=r'$\log_{10}(T \ [\mathrm{K}])$', color='k')
        cbar.ax.tick_params(colors='k')
        cbar.outline.set_edgecolor('k')
        ax[2, 2].set_xticks([])
        ax[2, 2].set_yticks([])

        ax[2, 2].plot([0], [0], c='k', marker='o', ms=10,
                      markeredgecolor='tab:cyan', markeredgewidth=2, ls='none')
        ax[2, 2].plot(bh_acc_radius * np.cos(np.linspace(0, 2*np.pi, 360)),
                      bh_acc_radius * np.sin(np.linspace(0, 2*np.pi, 360)),
                      ls='-', c='tab:cyan', lw=4)
        ax[2, 2].plot(bh_acc_radius * np.cos(np.linspace(0, 2*np.pi, 360)),
                      bh_acc_radius * np.sin(np.linspace(0, 2*np.pi, 360)),
                      ls='--', c='k', lw=2)

        # Temperature level 3 (edge-on):
        im = ax[3, 2].imshow(np.log10(im_gas_edge3['temp']),
                             vmin=vmin_temp_lvl3, vmax=vmax_temp_lvl3,
                             extent=im_gas_edge3['extent'],
                             origin='lower', cmap=cmap_temp)
        ax[3, 2].set_xticks([])
        ax[3, 2].set_yticks([])
        ax[3, 2].set_ylim(im_gas_edge3['extent'][1]/2,
                          im_gas_edge3['extent'][2]/2)

        ax[3, 2].plot([0], [0], c='k', marker='o', ms=10,
                      markeredgecolor='tab:cyan', markeredgewidth=2, ls='none')
        ax[3, 2].plot(bh_acc_radius * np.cos(np.linspace(0, 2*np.pi, 360)),
                      bh_acc_radius * np.sin(np.linspace(0, 2*np.pi, 360)),
                      ls='-', c='tab:cyan', lw=4)
        ax[3, 2].plot(bh_acc_radius * np.cos(np.linspace(0, 2*np.pi, 360)),
                      bh_acc_radius * np.sin(np.linspace(0, 2*np.pi, 360)),
                      ls='--', c='k', lw=2)

        # Time:
        ax[3, 2].text(1, -0.02, f'{time:.2f} Myr', fontsize=18, color='k',
                      ha='right', va='top', transform=ax[3, 2].transAxes)

        # Info:
        ax[0, 0].text(0, 1.02, label1, fontsize=18, color='k',
                      ha='left', va='bottom', transform=ax[0, 0].transAxes)
        ax[0, 1].text(0, 1.02, label2, fontsize=18, color='k',
                      ha='left', va='bottom', transform=ax[0, 1].transAxes)
        ax[0, 2].text(0, 1.02, label3, fontsize=18, color='k',
                      ha='left', va='bottom', transform=ax[0, 2].transAxes)

        fig.savefig(f'./vframes/zoom/zoom_{snap}.png')
        plt.close()

    return None

# -------------- End of file
