# Description: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-09-15


# -------------- Required packages
import numpy as np
import matplotlib.pyplot as plt

from ..read import read_hdf5
from ..constants import const
from ..get_data.get_image_data import get_image_data

import configv


# -------------- Declare function(s)
def plot_noctua(start, end):
    '''TODO
    '''
    plt.style.use(f'{configv.homedir}/vatpy/mpl/noctua.mplstyle')

    bins = 300

    vmin_dens, vmax_dens = -5, -2
    vmin_temp, vmax_temp = 2, 6

    cmap_dens = 'cmr.rainforest'
    cmap_temp = 'cmr.amber'
    cmap_age = 'cmr.guppy_r'

    # Loop over snapshots:
    for i in range(start, end+1):
        snap = '000'[:3-len(str(i))] + str(i)
        file = 'snap_' + snap + '.hdf5'
        print(f'  * Generating frame for snapshot {i}')

        im_dens_global = get_image_data(file=file, bins=bins, box=(48, 52),
                                        quantity=['mass', 'temp'])
        im_dens_local = get_image_data(file=file, bins=bins,
                                       box=(49.75, 50.25),
                                       quantity=['mass', 'temp'])
        im_dens_bh = get_image_data(file=file, bins=bins, box=(-25, 25),
                                    blackholefocus=True, unitlength='pc',
                                    quantity=['mass', 'temp'])

        fig = plt.figure(figsize=(12, 9.6), layout='constrained')
        ax1 = plt.subplot2grid((4, 5), (0, 0), rowspan=3, colspan=3)
        ax2 = plt.subplot2grid((4, 5), (0, 3), rowspan=2, colspan=2)
        ax3 = plt.subplot2grid((4, 5), (2, 3), rowspan=1, colspan=1)
        ax4 = plt.subplot2grid((4, 5), (2, 4), rowspan=1, colspan=1)
        ax5 = plt.subplot2grid((4, 5), (3, 0), rowspan=1, colspan=1)
        ax6 = plt.subplot2grid((4, 5), (3, 1), rowspan=1, colspan=1)
        ax7 = plt.subplot2grid((4, 5), (3, 2), rowspan=1, colspan=1)

        im1 = ax1.imshow(np.log10(im_dens_global['mass']), cmap=cmap_dens,
                         extent=im_dens_global['extent'],
                         vmin=vmin_dens, vmax=vmax_dens)
        ax2.imshow(np.log10(im_dens_local['mass']), cmap=cmap_dens,
                   extent=im_dens_local['extent'],
                   vmin=vmin_dens, vmax=vmax_dens)
        ax3.imshow(np.log10(im_dens_bh['mass']), cmap=cmap_dens,
                   extent=im_dens_bh['extent'], vmin=vmin_dens, vmax=vmax_dens)
        ax4.imshow(np.log10(im_dens_bh['mass']), cmap=cmap_dens,
                   extent=im_dens_bh['extent'], vmin=vmin_dens, vmax=vmax_dens)

        im5 = ax5.imshow(np.log10(im_dens_global['temp']), cmap=cmap_temp,
                         extent=im_dens_global['extent'],
                         vmin=vmin_temp, vmax=vmax_temp)
        ax6.imshow(np.log10(im_dens_local['temp']), cmap=cmap_temp,
                   extent=im_dens_local['extent'],
                   vmin=vmin_temp, vmax=vmax_temp)
        ax7.imshow(np.log10(im_dens_bh['temp']), cmap=cmap_temp,
                   extent=im_dens_bh['extent'], vmin=vmin_temp, vmax=vmax_temp)

        # Colorbars:
        cax1 = ax1.inset_axes([0.05, 0.95, 0.4, 0.03])
        fig.colorbar(im1, cax=cax1, orientation='horizontal',
                     label=r'$\log_{10}(\Sigma \ [\text{g} \ $' +
                           r'$\text{cm}^{-2}])$')

        ax8 = plt.subplot2grid((4, 5), (3, 3), rowspan=1, colspan=1)
        ax8.set_axis_off()
        cax8 = ax8.inset_axes([0, 0, 0.1, 1.0])
        fig.colorbar(im5, cax=cax8, orientation='vertical',
                     label=r'$\log_{10}(T \ [\text{K}])$')

        # Remove x and y ticks, but add scales:
        for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7]:
            ax.set_xticks([])
            ax.set_yticks([])

        ax1.set_title('4x4 kpc', loc='left')
        ax2.set_title('500x500 pc', loc='left')
        ax3.set_title('50x50 pc', loc='left')

        ax5.set_title('4x4 kpc', loc='left')
        ax6.set_title('500x500 pc', loc='left')
        ax7.set_title('50x50 pc', loc='left')

        # Zoom:
        inset_indicator1 = ax1.indicate_inset_zoom(ax2, edgecolor='k', lw=1,
                                                   alpha=1)
        inset_indicator1.connectors[0].set_visible(True)
        inset_indicator1.connectors[1].set_visible(True)
        inset_indicator1.connectors[2].set_visible(False)
        inset_indicator1.connectors[3].set_visible(False)

        inset_indicator5 = ax5.indicate_inset_zoom(ax6, edgecolor='k', lw=1,
                                                   alpha=1)
        inset_indicator5.connectors[0].set_visible(True)
        inset_indicator5.connectors[1].set_visible(True)
        inset_indicator5.connectors[2].set_visible(False)
        inset_indicator5.connectors[3].set_visible(False)

        # BH:
        h, iu = read_hdf5(file=file)
        coord_bh = h['PartType5']['Coordinates'] * iu['ulength'] / const['kpc']
        ax2.scatter(coord_bh[:, 0], coord_bh[:, 1], marker='.', s=140, c='k')
        ax6.scatter(coord_bh[:, 0], coord_bh[:, 1], marker='.', s=40, c='k')

        ax3.scatter([0], [0], marker='o', s=40, c='k', zorder=10)
        ax4.scatter([0], [0], marker='o', s=40, c='k', zorder=10)
        ax7.scatter([0], [0], marker='o', s=40, c='k', zorder=10)

        # NSC:
        if 'PartType3' in h.keys():
            coord_nsc = (h['PartType3']['Coordinates'] * iu['ulength'] /
                         const['pc'])
            coord_nsc -= coord_bh * 1e3

            ax3.scatter(coord_nsc[:, 0], coord_nsc[:, 1], marker='.', s=1,
                        c='w', alpha=0.2)
            ax3.set_xlim(im_dens_bh['extent'][0], im_dens_bh['extent'][1])
            ax3.set_ylim(im_dens_bh['extent'][2], im_dens_bh['extent'][3])

        # Newly formed stars:
        if 'PartType4' in h.keys():
            coord_stars = (h['PartType4']['Coordinates'] * iu['ulength'] /
                           const['pc'])
            coord_stars -= coord_bh * 1e3
            lifetime = ((h['PartType4']['StellarFormationTime'] * iu['utime'] -
                         h['Header'].attrs['Time']) / const['Myr'])

            mask = lifetime < 30
            ax4.scatter(coord_stars[:, 0][mask], coord_stars[:, 1][mask],
                        marker='.', s=5, c=lifetime[mask], cmap=cmap_age,
                        alpha=1.0)

        # Info panel:
        ax9 = plt.subplot2grid((4, 5), (3, 4), rowspan=1, colspan=1)
        ax9.set_axis_off()
        time = h['Header'].attrs['Time'] * iu['utime'] / const['Myr']
        m_bh = h['PartType5']['Masses'][0] * iu['umass'] / const['Msol']
        ax9.text(0.05, 0.95, f'Time: {round(time, 2)} Myr\n' +
                 r'$M_\text{BH}$: ' + f'{round(m_bh)} ' +
                 r'$\text{M}_\odot$',
                 ha='left', va='top', transform=ax9.transAxes)

        fig.savefig(f'./vframes/noctua/noctua_{snap}.png')
        plt.close()

    return None

# -------------- End of file
