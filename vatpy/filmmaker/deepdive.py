# Description: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-10-07


# -------------- Required packages
import numpy as np
import matplotlib.pyplot as plt

from ..read import read_hdf5
from ..constants import const
from ..get_data.get_image_data import get_image_data

import configv


# -------------- Declare function(s)
def plot_deepdive(start):
    '''TODO
    '''
    plt.style.use(f'{configv.homedir}/vatpy/mpl/noctua.mplstyle')

    bins = 200

    vmin_dens, vmax_dens = -5, -2
    vmin_temp, vmax_temp = 2, 6

    cmap_dens = 'cmr.rainforest'
    cmap_temp = 'cmr.amber'
    cmap_age = 'cmr.guppy_r'

    # Spherical coordinates:
    Nact = 4
    radius_list = []
    rot_list = []

    # Act 1:
    Nframes = 10
    radius_list.append(np.linspace(20, 2, Nframes))
    rot_list.append(np.full((Nframes, 3), 0))

    # Act 2:
    Nframes = 10
    radius_list.append(np.full(Nframes, 2))
    rot_array = np.full((Nframes, 3), 0)
    rot_array[:, 0] = np.linspace(0, 360, Nframes)
    rot_array[:, 2] = np.linspace(0, 360, Nframes)
    rot_list.append(rot_array)

    # Act 3:
    Nframes = 10
    radius_list.append(np.linspace(2, 0.1, Nframes))
    rot_list.append(np.full((Nframes, 3), 0))

    # Act 4:
    Nframes = 20
    radius_list.append(np.full(Nframes, 0.1))
    rot_array = np.full((Nframes, 3), 0)
    xrot = [np.linspace(0, 90, int(Nframes/4)),
            np.linspace(90, 90, int(Nframes/4)),
            np.linspace(90, 90, int(Nframes/4)),
            np.linspace(90, 90, int(Nframes/4))]
    print(xrot)
    rot_array[:, 0] = np.array(xrot).flatten()
    zrot = [np.linspace(0, 0, int(Nframes/4)),
            np.linspace(0, 360, int(Nframes/4)),
            np.linspace(0, 360, int(Nframes/4)),
            np.linspace(0, 360, int(Nframes/4))]
    rot_array[:, 2] = np.array(zrot).flatten()
    rot_list.append(rot_array)

    # Snapshot file:
    snap = '000'[:3-len(str(start))] + str(start)
    file = 'snap_' + snap + '.hdf5'
    print(f'  * Generating frames for snapshot {start}')

    frame_nr = 0
    for i in range(0, Nact):
        for j in range(0, len(radius_list[i])):
            print(f'  * Frame {frame_nr}')

            box = (-radius_list[i][j], radius_list[i][j])
            rotate = rot_list[i][j]

            im_dens = get_image_data(file=file, bins=bins, box=box,
                                     blackholefocus=True, unitlength='kpc',
                                     quantity=['mass'], axis='zyx',
                                     rotate=rotate)

            fig, ax = plt.subplots(figsize=(10, 10), layout='constrained')

            im = ax.imshow(np.log10(im_dens['mass']), cmap=cmap_dens,
                           extent=im_dens['extent'],
                           vmin=vmin_dens, vmax=vmax_dens)

            frame = '000'[:3-len(str(frame_nr))] + str(frame_nr)
            fig.savefig(f'./vframes/deepdive/deepdive_{frame}.png')
            plt.close()

            frame_nr += 1

        '''
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
        '''

    return frame_nr - 1

# -------------- End of file
