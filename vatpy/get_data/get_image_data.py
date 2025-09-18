'''
Description: Functions to get multiple images at the same time.

Last updated: 2023-09-27
'''

# -------------- Required packages
import numpy as np

from ..constants import const
from ..read import read_hdf5
from ..get_gas_property import number_density, temperature
from ..interpolation import interpolate_to_2d_kdtree, interpolate_to_2d


# -------------- Declare function(s)
def get_image_data(file, axis='z', rotate=0, quantity=['mass'], bins=100,
                   interpolation='kdtree', unitlength='kpc', blackholefocus=False,
                   xrange=None, yrange=None, zrange=None, box=None, cut=None):
    # Unit length:
    if unitlength == 'kpc':
        ulength = const['kpc']
    elif unitlength == 'pc':
        ulength = const['pc']
    else:
        ulength = 1

    # Read the HDF5-file:
    h, iu = read_hdf5(file=file)
    time = h['Header'].attrs['Time'] * iu['utime'] / const['Myr']
    boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / ulength
    pos = h['PartType0']['Coordinates'] * iu['ulength'] / ulength
    dens_gas = h['PartType0']['Density'] * iu['udens']
    num = number_density(h, iu)
    temp = temperature(h, iu)

    if blackholefocus:
        bh = (h['PartType5']['Coordinates'][0] * iu['ulength'] / ulength)
        pos -= bh

    # Coordinate ranges:
    if not box:
        if not xrange:
            if blackholefocus:
                xrange = (-boxsize/2, boxsize/2)
            else:
                xrange = (0, boxsize)
        if not yrange:
            if blackholefocus:
                yrange = (-boxsize/2, boxsize/2)
            else:
                yrange = (0, boxsize)
        if not zrange:
            if blackholefocus:
                zrange = (-boxsize/2, boxsize/2)
            else:
                zrange = (0, boxsize)
    else:
        xrange = (box[0], box[1])
        yrange = (box[0], box[1])
        zrange = (box[0], box[1])

    # Rotation:
    if rotate != 0:
        rotation = Rotation.from_euler(axis, rotate, degrees=True)
        if blackholefocus:
            pos = rotation.apply(pos)
        else:
            pos = rotation.apply(pos - boxsize/2)
            pos += boxsize/2

    # Selection of gas quantity:
    image_data = {}
    for i in range(0, len(quantity)):
        if ((quantity[i] != 'mass') and (quantity[i] != 'temp')):
            dens = num[quantity[i]]
        else:
            dens = dens_gas

        if quantity[i] == 'temp':
            values = temp
            weights = dens
        else:
            values = dens
            weights = None

        # Interpolation:
        if interpolation == 'kdtree':
            interpValues = interpolate_to_2d_kdtree(pos=pos, unit=ulength,
                                                    values=values, bins=bins,
                                                    xrange=xrange,
                                                    yrange=yrange,
                                                    zrange=zrange, cut=cut,
                                                    weights=weights)

        else:
            interpValues = interpolate_to_2d(pos=pos, unit=ulength,
                                             values=values, bins=bins,
                                             xrange=xrange, yrange=yrange,
                                             zrange=zrange, cut=cut,
                                             weights=weights)
        # Add the image data:
        image_data[quantity[i]] = interpValues

    image_data['extent'] = [xrange[0], xrange[1], yrange[0], yrange[1]]

    return image_data

# -------------- End of file
