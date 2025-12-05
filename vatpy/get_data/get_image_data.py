'''
Description: Functions to get multiple images at the same time.

Last updated: 2023-09-27
'''

# -------------- Required packages
import numpy as np
from scipy.spatial.transform import Rotation

from ..constants import const
from ..read import read_hdf5
from ..get_gas_property import number_density, temperature
from ..interpolation import interpolate_to_2d_kdtree, interpolate_to_2d


# -------------- Declare function(s)
def get_image_data(file, axis='z', rotate=0, quantity=['mass'], bins=100,
                   interpolation='kdtree', unitlength='kpc',
                   blackholefocus=False, halfboxfocus=False, xrange=None,
                   yrange=None, zrange=None, box=None, cut=None,
                   xraybins=0):
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
    mass = h['PartType0']['Masses'] * iu['umass']
    dens_gas = h['PartType0']['Density'] * iu['udens']
    num = number_density(h, iu)
    temp = temperature(h, iu)
    vol = mass / dens_gas

    if 'MagneticField' in h['PartType0'].keys():
        bfield = h['PartType0']['MagneticField'] * iu['umagfield']
        bfield_tot = np.linalg.norm(bfield, axis=1)
    else:
        bfield_tot = np.zeros(shape=np.shape(dens_gas))

    photonrates = {}
    if 'PhotonRates' in h['PartType0'].keys():
        rates = h['PartType0']['PhotonRates'][:]
        photonrates['HIPI'] = rates[:, 0] / vol
        photonrates['HIPH'] = rates[:, 1] / vol
        photonrates['H2PI'] = rates[:, 2] / vol
        photonrates['H2PH'] = rates[:, 3] / vol
        photonrates['H2PD'] = rates[:, 4] / vol
        photonrates['DUSTPE'] = rates[:, 5] / vol
        photonrates['HEPI'] = rates[:, 6] / vol
        photonrates['HEPH'] = rates[:, 7] / vol

    photonflux = {}
    if 'PhotonFlux' in h['PartType0'].keys():
        flux = h['PartType0']['PhotonFlux'][:]
        photonflux['F056'] = flux[:, 0] / vol
        photonflux['F112'] = flux[:, 1] / vol
        photonflux['F136'] = flux[:, 2] / vol
        photonflux['F152'] = flux[:, 3] / vol
        photonflux['F246'] = flux[:, 4] / vol
        if xraybins >= 1:
            photonflux['FXRAY0'] = flux[:, 5] / vol
        if xraybins >= 2:
            photonflux['FXRAY1'] = flux[:, 6] / vol
        if xraybins >= 4:
            photonflux['FXRAY2'] = flux[:, 7] / vol
            photonflux['FXRAY3'] = flux[:, 8] / vol

    if blackholefocus:
        bh = (h['PartType5']['Coordinates'][0] * iu['ulength'] / ulength)
        pos -= bh

    if halfboxfocus:
        pos -= boxsize / 2

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
    if axis != 'xyz':
        if rotate != 0:
            rotation = Rotation.from_euler(axis, rotate, degrees=True)
            if blackholefocus or halfboxfocus:
                pos = rotation.apply(pos)
            else:
                pos = rotation.apply(pos - boxsize/2)
                pos += boxsize/2
    else:
        rotation = Rotation.from_euler(axis, rotate, degrees=True)
        if blackholefocus or halfboxfocus:
            pos = rotation.apply(pos)
        else:
            pos = rotation.apply(pos - boxsize/2)
            pos += boxsize/2

    # Selection of gas quantity:
    image_data = {}
    for i in range(0, len(quantity)):
        if ((quantity[i] not in ['mass', 'temp', 'bfield']) and
                (quantity[i] not in photonrates.keys()) and
                (quantity[i] not in photonflux.keys())):
            dens = num[quantity[i]]
        else:
            dens = dens_gas

        if quantity[i] == 'temp':
            values = temp
            weights = dens
        elif quantity[i] == 'bfield':
            values = bfield_tot
            weights = dens
        elif quantity[i] in photonrates.keys():
            values = photonrates[quantity[i]]
            weights = dens
        elif quantity[i] in photonflux.keys():
            values = photonflux[quantity[i]]
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
    image_data['time'] = time

    return image_data

# -------------- End of file
