# Description: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-10-28


# -------------- Required packages
import numpy as np

from ..constants import const
from ..read import read_hdf5


# -------------- Declare function(s)
def get_surface_density_data(output_dir, n, N, Rmin, Rmax, Zcut, bins,
                             RmaxHalo=None, ZmaxHalo=None):
    '''TODO
    '''
    # Directory of surface densities:
    sigma_data = {}

    # Radius & area bins:
    r_bins = np.logspace(np.log10(Rmin), np.log10(Rmax), bins)
    a_bins = np.pi * (r_bins[1:]**2 - r_bins[:-1]**2)

    if RmaxHalo:
        r_bins_dm = np.logspace(np.log10(Rmin), np.log10(RmaxHalo), bins)
        a_bins_dm = np.pi * (r_bins_dm[1:]**2 - r_bins_dm[:-1]**2)
    else:
        r_bins_dm = np.logspace(np.log10(Rmin), np.log10(RmaxHalo), bins)
        a_bins_dm = np.pi * (r_bins_dm[1:]**2 - r_bins_dm[:-1]**2)

    sigma0 = []
    sigma1 = []
    sigma2 = []
    sigma3 = []
    sigma4 = []
    for i in range(n, N+1):
        snap = '000'[len(str(i)):] + str(i)
        print(f'Reading data of snapshot {snap}', end='\r')

        h, iu = read_hdf5(f'{output_dir}/snap_{snap}.hdf5')
        # boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / const['pc']

        pos_gas = h['PartType0']['Coordinates'] * iu['ulength'] / const['pc']
        mass_gas = h['PartType0']['Masses'] * iu['umass'] / const['Msol']
        pos_dm = h['PartType1']['Coordinates'] * iu['ulength'] / const['pc']
        mass_dm = h['PartType1']['Masses'] * iu['umass'] / const['Msol']

        if 'PartType2' in h.keys():
            pos_disc = (h['PartType2']['Coordinates'] * iu['ulength'] /
                        const['pc'])
            mass_disc = h['PartType2']['Masses'] * iu['umass'] / const['Msol']

        if 'PartType3' in h.keys():
            pos_nsc = (h['PartType3']['Coordinates'] * iu['ulength'] /
                       const['pc'])
            mass_nsc = h['PartType3']['Masses'] * iu['umass'] / const['Msol']

        if 'PartType4' in h.keys():
            pos_stars = (h['PartType4']['Coordinates'] * iu['ulength'] /
                         const['pc'])
            mass_stars = h['PartType4']['Masses'] * iu['umass'] / const['Msol']

        pos_bh = h['PartType5']['Coordinates'][0] * iu['ulength'] / const['pc']

        # Centre the data on the BH:
        pos_gas -= pos_bh
        pos_dm -= pos_bh
        if 'PartType2' in h.keys():
            pos_disc -= pos_bh
        if 'PartType3' in h.keys():
            pos_nsc -= pos_bh
        if 'PartType4' in h.keys():
            pos_stars -= pos_bh

        # Limit the data in Z:
        if Zcut:
            mass_gas = mass_gas[np.abs(pos_gas[:, 2]) < Zcut]
            pos_gas = pos_gas[np.abs(pos_gas[:, 2]) < Zcut]
            if ZmaxHalo:
                mass_dm = mass_dm[np.abs(pos_dm[:, 2]) < ZmaxHalo]
                pos_dm = pos_dm[np.abs(pos_dm[:, 2]) < ZmaxHalo]
            else:
                mass_dm = mass_dm[np.abs(pos_dm[:, 2]) < Zcut]
                pos_dm = pos_dm[np.abs(pos_dm[:, 2]) < Zcut]
            if 'PartType2' in h.keys():
                mass_disc = mass_disc[np.abs(pos_disc[:, 2]) < Zcut]
                pos_disc = pos_disc[np.abs(pos_disc[:, 2]) < Zcut]
            if 'PartType3' in h.keys():
                mass_nsc = mass_nsc[np.abs(pos_nsc[:, 2]) < Zcut]
                pos_nsc = pos_nsc[np.abs(pos_nsc[:, 2]) < Zcut]
            if 'PartType4' in h.keys():
                mass_stars = mass_stars[np.abs(pos_stars[:, 2]) < Zcut]
                pos_stars = pos_stars[np.abs(pos_stars[:, 2]) < Zcut]

        '''
        # Limit the data in R:
        r_gas = np.linalg.norm(pos_gas[:, :2], axis=1)
        mass_gas = mass_gas[r_gas < R]
        r_gas = r_gas[r_gas < R]
        r_dm = np.linalg.norm(pos_dm[:, :2], axis=1)
        if HaloRmax:
            mass_dm = mass_dm[r_dm < HaloRmax]
            r_dm = r_dm[r_dm < HaloRmax]
        else:
            mass_dm = mass_dm[r_dm < R]
            r_dm = r_dm[r_dm < R]
        if 'PartType2' in h.keys():
            r_disc = np.linalg.norm(pos_disc[:, :2], axis=1)
            mass_disc = mass_disc[r_disc < R]
            r_disc = r_disc[r_disc < R]
        if 'PartType3' in h.keys():
            r_nsc = np.linalg.norm(pos_nsc[:, :2], axis=1)
            mass_nsc = mass_nsc[r_nsc < R]
            r_nsc = r_nsc[r_nsc < R]
        if 'PartType4' in h.keys():
            r_stars = np.linalg.norm(pos_stars[:, :2], axis=1)
            mass_stars = mass_stars[r_stars < R]
            r_stars = r_stars[r_stars < R]
        '''

        r_gas = np.linalg.norm(pos_gas[:, :2], axis=1)
        r_dm = np.linalg.norm(pos_dm[:, :2], axis=1)
        if 'PartType2' in h.keys():
            r_disc = np.linalg.norm(pos_disc[:, :2], axis=1)
        if 'PartType3' in h.keys():
            r_nsc = np.linalg.norm(pos_nsc[:, :2], axis=1)
        if 'PartType4' in h.keys():
            r_stars = np.linalg.norm(pos_stars[:, :2], axis=1)

        # Surface densities:
        H_gas = np.histogram(r_gas, bins=r_bins, weights=mass_gas)[0] / a_bins
        sigma0.append(H_gas)
        H_dm = (np.histogram(r_dm, bins=r_bins_dm, weights=mass_dm)[0] /
                a_bins_dm)
        sigma1.append(H_dm)
        if 'PartType2' in h.keys():
            H_disc = (np.histogram(r_disc, bins=r_bins, weights=mass_disc)[0] /
                      a_bins)
            sigma2.append(H_disc)
        if 'PartType3' in h.keys():
            H_nsc = (np.histogram(r_nsc, bins=r_bins, weights=mass_nsc)[0] /
                     a_bins)
            sigma3.append(H_nsc)
        if 'PartType4' in h.keys():
            H_stars = (
                np.histogram(r_stars, bins=r_bins, weights=mass_stars)[0]
                / a_bins)
            sigma4.append(H_stars)
        else:
            sigma4.append(np.full(len(r_bins)-1, 1e-99))

    sigma_data['Radius'] = (r_bins[:-1] + r_bins[1:]) / 2
    if RmaxHalo:
        sigma_data['RadiusDM'] = (r_bins_dm[:-1] + r_bins_dm[1:]) / 2

    if (N - n) > 0:
        sigma_data['SigmaGas'] = np.mean(np.array(sigma0), axis=0)
        sigma_data['SigmaGas10th'] = np.percentile(np.array(sigma0), q=10,
                                                   axis=0)
        sigma_data['SigmaGas90th'] = np.percentile(np.array(sigma0), q=90,
                                                   axis=0)
        sigma_data['SigmaDM'] = np.mean(np.array(sigma1), axis=0)
        sigma_data['SigmaDM10th'] = np.percentile(np.array(sigma1), q=10,
                                                  axis=0)
        sigma_data['SigmaDM90th'] = np.percentile(np.array(sigma1), q=90,
                                                  axis=0)
        if 'PartType2' in h.keys():
            sigma_data['SigmaDisc'] = np.mean(np.array(sigma2), axis=0)
            sigma_data['SigmaDisc10th'] = np.percentile(np.array(sigma2), q=10,
                                                        axis=0)
            sigma_data['SigmaDisc90th'] = np.percentile(np.array(sigma2), q=90,
                                                        axis=0)
        if 'PartType3' in h.keys():
            sigma_data['SigmaNSC'] = np.mean(np.array(sigma3), axis=0)
            sigma_data['SigmaNSC10th'] = np.percentile(np.array(sigma3), q=10,
                                                       axis=0)
            sigma_data['SigmaNSC90th'] = np.percentile(np.array(sigma3), q=90,
                                                       axis=0)
        if 'PartType4' in h.keys():
            sigma_data['SigmaStars'] = np.mean(np.array(sigma4), axis=0)
            sigma_data['SigmaStars10th'] = np.percentile(np.array(sigma4),
                                                         q=10, axis=0)
            sigma_data['SigmaStars90th'] = np.percentile(np.array(sigma4),
                                                         q=90, axis=0)
    else:
        sigma_data['SigmaGas'] = np.mean(np.array(sigma0), axis=0)
        sigma_data['SigmaDM'] = np.mean(np.array(sigma1), axis=0)
        if 'PartType2' in h.keys():
            sigma_data['SigmaDisc'] = np.mean(np.array(sigma2), axis=0)
        if 'PartType3' in h.keys():
            sigma_data['SigmaNSC'] = np.mean(np.array(sigma3), axis=0)
        if 'PartType4' in h.keys():
            sigma_data['SigmaStars'] = np.mean(np.array(sigma4), axis=0)

    return sigma_data
