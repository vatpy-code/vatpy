# Description: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-10-28


# -------------- Required packages
import numpy as np

from ..constants import const
from ..read import read_hdf5


# -------------- Declare function(s)
def get_surface_density_data(output_dir, n, N, r_bins_gas=None, r_bins_dm=None,
                             r_bins_disc=None, r_bins_nsc=None,
                             r_bins_stars=None, z_cut_halo=None,
                             z_cut_galaxy=None):
    '''TODO
    '''
    # Directory of surface densities:
    sigma_data = {}

    # Area bins:
    if r_bins_gas is not None:
        a_bins_gas = np.pi * (r_bins_gas[1:]**2 - r_bins_gas[:-1]**2)
    if r_bins_dm is not None:
        a_bins_dm = np.pi * (r_bins_dm[1:]**2 - r_bins_dm[:-1]**2)
    if r_bins_disc is not None:
        a_bins_disc = np.pi * (r_bins_disc[1:]**2 - r_bins_disc[:-1]**2)
    if r_bins_nsc is not None:
        a_bins_nsc = np.pi * (r_bins_nsc[1:]**2 - r_bins_nsc[:-1]**2)
    if r_bins_stars is not None:
        a_bins_stars = np.pi * (r_bins_stars[1:]**2 - r_bins_stars[:-1]**2)

    # Surface densities:
    sigma0 = []
    sigma1 = []
    sigma2 = []
    sigma3 = []
    sigma4 = []
    for i in range(n, N+1):
        snap = '000'[len(str(i)):] + str(i)
        print(f'Reading data of snapshot {snap}', end='\r')

        h, iu = read_hdf5(f'{output_dir}/snap_{snap}.hdf5')
        if 'PartType0' in h.keys() and r_bins_gas is not None:
            pos_gas = (h['PartType0']['Coordinates'] * iu['ulength'] /
                       const['pc'])
            mass_gas = h['PartType0']['Masses'] * iu['umass'] / const['Msol']
        if 'PartType1' in h.keys() and r_bins_dm is not None:
            pos_dm = (h['PartType1']['Coordinates'] * iu['ulength'] /
                      const['pc'])
            mass_dm = h['PartType1']['Masses'] * iu['umass'] / const['Msol']
        if 'PartType2' in h.keys() and r_bins_disc is not None:
            pos_disc = (h['PartType2']['Coordinates'] * iu['ulength'] /
                        const['pc'])
            mass_disc = h['PartType2']['Masses'] * iu['umass'] / const['Msol']
        if 'PartType3' in h.keys() and r_bins_disc is not None:
            pos_nsc = (h['PartType3']['Coordinates'] * iu['ulength'] /
                       const['pc'])
            mass_nsc = h['PartType3']['Masses'] * iu['umass'] / const['Msol']
        if 'PartType4' in h.keys() and r_bins_stars is not None:
            pos_stars = (h['PartType4']['Coordinates'] * iu['ulength'] /
                         const['pc'])
            mass_stars = h['PartType4']['Masses'] * iu['umass'] / const['Msol']
        pos_bh = h['PartType5']['Coordinates'][0] * iu['ulength'] / const['pc']

        # Centre the data on the BH:
        if 'PartType0' in h.keys() and r_bins_gas is not None:
            pos_gas -= pos_bh
        if 'PartType1' in h.keys() and r_bins_dm is not None:
            pos_dm -= pos_bh
        if 'PartType2' in h.keys() and r_bins_disc is not None:
            pos_disc -= pos_bh
        if 'PartType3' in h.keys() and r_bins_nsc is not None:
            pos_nsc -= pos_bh
        if 'PartType4' in h.keys() and r_bins_stars is not None:
            pos_stars -= pos_bh

        # Limit the data in Z:
        if z_cut_halo:
            if 'PartType1' in h.keys() and r_bins_dm is not None:
                mass_dm = mass_dm[np.abs(pos_dm[:, 2]) < z_cut_halo]
                pos_dm = pos_dm[np.abs(pos_dm[:, 2]) < z_cut_halo]

        if z_cut_galaxy:
            if 'PartType0' in h.keys() and r_bins_gas is not None:
                mass_gas = mass_gas[np.abs(pos_gas[:, 2]) < z_cut_galaxy]
                pos_gas = pos_gas[np.abs(pos_gas[:, 2]) < z_cut_galaxy]
            if 'PartType2' in h.keys() and r_bins_disc is not None:
                mass_disc = mass_disc[np.abs(pos_disc[:, 2]) < z_cut_galaxy]
                pos_disc = pos_disc[np.abs(pos_disc[:, 2]) < z_cut_galaxy]
            if 'PartType3' in h.keys() and r_bins_nsc is not None:
                mass_nsc = mass_nsc[np.abs(pos_nsc[:, 2]) < z_cut_galaxy]
                pos_nsc = pos_nsc[np.abs(pos_nsc[:, 2]) < z_cut_galaxy]
            if 'PartType4' in h.keys() and r_bins_stars is not None:
                mass_stars = mass_stars[np.abs(pos_stars[:, 2]) < z_cut_galaxy]
                pos_stars = pos_stars[np.abs(pos_stars[:, 2]) < z_cut_galaxy]

        # Radius in x and y:
        if 'PartType0' in h.keys() and r_bins_gas is not None:
            r_gas = np.linalg.norm(pos_gas[:, :2], axis=1)
        if 'PartType1' in h.keys() and r_bins_dm is not None:
            r_dm = np.linalg.norm(pos_dm[:, :2], axis=1)
        if 'PartType2' in h.keys() and r_bins_disc is not None:
            r_disc = np.linalg.norm(pos_disc[:, :2], axis=1)
        if 'PartType3' in h.keys() and r_bins_nsc is not None:
            r_nsc = np.linalg.norm(pos_nsc[:, :2], axis=1)
        if 'PartType4' in h.keys() and r_bins_stars is not None:
            r_stars = np.linalg.norm(pos_stars[:, :2], axis=1)

        # Surface densities:
        if 'PartType0' in h.keys() and r_bins_gas is not None:
            H_gas = (np.histogram(r_gas, bins=r_bins_gas, weights=mass_gas)[0]
                     / a_bins_gas)
            sigma0.append(H_gas)
        if 'PartType1' in h.keys() and r_bins_dm is not None:
            H_dm = (np.histogram(r_dm, bins=r_bins_dm, weights=mass_dm)[0] /
                    a_bins_dm)
            sigma1.append(H_dm)
        if 'PartType2' in h.keys() and r_bins_disc is not None:
            H_disc = (np.histogram(r_disc, bins=r_bins_disc,
                                   weights=mass_disc)[0] / a_bins_disc)
            sigma2.append(H_disc)
        if 'PartType3' in h.keys() and r_bins_nsc is not None:
            H_nsc = (np.histogram(r_nsc, bins=r_bins_nsc,
                                  weights=mass_nsc)[0] / a_bins_nsc)
            sigma3.append(H_nsc)
        if 'PartType4' in h.keys() and r_bins_stars is not None:
            H_stars = (
                np.histogram(r_stars, bins=r_bins_stars,
                             weights=mass_stars)[0] / a_bins_stars)
            sigma4.append(H_stars)
        else:
            sigma4.append(np.full(len(r_bins_stars)-1, 1e-99))

    if 'PartType0' in h.keys() and r_bins_gas is not None:
        sigma_data['RadiusGas'] = (r_bins_gas[:-1] + r_bins_gas[1:]) / 2
    if 'PartType1' in h.keys() and r_bins_dm is not None:
        sigma_data['RadiusDM'] = (r_bins_dm[:-1] + r_bins_dm[1:]) / 2
    if 'PartType2' in h.keys() and r_bins_disc is not None:
        sigma_data['RadiusDisc'] = (r_bins_disc[:-1] + r_bins_disc[1:]) / 2
    if 'PartType3' in h.keys() and r_bins_nsc is not None:
        sigma_data['RadiusNSC'] = (r_bins_nsc[:-1] + r_bins_nsc[1:]) / 2
    if 'PartType4' in h.keys() and r_bins_stars is not None:
        sigma_data['RadiusStars'] = (r_bins_stars[:-1] + r_bins_stars[1:]) / 2

    if (N - n) > 0:
        if 'PartType0' in h.keys() and r_bins_gas is not None:
            sigma_data['SigmaGas'] = np.mean(np.array(sigma0), axis=0)
            sigma_data['SigmaGas10th'] = np.percentile(np.array(sigma0), q=10,
                                                       axis=0)
            sigma_data['SigmaGas90th'] = np.percentile(np.array(sigma0), q=90,
                                                       axis=0)
        if 'PartType1' in h.keys() and r_bins_dm is not None:
            sigma_data['SigmaDM'] = np.mean(np.array(sigma1), axis=0)
            sigma_data['SigmaDM10th'] = np.percentile(np.array(sigma1), q=10,
                                                      axis=0)
            sigma_data['SigmaDM90th'] = np.percentile(np.array(sigma1), q=90,
                                                      axis=0)
        if 'PartType2' in h.keys() and r_bins_disc is not None:
            sigma_data['SigmaDisc'] = np.mean(np.array(sigma2), axis=0)
            sigma_data['SigmaDisc10th'] = np.percentile(np.array(sigma2), q=10,
                                                        axis=0)
            sigma_data['SigmaDisc90th'] = np.percentile(np.array(sigma2), q=90,
                                                        axis=0)
        if 'PartType3' in h.keys() and r_bins_nsc is not None:
            sigma_data['SigmaNSC'] = np.mean(np.array(sigma3), axis=0)
            sigma_data['SigmaNSC10th'] = np.percentile(np.array(sigma3), q=10,
                                                       axis=0)
            sigma_data['SigmaNSC90th'] = np.percentile(np.array(sigma3), q=90,
                                                       axis=0)
        if 'PartType4' in h.keys() and r_bins_stars is not None:
            sigma_data['SigmaStars'] = np.mean(np.array(sigma4), axis=0)
            sigma_data['SigmaStars10th'] = np.percentile(np.array(sigma4),
                                                         q=10, axis=0)
            sigma_data['SigmaStars90th'] = np.percentile(np.array(sigma4),
                                                         q=90, axis=0)
    else:
        if 'PartType0' in h.keys() and r_bins_gas is not None:
            sigma_data['SigmaGas'] = np.mean(np.array(sigma0), axis=0)
        if 'PartType1' in h.keys() and r_bins_dm is not None:
            sigma_data['SigmaDM'] = np.mean(np.array(sigma1), axis=0)
        if 'PartType2' in h.keys() and r_bins_disc is not None:
            sigma_data['SigmaDisc'] = np.mean(np.array(sigma2), axis=0)
        if 'PartType3' in h.keys() and r_bins_nsc is not None:
            sigma_data['SigmaNSC'] = np.mean(np.array(sigma3), axis=0)
        if 'PartType4' in h.keys() and r_bins_stars is not None:
            sigma_data['SigmaStars'] = np.mean(np.array(sigma4), axis=0)

    return sigma_data
