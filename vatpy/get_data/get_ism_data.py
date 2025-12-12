# Description: TODO
# Last updated: 2025-09-08


# -------------- Required packages
import numpy as np

from ..constants import const
from ..get_gas_property import number_density, temperature, thermalpressure


# -------------- Declare function(s)
def get_ism_data(h, iu, bin_edges_dens, bin_edges_temp, norm=False, Rcut=None):
    '''TODO
    '''
    # Directory of ISM phase data:
    ism_data = {}

    # Gas cell properties:
    boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / const['kpc']
    pos = h['PartType0']['Coordinates'] * iu['ulength'] / const['kpc']
    mass = h['PartType0']['Masses'] * iu['umass'] / const['Msol']
    dens = h['PartType0']['Density'] * iu['udens']
    vol = mass * const['Msol'] / dens

    # Get number densities:
    num = number_density(h, iu)

    # Number density of H:
    num_H = num['HI'] + num['HII'] + 2 * num['H2']

    # Mass of H and H2:
    mass_H = num_H * vol * const['mp'] / const['Msol']
    mass_H2 = num['H2'] * vol * 2 * const['mp'] / const['Msol']

    # Gas temperature and abundances:
    temp = temperature(h, iu)
    abun = h['PartType0']['ChemicalAbundances'][:]
    xH2, xHII, xCO = abun[:, 0], abun[:, 1], abun[:, 2]
    xHI = 1 - 2*xH2 - xHII

    if Rcut:
        maskR = np.linalg.norm(pos - boxsize/2, axis=1) < Rcut
    else:
        maskR = np.full(len(pos), True)

    # ISM phases:
    ism_mask_list = {}

    ism_mask_list['CMM'] = maskR * (temp < 6e3) * (xH2 > 0.25)
    ism_mask_list['CNM'] = maskR * (temp < 5e2) * (xHI > 0.5)
    ism_mask_list['UNM'] = maskR * (temp > 5e2) * (temp < 6e3) * (xHI > 0.5)
    ism_mask_list['UIM'] = maskR * (temp < 6e3) * (xHII > 0.5)
    ism_mask_list['WNM'] = maskR * (temp > 6e3) * (temp < 3.5e4) * (xHI > 0.5)
    ism_mask_list['WPIM'] = (maskR * (temp > 6e3) * (temp < 1.5e4) *
                             (xHII > 0.5))
    ism_mask_list['WCIM'] = (maskR * (temp > 1.5e4) * (temp < 3.5e4) *
                             (xHII > 0.5))
    ism_mask_list['WHIM'] = maskR * (temp > 3.5e4) * (temp < 5e5)
    ism_mask_list['HIM'] = maskR * (temp > 5e5)

    ism_mask_list['CNM+CMM'] = (ism_mask_list['CNM'] + ism_mask_list['CMM'])
    ism_mask_list['WCIM+WPIM'] = (ism_mask_list['WCIM'] +
                                  ism_mask_list['WPIM'])
    ism_mask_list['HIM+WHIM'] = (ism_mask_list['HIM'] + ism_mask_list['WHIM'])

    for p in ism_mask_list.keys():
        mask = ism_mask_list[p]
        if len(num_H[mask]) > 0:
            phase_dens = num_H[mask]
            phase_temp = temp[mask]
            phase_weight = mass_H[mask]
            hist_dens = np.histogram(np.log10(phase_dens), bins=bin_edges_dens,
                                     weights=phase_weight)[0]
            hist_temp = np.histogram(np.log10(phase_temp), bins=bin_edges_temp,
                                     weights=phase_weight)[0]
            if norm:
                bin_widths_dens = bin_edges_dens[1:] - bin_edges_dens[:-1]
                bin_widths_temp = bin_edges_temp[1:] - bin_edges_temp[:-1]
                hist_dens = hist_dens / (np.sum(hist_dens) * bin_widths_dens)
                hist_temp = hist_temp / (np.sum(hist_temp) * bin_widths_temp)

            ism_data['hist_' + p + '_dens'] = hist_dens
            ism_data['hist_' + p + '_temp'] = hist_temp
        else:
            ism_data[p] = None

    ism_data['bin_edges_dens'] = bin_edges_dens
    ism_data['bin_edges_temp'] = bin_edges_temp

    ism_data['bin_mids_dens'] = (bin_edges_dens[:-1] + bin_edges_dens[1:]) / 2
    ism_data['bin_mids_temp'] = (bin_edges_temp[:-1] + bin_edges_temp[1:]) / 2

    return ism_data


def get_phase_diagram_data(h, iu, bins_dens, bins_temp, bins_pres, Rcut=None,
                           n_range=None, T_range=None):
    '''TODO
    '''
    # Directory of ISM phase data:
    diagram_data = {}

    # Gas cell properties:
    boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / const['kpc']
    pos = h['PartType0']['Coordinates'] * iu['ulength'] / const['kpc']
    mass = h['PartType0']['Masses'] * iu['umass'] / const['Msol']
    dens = h['PartType0']['Density'] * iu['udens']
    vol = mass * const['Msol'] / dens

    # Get number densities:
    num = number_density(h, iu)

    # Number density of H:
    num_H = num['HI'] + num['HII'] + 2 * num['H2']

    # Mass of H and H2:
    mass_H = num_H * vol * const['mp'] / const['Msol']
    mass_H2 = num['H2'] * vol * 2 * const['mp'] / const['Msol']

    # Gas temperature and abundances:
    temp = temperature(h, iu)

    # Gas pressure:
    thermpres = thermalpressure(h, iu)
    thermpres = thermpres / const['kb']

    # Abundances:
    abun = h['PartType0']['ChemicalAbundances'][:]
    xH2, xHII, xCO = abun[:, 0], abun[:, 1], abun[:, 2]
    xHI = 1 - 2*xH2 - xHII

    if Rcut:
        maskR = np.linalg.norm(pos - boxsize/2, axis=1) < Rcut
    else:
        maskR = np.full(len(pos), True)

    # Temperature vs density:
    if (n_range is not None) and (T_range is not None):
        H, xedges, yedges = np.histogram2d(np.log10(num_H[maskR]),
                                           np.log10(temp[maskR]),
                                           bins=[np.linspace(n_range[0],
                                                             n_range[1], bins_dens),
                                                np.linspace(T_range[0],
                                                            T_range[1], bins_temp)],
                                           weights=mass_H[maskR])
    else:
        H, xedges, yedges = np.histogram2d(np.log10(num_H[maskR]),
                                           np.log10(temp[maskR]),
                                           bins=(bins_dens, bins_temp),
                                           weights=mass_H[maskR])

    with np.errstate(divide='ignore'):
        diagram_data['hist2d_T_vs_rho'] = np.log10(H.T)
        diagram_data['hist2d_T_vs_rho_bin_xedges'] = xedges
        diagram_data['hist2d_T_vs_rho_bin_yedges'] = yedges

    # Pressure vs density:
    H, xedges, yedges = np.histogram2d(np.log10(num_H[maskR]),
                                       np.log10(thermpres[maskR]),
                                       bins=(bins_dens, bins_pres),
                                       weights=mass_H[maskR])

    with np.errstate(divide='ignore'):
        diagram_data['hist2d_Ptherm_vs_rho'] = np.log10(H.T)
        diagram_data['hist2d_Ptherm_vs_rho_bin_xedges'] = xedges
        diagram_data['hist2d_Ptherm_vs_rho_bin_yedges'] = yedges

    return diagram_data

# -------------- End of file
