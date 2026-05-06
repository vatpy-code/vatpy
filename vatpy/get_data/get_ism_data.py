# Description: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-09-08


# -------------- Required packages
import numpy as np

from ..read import read_hdf5
from ..write import write_datadict_to_file
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


def get_phase_diagram_data(output_dir, N, bins_dens, bins_temp, bins_pres=None,
                           bins_abun=None, Rcut=None, colorcode='H'):
    '''TODO
    '''
    # Snapshot selection:
    snap = '000'[len(str(N)):] + str(N)
    print(f'Reading data of snapshot {snap}', end='\r')

    # Read HDF5 & binary dump file:
    h, iu = read_hdf5(f'{output_dir}/snap_{snap}.hdf5')

    # Data directory:
    data = {}

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

    # Gas temperature:
    temp = temperature(h, iu)

    # Gas pressure:
    thermpres = thermalpressure(h, iu)
    thermpres = thermpres / const['kb']

    # Abundances:
    abun = h['PartType0']['ChemicalAbundances'][:]
    xH2, xHII, xCO = abun[:, 0], abun[:, 1], abun[:, 2]
    xHI = 1 - 2*xH2 - xHII

    if Rcut:
        mask_radius = np.linalg.norm(pos - boxsize/2, axis=1) < Rcut
    else:
        mask_radius = np.full(len(pos), True)

    # Weights:
    if colorcode == 'H':
        weights = mass_H
    elif colorcode == 'H2':
        weights = mass_H2

    # Temperature vs density:
    H, xedges, yedges = np.histogram2d(np.log10(num_H[mask_radius]),
                                       np.log10(temp[mask_radius]),
                                       bins=(bins_dens, bins_temp),
                                       weights=weights[mask_radius])
    with np.errstate(divide='ignore'):
        data['hist2d_T_vs_rho'] = np.log10(H.T)
        data['hist2d_T_vs_rho_bin_xedges'] = xedges
        data['hist2d_T_vs_rho_bin_yedges'] = yedges
        data['extent_T_vs_rho'] = [xedges[0], xedges[-1], yedges[0],
                                   yedges[-1]]

    # Pressure vs density:
    if bins_pres:
        H, xedges, yedges = np.histogram2d(np.log10(num_H[mask_radius]),
                                           np.log10(thermpres[mask_radius]),
                                           bins=(bins_dens, bins_pres),
                                           weights=weights[mask_radius])
        with np.errstate(divide='ignore'):
            data['hist2d_Ptherm_vs_rho'] = np.log10(H.T)
            data['hist2d_Ptherm_vs_rho_bin_xedges'] = xedges
            data['hist2d_Ptherm_vs_rho_bin_yedges'] = yedges
            data['extent_Ptherm_vs_rho'] = [xedges[0], xedges[-1], yedges[0],
                                            yedges[-1]]

    # Hydrogen abundance vs temperature:
    if bins_abun:
        H, xedges, yedges = np.histogram2d(np.log10(temp[mask_radius]),
                                           xHI[mask_radius],
                                           bins=(bins_temp, bins_abun),
                                           weights=weights[mask_radius])
        with np.errstate(divide='ignore'):
            data['hist2d_xHI_vs_T'] = np.log10(H.T)
            data['hist2d_xHI_vs_T_bin_xedges'] = xedges
            data['hist2d_xHI_vs_T_bin_yedges'] = yedges
            data['extent_xHI_vs_T'] = [xedges[0], xedges[-1], yedges[0],
                                       yedges[-1]]

    return data


def get_chem_evolve_data(output_dir, n, N, Rcut=None, Zmax=None, write=False,
                         name=None):
    data = {
        'Time': [],
        'Temp': [],
        'MassFracH2': [],
        'MassFracHII': [],
        'MassFracCold': [],
        'MassFracHot': [],
        'MassFracWNM': [],
        'MassFracWIM': [],
    }

    if Rcut:
        [data['Temp'].append([]) for i in range(0, len(Rcut))]
        [data['MassFracH2'].append([]) for i in range(0, len(Rcut))]
        [data['MassFracHII'].append([]) for i in range(0, len(Rcut))]
        [data['MassFracCold'].append([]) for i in range(0, len(Rcut))]
        [data['MassFracHot'].append([]) for i in range(0, len(Rcut))]
        [data['MassFracWNM'].append([]) for i in range(0, len(Rcut))]
        [data['MassFracWIM'].append([]) for i in range(0, len(Rcut))]

    for i in range(n, N+1):
        # Snapshot selection:
        snap = '000'[len(str(i)):] + str(i)
        print(f'Reading data of snapshot {snap}', end='\r')

        # Read HDF5 & binary dump file:
        h, iu = read_hdf5(f'{output_dir}/snap_{snap}.hdf5')

        # Time:
        time = h['Header'].attrs['Time'] * iu['utime'] / const['Myr']
        data['Time'].append(time)

        # MBH:
        bh = h['PartType5']['Coordinates'][0] * iu['ulength'] / const['kpc']

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
        mass_HII = num['HII'] * vol * const['mp'] / const['Msol']
        mass_H2 = num['H2'] * vol * 2 * const['mp'] / const['Msol']

        # Gas temperature and abundances:
        temp = temperature(h, iu)
        abun = h['PartType0']['ChemicalAbundances'][:]
        xH2, xHII, xCO = abun[:, 0], abun[:, 1], abun[:, 2]
        xHI = 1 - 2*xH2 - xHII

        # Masks:
        maskCOLD = temp < 300
        maskHOT = temp > 5e5
        maskWNM = ((temp > 300) * (temp < 5e5) * (xHI > 0.5))
        maskWIM = ((temp > 300) * (temp < 5e5) * (xHII > 0.5))

        # Masks in R bins:
        if Rcut:
            for j in range(0, len(Rcut)):
                # Geometric mask (if R > Zmax, mask in a cylinder):
                if Zmax:
                    if Rcut[j] > Zmax:
                        maskR = ((np.linalg.norm(pos[:, :2] - bh[:2], axis=1) <
                                  Rcut[j]) *
                                 (np.abs(pos[:, 2] - bh[2]) < Zmax))
                    else:
                        maskR = np.linalg.norm(pos - bh, axis=1) < Rcut[j]
                else:
                    maskR = np.linalg.norm(pos - bh, axis=1) < Rcut[j]

                # Temp, H2 & HII:
                if len(pos[maskR]) > 0:
                    data['Temp'][j].append(
                        np.average(temp[maskR], weights=mass[maskR]))
                    data['MassFracH2'][j].append(
                        np.sum(mass_H2[maskR]) / np.sum(mass_H[maskR]))
                    data['MassFracHII'][j].append(
                        np.sum(mass_HII[maskR]) / np.sum(mass_H[maskR]))
                else:
                    data['Temp'][j].append(0)
                    data['MassFracH2'][j].append(0)
                    data['MassFracHII'][j].append(0)

                # Phases:
                # Cold:
                if len(pos[maskR * maskCOLD]) > 0:
                    data['MassFracCold'][j].append(
                        np.sum(mass_H[maskR * maskCOLD]) /
                        np.sum(mass_H[maskR]))
                else:
                    data['MassFracCold'][j].append(0)

                # Hot:
                if len(pos[maskR * maskHOT]) > 0:
                    data['MassFracHot'][j].append(
                        np.sum(mass_H[maskR * maskHOT]) /
                        np.sum(mass_H[maskR]))
                else:
                    data['MassFracHot'][j].append(0)

                # WNM:
                if len(pos[maskR * maskWNM]) > 0:
                    data['MassFracWNM'][j].append(
                        np.sum(mass_H[maskR * maskWNM]) /
                        np.sum(mass_H[maskR]))
                else:
                    data['MassFracWNM'][j].append(0)

                # WIM:
                if len(pos[maskR * maskWIM]) > 0:
                    data['MassFracWIM'][j].append(
                        np.sum(mass_H[maskR * maskWIM]) /
                        np.sum(mass_H[maskR]))
                else:
                    data['MassFracWIM'][j].append(0)

        else:
            maskR = np.full(len(pos), True)

            data['Temp'].append(np.average(temp[maskR], weights=mass[maskR]))
            data['MassFracH2'].append(np.sum(mass_H2[maskR]) /
                                      np.sum(mass_H[maskR]))
            data['MassFracHII'].append(np.sum(mass_HII[maskR]) /
                                       np.sum(mass_H[maskR]))
            data['MassFracCold'].append(np.sum(mass_H[maskR * maskCOLD]) /
                                        np.sum(mass_H[maskR]))
            data['MassFracHot'].append(np.sum(mass_H[maskR * maskHOT]) /
                                       np.sum(mass_H[maskR]))
            data['MassFracWNM'].append(np.sum(mass_H[maskR * maskWNM]) /
                                       np.sum(mass_H[maskR]))
            data['MassFracWIM'].append(np.sum(mass_H[maskR * maskWIM]) /
                                       np.sum(mass_H[maskR]))

    if write is True:
        write_datadict_to_file(data, name=name)

    return data

# -------------- End of file
