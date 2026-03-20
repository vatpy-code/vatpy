'''
Description: TODO
Authour(s): Jonathan Petersson
Last updated: 2026-03-17
'''


# -------------- Required packages
import numpy as np

from ..read import read_hdf5
from ..constants import const


# -------------- Declare function(s)
def get_star_formation_data(output_dir, n, N, dt, H=None, rbins=None):
    data = {}

    time_list = []
    birthtime_list = []
    birthmass_list = []
    totalgasmass_list = []
    if rbins is not None:
        [birthtime_list.append([]) for i in rbins]
        [birthmass_list.append([]) for i in rbins]
        [totalgasmass_list.append([]) for i in rbins]
    halfmassradius_list = []

    counter = 0
    for i in range(n, N+1):
        # Snapshot selection:
        snap = '000'[len(str(i)):] + str(i)
        print(f'Reading data of snapshot {snap}', end='\r')

        # Read HDF5 & binary dump file:
        h, iu = read_hdf5(f'{output_dir}/snap_{snap}.hdf5')
        time = h['Header'].attrs['Time'] * iu['utime'] / const['Myr']
        boxsize = h['Header'].attrs['BoxSize'] * iu['ulength'] / const['kpc']

        # Counter:
        if i == n:
            # counter = time
            data['Time0'] = time

        # Half-mass radius:
        if H is not None:
            N2 = np.sum(h['Header'].attrs['NumPart_Total'][2])
            N4 = np.sum(h['Header'].attrs['NumPart_Total'][4])
            totalstellarmass = np.zeros(int(N2 + N4))
            totalstellarpos = np.zeros((int(N2 + N4), 3))

            if 'PartType2' in h.keys():
                totalstellarmass[:N2] = (h['PartType2']['Masses'] * iu['umass']
                                         / const['Msol'])
                totalstellarpos[:N2] = (h['PartType2']['Coordinates'] *
                                        iu['ulength'] / const['kpc'])

            if 'PartType4' in h.keys():
                totalstellarmass[N2:] = (h['PartType4']['Masses'] * iu['umass']
                                         / const['Msol'])
                totalstellarpos[N2:] = (h['PartType4']['Coordinates'] *
                                        iu['ulength'] / const['kpc'])

            # Half-mass radius calculation:
            if (N2 + N4) > 0:
                totalstellarradius = np.linalg.norm(totalstellarpos -
                                                    boxsize/2, axis=1)
                tsm_sorted = totalstellarmass[np.argsort(totalstellarradius)]
                csm = np.cumsum(tsm_sorted)
                tsr_sorted = totalstellarradius[np.argsort(totalstellarradius)]
                halfmassindex = np.argmin(np.abs(csm - np.sum(tsm_sorted) / 2))
                halfmassradius = tsr_sorted[halfmassindex]
            else:
                halfmassradius = 0

            # Append half mass radius:
            halfmassradius_list.append(halfmassradius)

        # Total gas mass:
        mass_gas = (h['PartType0']['Masses'] * iu['umass'] / const['Msol'])
        pos_gas = h['PartType0']['Coordinates'] * iu['ulength'] / const['kpc']

        # Within radial bins:
        if rbins is not None:
            for j in range(0, len(rbins)):
                mask_radius = (np.linalg.norm(pos_gas - boxsize/2, axis=1) <
                               rbins[j])
                if len(mass_gas[mask_radius] > 0):
                    totalgasmass = np.sum(mass_gas[mask_radius])
                else:
                    totalgasmass = 0
                totalgasmass_list[j].append(totalgasmass)

        # Within H times the half mass radius:
        if H is not None:
            mask_radius = (np.linalg.norm(pos_gas - boxsize/2, axis=1) <
                           (H * halfmassradius))
            if len(mass_gas[mask_radius] > 0):
                totalgasmass = np.sum(mass_gas[mask_radius])
            else:
                totalgasmass = 0
            totalgasmass_list.append(totalgasmass)

        # Star formation:
        if 'PartType4' in h.keys():
            pos_stars = (h['PartType4']['Coordinates'] * iu['ulength'] /
                         const['kpc'])
            birthtime = (h['PartType4']['StellarFormationTime'] * iu['utime'] /
                         const['Myr'])

            if 'InitialMass' in h['PartType4'].keys():
                birthmass = (h['PartType4']['InitialMass'] * iu['umass'] /
                             const['Msol'])
            else:
                birthmass = (h['PartType4']['Masses'] * iu['umass'] /
                             const['Msol'])

            mask_newstars = birthtime > counter

            if rbins is not None:
                for j in range(0, len(rbins)):
                    mask_radius = (np.linalg.norm(pos_stars - boxsize/2,
                                                  axis=1) < rbins[j])
                    if len(birthtime[mask_newstars * mask_radius] > 0):
                        [birthtime_list[j].append(k) for k in
                         birthtime[mask_newstars * mask_radius]]
                        [birthmass_list[j].append(k) for k in
                         birthmass[mask_newstars * mask_radius]]

            if H is not None:
                mask_radius = (np.linalg.norm(pos_stars - boxsize/2, axis=1) <
                               (H * halfmassradius))
                if len(birthtime[mask_newstars * mask_radius] > 0):
                    [birthtime_list.append(j) for j in
                     birthtime[mask_newstars * mask_radius]]
                    [birthmass_list.append(j) for j in
                     birthmass[mask_newstars * mask_radius]]

        time_list.append(time)
        counter = time

    # SFR calculation (with or without radial bins):
    bins_time = np.arange(0, counter+dt, dt)
    if rbins is not None:
        SFR_list = []
        CSF_list = []
        for i in range(0, len(rbins)):
            H_sfr, bin_edges = np.histogram(birthtime_list[i], bins=bins_time,
                                            weights=birthmass_list[i])
            bin_mids = (bin_edges[:-1] + bin_edges[1:]) / 2
            SFR_list.append(H_sfr / (dt * 1e6))
            CSF_list.append(np.cumsum(H_sfr[bin_mids > data['Time0']]))
        data['SFR'] = SFR_list
        data['CSF'] = CSF_list
    else:
        H_sfr, bin_edges = np.histogram(birthtime_list, bins=bins_time,
                                        weights=birthmass_list)
        bin_mids = (bin_edges[:-1] + bin_edges[1:]) / 2
        SFR = H_sfr / (dt * 1e6)
        data['SFR'] = SFR

    data['Time'] = bin_mids

    # sSFR calculation:
    if rbins is not None:
        sSFR_list = []
        totalgasmass_mean_list = []
        for i in range(0, len(rbins)):
            totalgasmass_mean_list.append(
                np.histogram(time_list, bins=bins_time,
                             weights=totalgasmass_list[i])[0]
                / np.histogram(time_list, bins=bins_time)[0])
            sSFR_list.append(SFR_list[i] / totalgasmass_mean_list[i])
        data['sSFR'] = sSFR_list
        data['GasMass'] = totalgasmass_mean_list
    else:
        totalgasmass_mean = (np.histogram(time_list, bins=bins_time,
                                          weights=totalgasmass_list)[0] /
                             np.histogram(time_list, bins=bins_time)[0])
        sSFR = SFR / totalgasmass_mean
        data['sSFR'] = sSFR
        data['GasMass'] = totalgasmass_mean

    if H is not None:
        halfmassradius_mean = (np.histogram(time_list, bins=bins_time,
                                            weights=halfmassradius_list)[0] /
                               np.histogram(time_list, bins=bins_time)[0])
        data['HalfMassRadius'] = halfmassradius_mean

    return data


def get_star_formation_data_simple(output_dir, N, dt):
    data = {}

    # Snapshot selection:
    snap = '000'[len(str(N)):] + str(N)
    print(f'Reading data of snapshot {snap}', end='\r')

    # Read HDF5 & binary dump file:
    h, iu = read_hdf5(f'{output_dir}/snap_{snap}.hdf5')
    time = h['Header'].attrs['Time'] * iu['utime'] / const['Myr']

    # Star particles:
    birthtime = (h['PartType4']['StellarFormationTime'] * iu['utime'] /
                 const['Myr'])
    if 'InitialMass' in h['PartType4'].keys():
        birthmass = (h['PartType4']['InitialMass'] * iu['umass'] /
                     const['Msol'])
    else:
        birthmass = h['PartType4']['Masses'] * iu['umass'] / const['Msol']

    # SFR calculation:
    bins_time = np.arange(0, time+dt, dt)
    H, bin_edges = np.histogram(birthtime, bins=bins_time, weights=birthmass)
    bin_mids = (bin_edges[:-1] + bin_edges[1:]) / 2
    SFR = H / (dt * 1e6)
    data['Time'] = bin_mids
    data['SFR'] = SFR

    return data

# -------------- End of file
