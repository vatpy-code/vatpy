'''
Description:
Authour(s): Jonathan Petersson
Last updated: 2024-11-12
'''


# -------------- Required packages
import numpy as np

from ..read import read_hdf5, read_dump
from ..constants import const


# -------------- Declare function(s)
def get_black_hole_data(output_dir, n, N, spin=False, hm=False, rcirc=False,
                        vcr=False, rad=False, nfreq=5,
                        convert_to_arrays=False):
    # Data dictionary:
    black_hole_data = {
        'Time': [],
        'MassBH': [],
        'MassDisk': [],
        'MassReserv': [],
        'MassSink': [],
        'TimeAcc': [],
        'MdotBH': [],
        'MdotEdd': [],
        'MdotSink': [],
        'AngMom': [],
        'CircRadius': []
    }

    # Simple flag for when we use a varibale circ. radius:
    if vcr is True:
        spin, hm, rcirc = True, True, True

    # Loop over snapshots:
    for i in range(n, N+1):
        # Snapshot selection:
        snap = '000'[len(str(i)):] + str(i)
        if i < N:
            print(f'  * Reading data of snapshot {snap}', end='\r')
        else:
            print(f'  * Reading data of snapshot {snap}', end='\n')

        # Read HDF5 & binary dump file:
        h, iu = read_hdf5(f'{output_dir}/snap_{snap}.hdf5')
        dump = read_dump(f'{output_dir}/sink_snap_{snap}', feedback=False,
                         spin=spin, bh=True, hm=hm, rcirc=rcirc, rad=rad,
                         nfreq=nfreq)[2]

        # Append data to the dictionary:
        black_hole_data['Time'].append(h['Header'].attrs['Time']
                                       * iu['utime'] / const['Myr'])
        black_hole_data['MassBH'].append(dump['BlackHoleMass'][0]
                                         * iu['umass'] / const['Msol'])
        black_hole_data['MassDisk'].append(dump['BlackHoleDiskMass'][0]
                                           * iu['umass'] / const['Msol'])
        black_hole_data['MassReserv'].append(dump['BlackHoleReservoir'][0]
                                             * iu['umass'] / const['Msol'])
        black_hole_data['MassSink'].append(dump['Mass'][0] * iu['umass']
                                           / const['Msol'])
        if vcr == True:
            black_hole_data['AngMom'].append(dump['AngularMomentum'][0]
                                             * iu['uangmom'])
            black_hole_data['CircRadius'].append(dump['BlackHoleCircRadius'][0]
                                                 * iu['ulength'] / const['pc'])

    # Accretion rates:
    dt = (np.array(black_hole_data['Time'])[1:]
          - np.array(black_hole_data['Time'])[:-1])
    mt = ((np.array(black_hole_data['Time'])[1:]
           + np.array(black_hole_data['Time'])[:-1]) / 2)
    black_hole_data['TimeAcc'] = list(mt)
    mdot_sink = ((np.array(black_hole_data['MassSink'])[1:]
                  - np.array(black_hole_data['MassSink'])[:-1]) / (dt * 1e6))
    black_hole_data['MdotSink'] = list(mdot_sink)
    mdot_bh = ((np.array(black_hole_data['MassBH'])[1:]
                - np.array(black_hole_data['MassBH'])[:-1]) / (dt * 1e6))
    black_hole_data['MdotBH'] = list(mdot_bh)

    # Eddington accretion limit:
    mean_bh_mass = ((np.array(black_hole_data['MassBH'])[:-1]
                     + np.array(black_hole_data['MassBH'])[1:]) / 2)
    limit_Edd = (4 * np.pi * const['G'] * mean_bh_mass * const['Msol'] *
                 const['mp']) / (0.1 * const['thom'] * const['c'])
    black_hole_data['MdotEdd'] = list(limit_Edd / (const['Msol'] / const['yr']))

    # Convert the data into arrays instead of lists:
    if convert_to_arrays is True:
        for key in black_hole_data.keys():
            black_hole_data[key] = np.array(black_hole_data[key])

    return black_hole_data

# -------------- End of file
