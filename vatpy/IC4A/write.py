'''
Description: TODO
Authour(s): Jonathan Petersson
Last updated: 2025-05-16
'''


# -------------- Required packages
import numpy as np
import h5py


# -------------- Declare function(s)
def write_ic_file(filename, savepath, boxsize, parttypes, masstable=None):
    # Number of particle types:
    NumPartTypes = 6

    with h5py.File(f'{savepath}/{filename}.hdf5', 'w') as f:
        # Write particle types and their respective fields to HDF5 file:
        for pt in parttypes.keys():
            group = f.create_group(pt)
            for field in parttypes[pt]:
                group[field] = parttypes[pt][field]

        # Set particle counts:
        CountPart = np.zeros(NumPartTypes, dtype='int64')
        for pt in parttypes.keys():
            n = int(pt[-1])
            CountPart[n] = len(parttypes[pt]['ParticleIDs'])

        # Header:
        h = f.create_group('Header')
        h.attrs['BoxSize'] = boxsize
        h.attrs['NumFilesPerSnapshot'] = 1
        h.attrs['NumPart_ThisFile'] = CountPart
        h.attrs['NumPart_Total'] = CountPart
        h.attrs['NumPart_Total_HighWord'] = np.zeros(NumPartTypes)

        for k in ['Time', 'Redshift', 'Omega0', 'OmegaLambda', 'HubbleParam']:
            h.attrs[k] = 0.0
        for k in ['Sfr', 'Cooling', 'StellarAge', 'Metals', 'Feedback']:
            h.attrs[f'Flag_{k}'] = 0
        h.attrs['Flag_DoublePrecision'] = 1

        if masstable is not None:
            h.attrs['MassTable'] = masstable
        else:
            h.attrs['MassTable'] = np.zeros(NumPartTypes)

    print('  * IC generated!')
    print(f'  * File name: {filename}.hdf5')
    print(f'  * Save path: {savepath}')

    return None


# -------------- End of file
