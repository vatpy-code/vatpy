'''
Description:
Authour(s): Jonathan Petersson
Last updated: 2024-11-12
'''

# -------------- Required packages
import numpy as np
import h5py

from scipy.interpolate import NearestNDInterpolator

from .read import read_hdf5, read_dump

# -------------- Declare function(s)
def convert_pNbody_IC_to_Arepo(source, dest, addbackgroundgrid=False, N=int(1e5), backgrounddensity_in_cgs=1e-30, addBH = False, mass=1e5):
    '''
    Description: 
    '''
    # Read the source file:
    print('  * Copying source file to assigned destination')
    with h5py.File(source, 'r') as src:
        # Create the destination file:
        with h5py.File(dest, 'w') as dst:
            # Copy the data from source to destination:
            for key in src.keys():
                src.copy(key, dst)

    # Open the destination file and make edits:
    print('  * Opening copied file to make edits')
    with h5py.File(dest, 'r+') as f:
        # Delete unnecessary groups:
        print('  * Deleting unnecessary data groups')
        if 'Code' in f:
            del f['Code']
        if 'Cosmology' in f:
            del f['Cosmology']
        if 'HydroScheme' in f:
            del f['HydroScheme']
        if 'RuntimePars' in f:
            del f['RuntimePars']

        # Get the units and then delete the group:
        print('  * Retriving internal units')
        ulength = float(f['Units'].attrs['Unit length in cgs (U_L)'])
        umass = float(f['Units'].attrs['Unit mass in cgs (U_M)'])
        utime = float(f['Units'].attrs['Unit time in cgs (U_t)'])
        uvel = ulength / utime
        udens = umass / np.power(ulength, 3)
        del f['Units']
        
        # Get important header info:
        print('  * Obtaining important header information')
        time = float(f['Header'].attrs['Time'])
        boxsize = int(f['Header'].attrs['BoxSize'])
        masstable = f['Header'].attrs['MassTable']
        num_part_thisfile = f['Header'].attrs['NumPart_ThisFile']
        num_part_total = f['Header'].attrs['NumPart_Total']

        # Delete the header, and create a new one:
        print('  * Creating a new header')
        del f['Header']
        h = f.create_group('Header')

        h.attrs['Time'] = time
        h.attrs['BoxSize'] = boxsize
        h.attrs['NumPart_ThisFile'] = num_part_thisfile
        h.attrs['NumPart_Total'] = num_part_total

        h.attrs['MassTable'] = np.zeros(6)
        h.attrs['NumFilesPerSnapshot'] = 1
        h.attrs['NumPart_Total_HighWord'] = np.zeros(6)

        for i in ['Redshift', 'Omega0', 'OmegaLambda', 'HubbleParam']:
            h.attrs[i] = 0.0
        for i in ['Sfr', 'Cooling', 'StellarAge', 'Metals', 'Feedback']:
    	    h.attrs[f'Flag_{i}'] = 0
        h.attrs['Flag_DoublePrecision'] = 1

        # Added the units to the header:
        h.attrs['UnitLength_in_cm'] = ulength
        h.attrs['UnitMass_in_g'] = umass
        h.attrs['UnitVelocity_in_cm_per_s'] = uvel
        
        # Loop over the particle types:
        print('  * Looping over particle types to remove unnecessary data')
        for i in range(0, 6):
            if f'PartType{i}' in f:
                if 'InternalEnergy' in f[f'PartType{i}']:
                    del f[f'PartType{i}']['InternalEnergy']
                if 'SmoothingLength' in f[f'PartType{i}']:
                    del f[f'PartType{i}']['SmoothingLength']
                if 'StellarParticleType' in f[f'PartType{i}']:
                    del f[f'PartType{i}']['StellarParticleType']
                if 0 in f[f'PartType{i}']['ParticleIDs']:
                    partIDs = f[f'PartType{i}']['ParticleIDs']
                    partIDs[partIDs == 0] = np.sum(num_part_total)

        if ('PartType4' in f) and ('PartType2' not in f):
            print('  * Making particles assigned to PartType4, into PartType2 instead')
            pt2 = f.create_group('PartType2')
            for i in f['PartType4'].keys():
                pt2[i] = f['PartType4'][i][:]
            del f['PartType4']

            n4 = num_part_total[4]
            num_part_thisfile[2] = n4
            num_part_total[2] = n4
            num_part_thisfile[4] = 0
            num_part_total[4] = 0

            h.attrs['NumPart_ThisFile'] = num_part_thisfile
            h.attrs['NumPart_Total'] = num_part_total

        if addbackgroundgrid == True:
            print('  * Adding a background grid of gas cells')
            grid_pos = boxsize * np.random.rand(N, 3) - boxsize/2
            interp = NearestNDInterpolator(f['PartType1']['Coordinates'][:], 
                                           f['PartType1']['Velocities'][:])
            grid_vel = interp(grid_pos)
            avg_vol = np.power(boxsize, 3) / N
            grid_mass = np.full(N, (backgrounddensity_in_cgs / udens) * avg_vol)
            grid_ids = np.arange(0, N, 1) + (np.sum(num_part_total) + 1)
    
            pos = f['PartType0']['Coordinates'][:]
            pos_new = np.append(pos, grid_pos, axis=0)
            del f['PartType0']['Coordinates']
            f['PartType0'].create_dataset('Coordinates', data=pos_new)

            vel = f['PartType0']['Velocities'][:]
            vel_new = np.append(vel, grid_vel, axis=0)
            del f['PartType0']['Velocities']
            f['PartType0'].create_dataset('Velocities', data=vel_new)
    
            mass = f['PartType0']['Masses'][:]
            mass_new = np.append(mass, grid_mass, axis=0)
            del f['PartType0']['Masses']
            f['PartType0'].create_dataset('Masses', data=mass_new)
    
            ids = f['PartType0']['ParticleIDs'][:]
            ids_new = np.append(ids, grid_ids, axis=0)
            del f['PartType0']['ParticleIDs']
            f['PartType0'].create_dataset('ParticleIDs', data=ids_new)
    
            num_part_thisfile[0] += N
            num_part_total[0] += N
            
            h.attrs['NumPart_ThisFile'] = num_part_thisfile
            h.attrs['NumPart_Total'] = num_part_total        

        if addBH == True:
            print('  * Adding a central black hole,')
            pt5 = f.create_group('PartType5')

            id_BH = np.sum(num_part_total) + 1

            pt5.create_dataset('ParticleIDs', data=np.array([id_BH]))
            pt5.create_dataset('Coordinates', data=np.array([[0.0], [0.0], [0.0]]).T)
            pt5.create_dataset('Velocities', data=np.array([[0.0], [0.0], [0.0]]).T)
            pt5.create_dataset('Masses', data=np.array([mass]))

            num_part_thisfile[5] = 1
            num_part_total[5] = 1

            h.attrs['NumPart_ThisFile'] = num_part_thisfile
            h.attrs['NumPart_Total'] = num_part_total
            
    print('  * Done!')
    
    return None
    
# -------------- End of file
