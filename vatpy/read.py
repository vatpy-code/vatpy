'''
Description:
Authour(s):
Last updated: 2023-09-27
'''

# -------------- Required packages
import numpy as np
import h5py
import pycstruct


# -------------- Declare function(s)
def read_hdf5(file):
    '''
    Description: Read a HDF5 file output from Arepo and return the data +
                 internal units (in cgs).
    '''
    # Read hdf5-file:
    h = h5py.File(file, 'r')

    # Internal units:
    if 'Parameters' in h.keys():
        if 'UnitLength_in_cm' in h['Parameters'].attrs:
            group = 'Parameters'
        else:
            group = 'Header'
    else:
        group = 'Header'

    ulength = h[group].attrs['UnitLength_in_cm']
    umass = h[group].attrs['UnitMass_in_g']
    uvel = h[group].attrs['UnitVelocity_in_cm_per_s']
    utime = ulength/uvel
    udens = umass/(ulength**3)
    uaccel = uvel/utime
    uinterg = uvel**2
    uangmom = ulength * uvel * umass
    iu = {
        'ulength': ulength,
        'umass': umass,
        'uvel': uvel,
        'utime': utime,
        'udens': udens,
        'uaccel': uaccel,
        'uinterg': uinterg,
        'uangmom': uangmom
    }

    return h, iu

def read_dump(file, feedback=False, spin=False, bh=False, hm=False,
              rcirc=False, rad=False, nfreq=1):
    '''
    Description: Read a sink particle (binary) file output from Arepo. Please
                 note that there is high risk of mismatch between data fields
                 unless you know what the expected data structure is of your
                 simulation output.
    '''
    f = open(file, 'rb')

    time = np.fromfile(f, np.float64, 1)
    NSinksAllTasks = np.fromfile(f, np.uint32, 1)
    sinks = {}

    fields = ['Pos', 'Vel', 'Accel', 'Mass', 'FormationMass', 'FormationTime',
              'ID', 'HomeTask', 'Index', 'FormationOrder']

    if feedback is True:
        fields += ['N_sne', 'StellarMass', 'explosion_time',
                   'MassStillToConvert', 'AccretionTime']

    if spin is True:
        fields += ['AngularMomentum']

    if bh is True:
        fields += ['BlackHole']
        if hm is True:
            fields += ['BlackHoleHotMode']
        fields += ['BlackHoleAccRadius', 'BlackHoleMass', 'BlackHoleDiskMass',
                   'BlackHoleReservoir', 'BlackHoleSinkAccRate',
                   'CellsTotalMassBuffer']
        if rcirc is True:
            fields += ['BlackHoleCircRadius']
        if rad is True:
            fields += ['BlackHoleEddFrac', 'PhotoIonRate']

    for i in fields:
        sinks[i] = []

    for i in range(NSinksAllTasks[0]):
        struct = pycstruct.StructDef(alignment=8)
        struct.add('float64', 'Pos', shape=3)
        struct.add('float64', 'Vel', shape=3)
        struct.add('float64', 'Accel', shape=3)
        struct.add('float64', 'Mass')
        struct.add('float64', 'FormationMass')
        struct.add('float64', 'FormationTime')
        struct.add('uint64', 'ID')
        struct.add('uint32', 'HomeTask')
        struct.add('uint32', 'Index')
        struct.add('uint32', 'FormationOrder')
        if feedback is True:
            struct.add('uint32', 'N_sne')
            struct.add('float64', 'StellarMass')
            struct.add('float64', 'explosion_time', shape=800)
            struct.add('float64', 'MassStillToConvert', shape=50)
            struct.add('float64', 'AccretionTime', shape=50)
        if spin is True:
            struct.add('float64', 'AngularMomentum', shape=3)
        if bh is True:
            struct.add('uint32', 'BlackHole')
            if hm is True:
                struct.add('uint32', 'BlackHoleHotMode')
            struct.add('float64', 'BlackHoleAccRadius')
            struct.add('float64', 'BlackHoleMass')
            struct.add('float64', 'BlackHoleDiskMass')
            struct.add('float64', 'BlackHoleReservoir')
            struct.add('float64', 'BlackHoleSinkAccRate')
            struct.add('float64', 'CellsTotalMassBuffer')
            if rcirc is True:
                struct.add('float64', 'BlackHoleCircRadius')
            if rad is True:
                struct.add('float64', 'BlackHoleEddFrac')
                struct.add('float64', 'PhotoIonRate', shape=nfreq)

        inbytes = f.read(struct.size())
        data = struct.deserialize(inbytes)
        for field in fields:
            sinks[field] += [data[field]]

    for field in fields:
        sinks[field] = np.array(sinks[field])

    f.close()

    return time, NSinksAllTasks, sinks

# -------------- End of file
