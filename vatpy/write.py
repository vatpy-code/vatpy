'''
Description:
Authour(s): Jonathan Petersson
Last updated: 2024-11-12
'''

# -------------- Required packages
import numpy as np
import pycstruct


# -------------- Declare function(s)
def write_dump(filename, ic, feedback=False, spin=False, bh=False, hm=False,
               rcirc=False, sgs=False, rad=False, nfreq=1):
    '''
    Description:
    '''
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
            if sgs is True:
                struct.add('float64', 'BlackHoleDiskSpin', shape=3)
        if rad is True:
            struct.add('float64', 'BlackHoleEddFrac')
            struct.add('float64', 'PhotoIonRate', shape=nfreq)

    data = {}
    data['Pos'] = ic['sink_pos']
    data['Vel'] = ic['sink_vel']
    data['Accel'] = [0, 0, 0]
    data['Mass'] = ic['sink_mass']
    data['FormationMass'] = ic['sink_mass']
    data['FormationTime'] = 0
    data['ID'] = ic['sink_ID']
    data['HomeTask'] = 0
    data['Index'] = 0
    data['FormationOrder'] = 0
    if feedback is True:
        data['N_sne'] = 0
        data['StellarMass'] = 0
        data['explosion_time'] = np.zeros(800)
        data['MassStillToConvert'] = np.zeros(50)
        data['AccretionTime'] = np.zeros(50)
    if spin is True:
        data['AngularMomentum'] = ic['sink_angmom']
    if bh is True:
        data['BlackHole'] = 1
        if hm is True:
            data['BlackHoleHotMode'] = 0
        data['BlackHoleAccRadius'] = ic['bh_acc_radius']
        data['BlackHoleMass'] = ic['bh_mass']
        data['BlackHoleDiskMass'] = ic['bh_disk_mass']
        data['BlackHoleReservoir'] = ic['bh_reservoir']
        data['BlackHoleSinkAccRate'] = 0
        data['CellsTotalMassBuffer'] = 0
        if rcirc is True:
            data['BlackHoleCircRadius'] = 0
            if sgs is True:
                data['BlackHoleDiskSpin'] = [0, 0, 0]
        if rad is True:
            data['BlackHoleEddFrac'] = 0
            data['PhotoIonRate'] = np.zeros(nfreq)
    buffer = struct.serialize(data)

    f = open(filename, 'wb')
    np.array([0], dtype=np.double).tofile(f)
    np.array([1], dtype=np.intc).tofile(f)
    f.write(buffer)
    f.close()

    return None

# -------------- End of file
