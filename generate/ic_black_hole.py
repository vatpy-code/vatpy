'''
Description:
Authour(s): Jonathan Petersson
Last updated: 2025-05-13
'''

# -------------- Required Packages
import sys
import yaml

# -------------- Import write_dump
from vatpy import write_dump, read_dump


# -------------- Run script
print('\nWelcome to this IC generator script (powered by Vatpy):')
print('  * Generating a binary file for black hole sink particle(s)')

# Parameter file & group:
paramfile = sys.argv[1]
paramgroup = sys.argv[2]

print(f'  * Parameters: {paramgroup} in {paramfile}')

# Read parameters:
with open(paramfile, 'r') as file:
    f = yaml.safe_load(file)
    param = f[paramgroup]

    feedback = param['feedback']
    spin = param['spin']
    bh = param['bh']
    hm = param['hm']
    rcirc = param['rcirc']
    rad = param['rad']
    nfreq = param['nfreq']

filename = param['filename']
print(f'  * Writing data to: {filename}')
write_dump(filename=filename, ic=param, feedback=feedback, spin=spin, bh=bh,
           hm=hm, rcirc=rcirc, rad=rad, nfreq=nfreq)

print('  * Binary file generated!')
print('  * Reading newly created binary file to check that all data have' +
      ' been generated correctly')
dump = read_dump(param['filename'], feedback=feedback, spin=spin, bh=bh, hm=hm,
                 rcirc=rcirc, rad=rad, nfreq=nfreq)
print('  * Binary file content:')
for i in dump[2].keys():
    print(f'      {i}: {dump[2][i]}')

print('  * Done!')

# -------------- End of file
