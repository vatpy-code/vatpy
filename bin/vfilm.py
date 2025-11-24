# Description: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-09-15


# -------------- Required Packages
import sys


# -------------- Import Sweep post-processing function
from vatpy import FilmMaker


# -------------- Run script
if sys.argv[1] == 'noctua' or 'mosaic' or 'zoom' or 'ionflux' or 'ionrate':
    snap_min = int(sys.argv[2])
    snap_max = int(sys.argv[3])
    snap_range = (snap_min, snap_max)
    film = sys.argv[1]
    fm = FilmMaker(film=film)
    fm.generate(snap_range=snap_range)
elif sys.argv[1] == 'deepdive':
    snap = int(sys.argv[2])
    film = sys.argv[1]
    fm = FilmMaker(film=film)
    fm.generate_single_snap(snap=snap)
else:
    print('  * No matching film...')

# -------------- End of file
