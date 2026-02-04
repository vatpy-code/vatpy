# Description: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-09-15


# -------------- Required Packages
import sys


# -------------- Import Sweep post-processing function
from vatpy import FilmMaker


# -------------- Run script
avail_films = ['noctua', 'mosaic', 'zoom', 'ionflux', 'ionrate', 'ionfluxhuv']
avail_films_single_snapshot = ['deepdive']

if sys.argv[1] in avail_films:
    film = sys.argv[1]
    snap_min = int(sys.argv[2])
    snap_max = int(sys.argv[3])
    snap_range = (snap_min, snap_max)
    fm = FilmMaker(film=film)
    fm.generate(snap_range=snap_range)
elif sys.argv[1] in avail_films_single_snapshot:
    film = sys.argv[1]
    snap = int(sys.argv[2])
    fm = FilmMaker(film=film)
    fm.generate_single_snap(snap=snap)
else:
    print('  * No matching film...')

# -------------- End of file
