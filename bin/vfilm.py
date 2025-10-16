# Description: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-09-15


# -------------- Required Packages
import sys


# -------------- Import Sweep post-processing function
from vatpy import FilmMaker


# -------------- Run script
if sys.argv[1] == 's':
    snap = int(sys.argv[2])
    film = sys.argv[3]
    fm = FilmMaker(film=film)
    fm.generate_single_snap(snap=snap)
else:
    snap_min = int(sys.argv[1])
    snap_max = int(sys.argv[2])
    snap_range = (snap_min, snap_max)
    film = sys.argv[3]

    fm = FilmMaker(film=film)
    fm.generate(snap_range=snap_range)

# -------------- End of file
