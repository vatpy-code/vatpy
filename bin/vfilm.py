# Description: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-09-15


# -------------- Required Packages
import sys


# -------------- Import Sweep post-processing function
from vatpy import FilmMaker


# -------------- Run script
snap_min = sys.argv[1]
snap_max = sys.argv[2]
snap_range = (snap_min, snap_max)
film = sys.argv[3]

fm = FilmMaker(film=film, snap_range=snap_range)
fm.generate()

# -------------- End of file
