'''
Description: TODO
Authour(s): Jonathan Petersson
Last updated: 2025-09-01
'''

# -------------- Required Packages
import sys

# -------------- Import Sweep post-processing function
from vatpy import do_sweep_post_process


# -------------- Run script
source_dir = sys.argv[1]
snap_min = sys.argv[2]
snap_max = sys.argv[3]
snap_range = (snap_min, snap_max)

do_sweep_post_process(source_dir=source_dir, snap_range=snap_range)

# -------------- End of file
