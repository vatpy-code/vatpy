# Description: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-09-12


# -------------- Required packages
# N/A

# -------------- Declare function(s)
def mplfigsize(option, ratio=None):
    '''TODO
    '''
    # Figure size options:
    mm = 1 / 25.4  # [in]
    figsizes = {
        'A&A_1col': (88 * mm, 88 * mm),
        'A&A_2col': (120 * mm, 120 * mm)
    }

    # Select figure size:
    width, height = figsizes[option]

    if ratio:
        height /= ratio

    return (width, height)

# -------------- End of file
