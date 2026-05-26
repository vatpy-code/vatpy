# Description: TODO
# Authour(s): Jonathan Petersson
# Last updated: 2025-09-12


# -------------- Required packages
# N/A

# -------------- Declare function(s)
def mplfigsize(option='A&A_1col', ratio=None):
    '''TODO
    '''
    # Figure size options:
    mm = 1 / 25.4  # [in]
    pt = 1 / 64    # [in]
    figsizes = {
        'A&A_1col': (88 * mm, 88 * mm),
        'A&A_2col': (170 * mm, 170 * mm),
        'A&A_2col_sc': (120 * mm, 120 * mm),
        'ApJ_1col': (242.3 * pt, 242.3 * pt),
        'ApJ_2col': (2 * 242.3 * pt, 2 * 242.3 * pt),
    }

    # Select figure size:
    width, height = figsizes[option]

    if ratio:
        height /= ratio

    return (width, height)


def mplstyle(plt, option='AA'):
    '''TODO
    '''
    plt.style.use(f'~/vatpy/mpl/{option}.mplstyle')

    return None

# -------------- End of file
