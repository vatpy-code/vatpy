# Description: Functions to get various gas properties, e.g. the number density
#              of different chemical species, as well as gas temperature &
#              thermal pressure.
# Authour(s): Jonathan Petersson
# Last updated: 2025-09-10


# -------------- Required packages
import numpy as np


# -------------- Declare function(s)
def number_density(h, iu):
    # Constants:
    mp = 1.6726e-24  # [g]

    # Gas density:
    rho = h['PartType0']['Density'] * iu['udens']

    # Chemical abundances:
    X = h['PartType0']['ChemicalAbundances'][:]
    xH2, xHII, xCO = X[:, 0], X[:, 1], X[:, 2]

    # Number densities from chemical abundances:
    xHe = 0.1
    mu = (1 + 4 * xHe)
    nH = rho / (mu * mp)

    num = {
        'HII': xHII * nH,
        'H2': xH2 * nH,
        'HI': (1 - xHII - 2*xH2) * nH,
        'CO': xCO * nH,
        'He': xHe * nH,
        'e': xHII * nH
    }

    # Total number density:
    # num['n'] = np.sum(np.array(list(num.values())), axis=0)
    num['n'] = num['HI'] + num['HII'] + num['e'] + num['H2'] + num['He']

    return num


def temperature(h, iu):
    mp = 1.6726e-24  # [g]
    kb = 1.380658e-16  # [erg K^-1]

    if 'ChemicalAbundances' in h['PartType0'].keys():
        # Convert mass densities into number densities:
        num = number_density(h, iu)

        # Estimate the temperature of each gas cell:
        interg = h['PartType0']['InternalEnergy'] * iu['uinterg']
        mu = ((num['HII'] * 1 + num['HI'] * 1 + num['H2'] * 2 + num['He'] * 4 +
               num['CO'] * 14) / num['n'])
        temp = (2 * mu * mp * interg) / (3 * kb)
    else:
        xHe = 0.1
        mu = (1 + 4 * xHe)
        interg = h['PartType0']['InternalEnergy'] * iu['uinterg']
        temp = (2 * mu * mp * interg) / (3 * kb)

    return temp


def thermalpressure(h, iu):
    # Gas density:
    rho = h['PartType0']['Density'] * iu['udens']

    # Internal energy:
    interg = h['PartType0']['InternalEnergy'] * iu['uinterg']

    # Thermal pressure:
    gamma = 5 / 3
    cs2 = gamma * (gamma - 1) * interg
    Ptherm = cs2 / gamma * rho

    return Ptherm

# -------------- End of file
