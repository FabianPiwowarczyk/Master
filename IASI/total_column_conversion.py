import numpy as np

from constants import *


def total_column(gas_lay, dry_col):
    """
    compute total column (average dry air mole fraction)
    input:
        gas_lay : ?
        dry_col : dry column
    """
    tc = np.sum(gas_lay * dry_col, axis=-1) / np.sum(dry_col, axis=-1)  # last dim = altitude

    return tc


def dry_column(pre_lev, h2o_dry, gra_lay):
    """
    calculate dry column
    input:
        pre_lev : pressure on levels
        h2o_dry : dry air mole fraction of h2o in ppm
        gra_lay : gravitational acceleration at layers
    """

    # calculate delta pressure in Pa divided by 10^4 in order to calculate number of particles over 1cm^2 (not 1m^2)
    pre_del = (pre_lev[:, :, 0:-1] - pre_lev[:, :, 1:])
    pre_del /= 100.
    # convert dry mole fraction of h2o to wet mole fraction (not in ppm)
    hdry = h2o_dry * 1e-6  # ppm -> 1
    h2o_wet = hdry / (1. + hdry)
    # calculate molar mass of wet air
    mwet = mdry * (1. - h2o_wet) + mh2o * h2o_wet
    # calculate number of wet air particles per layer
    wet_col = pre_del * avogadro / mwet / gra_lay
    # calculate number of dry air particles per layer
    dry_col = wet_col * (1. - h2o_wet)

    return dry_col


def grav_acc(gra_geo, gmh_lay):
    """
    output:
        gra_lay : gravitational acceleration at layers
    """
    gra3d = np.zeros_like(gmh_lay)
    gra3d[:,:,:] = gra_geo[:, :,None]
    gra_lay = gra3d + fac * gmh_lay

    return gra_lay


def geometric_height(phi_lay, gra_geo):
    """
    output:
        gmh_lay : geometric height at layers
    """

    print(gra_geo.shape, gra_geo[:, :, None].shape)
    gra3d = np.zeros_like(phi_lay)
    gra3d[:, :, :] = gra_geo[:, :, None]

    gmh_lay = -gra3d / fac - np.sqrt(np.clip((gra3d ** 2.) / (fac ** 2.) + (2. / fac) * phi_lay, 0, None))

    return gmh_lay


def grav_geoid(lat):
    """
    input:
        lat : latitude
    output:
        gra_geo : gravity constant at geoid (NN)
    """
    gra_geo = 9.780327 * (1. + 0.0053024) * (np.sin(lat / 180. * np.pi) ** 2. - 0.0000058 * (np.sin(2. * lat / 180. * np.pi)) ** 2.)

    return gra_geo


