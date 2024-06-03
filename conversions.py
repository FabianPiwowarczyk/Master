# routines for total column etc. conversions

import numpy as np

from constants import *

def hum2mr(hum_lay):
    """
    convert specific humidity (mass h2o / mass wet air) to dry air mole fraction (n_h2o / n_dryair) in ppm
    """
    h2o_lay = hum_lay / (1. - hum_lay) * mdry / mh2o * 1.e6 # reversed, TBC

    return(h2o_lay)

def gas2mr(gas_lay, mgas, h2o_dry):
    """
    convert gas profiles (mass gas / mass wet air) to dry air mole fraction (n_gas / n_dryair) in ppm

    gas_lay = input profile in kg/kg
    mgas = molar mass of gas
    h2o_dry = dry air mole fraction of h2o in ppm
    """

    hdry = h2o_dry * 1.e-6 # ppm -> 1
    
    #calculate molar mass of wet air
    h2o_wet = hdry / (1. + hdry)
    mwet = mdry * (1. - h2o_wet) + mh2o * h2o_wet

    # dry air mole fraction
    gas_dry = gas_lay * mwet/mgas * (1. + hdry) * 1.e6 # in ppm
    
    return gas_dry

def total_column(gas_lay, dry_col):
    """
    compute total column (average dry air mole fraction)
    """
    tc = np.sum(gas_lay*dry_col, axis=-1)/np.sum(dry_col, axis=-1) # last dim = altitude
    
    return tc

def grav_geoid(lat):
    """
    calculate gravity constant at geoid (NN)
    """
    gra_geo = 9.780327 * (1. + 0.0053024) * (np.sin(lat / 180. * np.pi) ** 2. - 0.0000058 * (np.sin(2. * lat / 180. * np.pi)) ** 2.)

    return gra_geo

def pressure_levels(pre_lay, pre_sur):
    """
    calculate pressure levels from layers
    """
    nsc = pre_lay.shape[0]
    ngp = pre_lay.shape[1]
    nlev = pre_lay.shape[2] + 1

    pre_lev = np.zeros((nsc,ngp,nlev),dtype='d')
    pre_lev[:,:,0] = pre_sur[:,:]

    for i in range(nlev-1):
        #pre_lev[:,:,i+1] = pre_lev[:,:,i] - 2. * (pre_lev[:,:,i] - pre_lay[:,:,i])
        pre_lev[:,:,i+1] = pre_lev[:,:,i] + 2. * (pre_lay[:,:,i] - pre_lev[:,:,i])

    return pre_lev

def geopot_levels(phi_sur, tem_lay, hum_lay, pre_lev):
    """
    calculate geopotential at levels
    """

    phi_lev = np.zeros_like(pre_lev)
    phi_lev[:,:,0] = phi_sur
    phi_lev[:,:,1:-1] = gasc_air * tem_lay[:,:,0:-1] * (1. + (gasc_h2o / gasc_air - 1.) * hum_lay[:,:,0:-1]) * (np.log(pre_lev[:,:,0:-2] / pre_lev[:,:,1:-1]))
    phi_lev = phi_lev.cumsum(axis=2)
    phi_lev[:,:,-1] = np.inf                                                                

    return phi_lev

def geopot_layers(phi_lev, tem_lay, hum_lay, pre_lev):
    """
    calculate geopotential at layers
    """

    phi_lay = np.zeros_like(tem_lay)
    phi_lay[:,:,0:-1] = phi_lev[:,:,0:-2] + gasc_air * tem_lay[:,:,0:-1] * (1. + (gasc_h2o / gasc_air - 1.) * hum_lay[:,:,0:-1]) * (1. - pre_lev[:,:,1:-1] / (pre_lev[:,:,0:-2] - pre_lev[:,:,1:-1]) * (np.log(pre_lev[:,:,0:-2]) - np.log(pre_lev[:,:,1:-1])))
    
    phi_lay[:,:,-1] = phi_lev[:,:,-2] + np.log(2.) * gasc_air * tem_lay[:,:,-1] * (1. + (gasc_h2o / gasc_air - 1.) * hum_lay[:,:,-1])

    return phi_lay
       
def geometric_height(phi_lay, gra_geo):
    """
    calculate geometric height at layers
    """

    print(gra_geo.shape, gra_geo[:, :,None].shape)
    gra3d = np.zeros_like(phi_lay)
    gra3d[:,:,:] = gra_geo[:, :,None]
    
    gmh_lay = -gra3d / fac - np.sqrt(np.clip((gra3d**2.)/(fac**2.) + (2./fac)*phi_lay,0,None))

    return gmh_lay

def grav_acc(gra_geo, gmh_lay):
    """
    calculate gravitational acceleration at layers
    """
    gra3d = np.zeros_like(gmh_lay)
    gra3d[:,:,:] = gra_geo[:, :,None]
    gra_lay = gra3d + fac * gmh_lay

    return gra_lay

def dry_column(pre_lev, h2o_dry, gra_lay):
    """
    calculate dry column
    pre_lev = pressure on levels
    h2o_dry =  dry air mole fraction of h2o in ppm
    gra_lay = gravitational acceleration at layers
    """
    
    #calculate delta pressure in Pa devided by 10^4 in order to calculate number of particles over 1cm^2 (not 1m^2)
    pre_del = (pre_lev[:,:,0:-1] - pre_lev[:,:,1:])
    pre_del /= 100.
    #convert dry mole fraction of h2o to wet mole fraction (not in ppm)
    hdry = h2o_dry * 1e-6 # ppm -> 1
    h2o_wet = hdry / (1. + hdry)
    #calculate molar mass of wet air
    mwet = mdry * (1. - h2o_wet) + mh2o * h2o_wet
    #calculate number of wet air particles per layer
    wet_col = pre_del * avogadro / mwet / gra_lay
    #calculate number of dry air particles per layer
    dry_col = wet_col * (1. - h2o_wet)

    return dry_col
