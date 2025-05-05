import glob2
from os.path import join
import netCDF4 as nc
import numpy as np
import pandas as pd


def cams_tc(month):
    dir = '/misc/ghgcci3/data/CAMS_N2O_v21r1/'
    paths = glob2.glob(join(f'{dir}', 'cams73_v21r1_n2o_conc_surface_inst_2020*.nc'))
    paths.sort()

    dataset = nc.Dataset(paths[month-1], mode='r')

    N2O = dataset.variables['N2O'][:, :, :, :]
    latitude = dataset.variables['latitude'][:]
    longitude = dataset.variables['longitude'][:]

    A = dataset.variables['ap'][:]
    B = dataset.variables['bp'][:]
    Psurf = dataset.variables['Psurf'][:, :, :]

    # Expand A and B to allow broadcasting
    A_expanded = A[None, :, None, None]  # Shape (1, lev, 1, 1)
    B_expanded = B[None, :, None, None]  # Shape (1, lev, 1, 1)

    # Expand Psurf to match
    Psurf_expanded = Psurf[:, None, :, :]  # Shape (time, 1, lat, lon)

    P = A_expanded + B_expanded * Psurf_expanded  # Shape (time, lev, lat, lon)

    del_P = P[:, :-1, :, :] - P[:, 1:, :, :]

    tc = np.sum(del_P * N2O, axis=1) / Psurf  # Shape (time, lat, lon)

    tc_reshape = np.zeros((tc.shape[0] * tc.shape[1] * tc.shape[2], 3))

    idx_count = 0
    for time in range(tc.shape[0]):
        for lat in range(tc.shape[1]):
            for lon in range(tc.shape[2]):
                tc_reshape[idx_count, 0] = tc[time, lat, lon]
                tc_reshape[idx_count, 1] = longitude[lon]
                tc_reshape[idx_count, 2] = latitude[lat]
                idx_count += 1

    df = pd.DataFrame({
        'tot_col': tc_reshape[:, 0],
        'lon': tc_reshape[:, 1],
        'lat': tc_reshape[:, 2]
    })

    return df
