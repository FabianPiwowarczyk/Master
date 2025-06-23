import glob2
from os.path import join
import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta


def cams_tc(month):
    dir = '/misc/ghgcci3/data/CAMS_N2O_v21r1/'
    paths = glob2.glob(join(f'{dir}', 'cams73_v21r1_n2o_conc_surface_inst_2020*.nc'))
    paths.sort()

    dataset = nc.Dataset(paths[month-1], mode='r')

    N2O = dataset.variables['N2O'][:, :, :, :]  # Shape (time, lay, lat, lon)
    latitude = dataset.variables['latitude'][:]
    longitude = dataset.variables['longitude'][:]
    time = dataset.variables['time'][:]

    find_nearest(time, latitude, longitude)

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


def conv_iasi_time_readable(time):
    # Define the epoch start date
    epoch_start = datetime(2000, 1, 1, 0, 0, 0)

    # Vectorized calculation of datetime objects
    utc_times = np.array([epoch_start + timedelta(seconds=int(seconds)) for seconds in time])

    return utc_times


def conv_cams_time_readable(time):
    # Define the epoch start date
    epoch_start = datetime(2020, 1, 1, 0, 0, 0)

    # Vectorized calculation of datetime objects
    utc_times = np.array([epoch_start + timedelta(hours=int(hours)) for hours in time])

    return utc_times


def find_nearest(time, lat, lon):

    utc_cams = conv_cams_time_readable(time)

    for time_step in utc_cams:
        iasi_data = read_iasi_day(time_step.month, time_step.day)

        print(iasi_data)

        import sys
        sys.exit()


def read_iasi_day(month, day):
    path = f'finished_data/2020_{month:02d}_{day:02d}_qf3.nc'
    ds = xr.open_dataset(path)
    return ds