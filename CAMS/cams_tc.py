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
    for t in range(tc.shape[0]):
        for lat in range(tc.shape[1]):
            for lon in range(tc.shape[2]):
                tc_reshape[idx_count, 0] = tc[t, lat, lon]
                tc_reshape[idx_count, 1] = longitude[lon]
                tc_reshape[idx_count, 2] = latitude[lat]
                idx_count += 1

    df = pd.DataFrame({
        'tot_col': tc_reshape[:, 0],
        'lon': tc_reshape[:, 1],
        'lat': tc_reshape[:, 2]
    })

    iasi_avk_tc(time, latitude, longitude, N2O, P)

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


def iasi_avk_tc(time, lat, lon, n2o_lay, cams_pre_lev):

    utc_cams = conv_cams_time_readable(time)

    for idx_t, time_step in enumerate(utc_cams):
        iasi_data = read_iasi_day(time_step.month, time_step.day)

        # Extract lat and lon as 1D arrays
        iasi_lat = iasi_data['lat'].values
        iasi_lon = iasi_data['lon'].values

        df_iasi = pd.DataFrame({
            'lat': iasi_lat,
            'lon': iasi_lon,
            'distance': np.nan  # Fill with NaN
        })

        for idx_lat, lat in enumerate(lat):
            for idx_lon, lon in enumerate(lon):
                coord = (lat, lon)

                # calc distance to coord point to all iasi points
                df_iasi['distance'] = df_iasi.apply(lambda row: haversine_distance((row['lat'], row['lon']),
                                                                                   coord), axis=1)
                min_idx = df_iasi['distance'].idxmin()  # index with minimal distance

                # prep iasi data
                nan_count = sum(np.isnan(iasi_data['apri'].values[min_idx, :]))
                avk = iasi_data['avk'].values[min_idx, nan_count:, nan_count:]
                apri = iasi_data['apri'].values[min_idx, nan_count:]
                iasi_pre = iasi_data['pre_lev'].values[min_idx, nan_count:]
                iasi_dry_col = iasi_data['dry_col'].values[min_idx, nan_count:]

                # prep cams N2O and pressure
                cams_pre_lay = layer_mid_pressure(cams_pre_lev[idx_t, :, idx_lat, idx_lon])
                cams_pre_lay /= 100  # cams to hPa to fit iasi
                cams_n2o_lay = n2o_lay[idx_t, :, idx_lat, idx_lon]
                cams_n2o_lay /= 1000  # cams to ppm to fit iasi

                # prep array for cams n2o on iasi levels
                n2o_iasi_lev = np.full((iasi_pre.shape[0]), np.nan)
                # edge cases for iasi and cams ground and top of atmos. not aligning
                n2o_iasi_lev[np.where(cams_pre_lay[0] < iasi_pre)[0]] = cams_n2o_lay[0]
                n2o_iasi_lev[np.where(cams_pre_lay[-1] > iasi_pre)[0]] = cams_n2o_lay[-1]


                # interpolate n2o ppm values for iasi lev. All flipped bc np.interp needs ascending values
                #gosat_n2o_lev[1:-1] = np.interp(pre_lvl[1:-1][::-1], pre_lay[::-1], n2o_gosat_lay[::-1])[::-1]

                import matplotlib.pyplot as plt

                plt.scatter(cams_n2o_lay, cams_pre_lay)
                plt.scatter(apri, iasi_pre)
                plt.show()

                import sys
                sys.exit()


def read_iasi_day(month, day):
    path = f'finished_data/2020_{month:02d}_{day:02d}_qf3.nc'
    ds = xr.open_dataset(path)
    return ds


def haversine_distance(coord1, coord2):
    """
    Calculate the great-circle distance between two (lat, lon) points in kilometers.

    Parameters:
        coord1 (tuple): (latitude, longitude) in degrees
        coord2 (tuple): (latitude, longitude) in degrees

    Returns:
        float: distance in kilometers
    """
    # Earth radius in kilometers
    R = 6371.0

    lat1, lon1 = np.radians(coord1)
    lat2, lon2 = np.radians(coord2)

    dlat = lat2 - lat1
    dlon = lon2 - lon1

    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arcsin(np.sqrt(a))

    return R * c


def layer_mid_pressure(P):
    """
    Compute pressure at the center of each layer, assuming pressure drops exponentially with height
    and constant temperature in a layer.

    Parameters:
        P (array-like): Pressure at level boundaries (length n+1)

    Returns:
        np.ndarray: Pressure at layer centers (length n)
    """

    P = np.where(P == 0, 1e-6, P)  # don't let a log(0) stop you!
    P_top = P[:-1]
    P_bottom = P[1:]

    # Geometric/logarithmic mean (exponential profile assumption)
    P_mid = np.exp(0.5 * (np.log(P_top) + np.log(P_bottom)))

    return P_mid