import glob2
from os.path import join
import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta
from tqdm import tqdm


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

    tc_reshape = np.full((tc.shape[0] * tc.shape[1] * tc.shape[2], 3), np.nan)

    idx_count = 0

    print(f'Reshaping cams tc{month}:')
    for t in tqdm(list(range(tc.shape[0]))):  # loading bar
        for lat in range(tc.shape[1]):
            for lon in range(tc.shape[2]):
                tc_reshape[idx_count, 0] = tc[t, lat, lon]
                tc_reshape[idx_count, 1] = longitude[lon]
                tc_reshape[idx_count, 2] = latitude[lat]
                idx_count += 1

    iasi_tc = iasi_avk_tc(time, latitude, longitude, N2O, P, month)

    df_cams = pd.DataFrame({
        'tot_col': tc_reshape[:, 0],
        'lon': tc_reshape[:, 1],
        'lat': tc_reshape[:, 2]
    })

    df_iasi = pd.DataFrame({
        'tot_col': iasi_tc[:],
        'lon': tc_reshape[:, 1],
        'lat': tc_reshape[:, 2]
    })

    return df_cams, df_iasi


def conv_iasi_time_readable(time):
    # Define the epoch start date
    epoch_start = datetime(2000, 1, 1, 0, 0, 0)

    # Vectorized calculation of datetime objects
    utc_times = np.array([epoch_start + timedelta(seconds=int(seconds)) for seconds in time])

    return utc_times


def conv_cams_time_readable(time, month):
    # Define the epoch start date
    epoch_start = datetime(2020, month, 1, 0, 0, 0)

    # Vectorized calculation of datetime objects
    utc_times = np.array([epoch_start + timedelta(hours=int(hours)) for hours in time])

    return utc_times


def iasi_avk_tc(time, lat, lon, n2o_lay, cams_pre_lev, m):

    utc_cams = conv_cams_time_readable(time, m)

    iasi_tc = np.full((time.shape[0] * lat.shape[0] * lon.shape[0]), np.nan)

    idx_count = 0
    for idx_t, time_step in tqdm(enumerate(utc_cams), total=len(utc_cams)):  # create loading bar
        # print(idx_t)
        # print(time_step)
        # if idx_t < 199:
        #     continue
        iasi_data = read_iasi_day(time_step.month, time_step.day)

        # Extract lat and lon as 1D arrays
        iasi_lat = iasi_data['lat'].values
        iasi_lon = iasi_data['lon'].values

        df_iasi = pd.DataFrame({
            'lat': iasi_lat,
            'lon': iasi_lon,
            'distance': np.nan  # Fill with NaN
        })

        for idx_lat, lat_val in enumerate(lat):
            for idx_lon, lon_val in enumerate(lon):

                # calc distance to coord point to all iasi points
                df_iasi['distance'] = haversine_distance_vectorized(df_iasi['lat'].values,
                                                                    df_iasi['lon'].values,
                                                                    lat_val,
                                                                    lon_val)
                min_idx = df_iasi['distance'].idxmin()  # index with minimal distance

                # prep iasi data
                nan_count = sum(np.isnan(iasi_data['apri'].values[min_idx, :]))
                avk = iasi_data['avk'].isel(index=min_idx).values[nan_count:, nan_count:]

                apri = iasi_data['apri'].values[min_idx, nan_count:]
                iasi_pre = iasi_data['pre_lev'].values[min_idx, nan_count:]

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
                n2o_iasi_lev[np.where(np.isnan(n2o_iasi_lev))[0][::-1]] = np.interp(
                    iasi_pre[np.where(np.isnan(n2o_iasi_lev))[0]][::-1],  # new levels on which to interpolate
                    cams_pre_lay[::-1], cams_n2o_lay[::-1])  # cams data

                # how would iasi see the cams atmosphere
                x_sim = np.exp(np.matmul(avk, np.log(n2o_iasi_lev[::-1])) +
                                   np.matmul((np.eye(avk.shape[0]) - avk), np.log(apri[::-1])))[::-1]
                x_sim_lay = lev2lay(x_sim)

                # calc cams tc as iasi would have seen it, with cams tc method
                del_P = iasi_pre[:-1] - iasi_pre[1:]
                x_sim_tc = np.sum(del_P * x_sim_lay) / iasi_pre[0]

                # calc total column as iasi does, with iasi dry column
                #x_sim_tc = total_column(x_sim_lay, iasi_dry_col)

                iasi_tc[idx_count] = x_sim_tc * 1000  # back to ppb
                idx_count += 1

    return iasi_tc


def read_iasi_day(month, day):
    #path = f'/data/Data/IASI_N2O/2020_{month:02d}_{day:02d}_qf3.nc'
    path = f'/misc/ghgcci7/fabian/iasi_datav1/2020_{month:02d}_{day:02d}_qf3.nc'
    ds = xr.open_dataset(path)
    return ds


def haversine_distance_vectorized(lat1, lon1, lat2, lon2):
    """
    Vectorized haversine distance (in km) between arrays of (lat1, lon1) and a single (lat2, lon2).
    """
    R = 6371.0  # Earth radius in km
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])

    dlat = lat2 - lat1
    dlon = lon2 - lon1

    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    a = np.clip(a, 0, 1)  # negates rounding error by python
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


def lev2lay(x_lev):
    return (x_lev[1:] + x_lev[:-1]) / 2


def total_column(gas_lay, dry_col):
    """
    compute total column (average dry air mole fraction)
    input:
        gas_lay : gas mixing ration in layers
        dry_col : dry column
    """

    tc = np.sum(gas_lay * dry_col) / np.sum(dry_col)  # last dim = altitude

    return tc