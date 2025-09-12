import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta
from tqdm import tqdm


def cams_iasi_avk(time, lat, lon, n2o_lay, cams_pre_lev, m, coord):

    utc_cams = conv_cams_time_readable(time, m)

    lon_min, lon_max = coord[0]
    lat_min, lat_max = coord[1]

    cams_avk = []
    pressure = []

    for idx_t, time_step in tqdm(enumerate(utc_cams), total=len(utc_cams)):  # create loading bar

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
            if not (lat_min <= lat_val <= lat_max):
                continue

            for idx_lon, lon_val in enumerate(lon):
                if not (lon_min <= lon_val <= lon_max):
                    continue

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

                # Replace *all* zeros in cams_n2o_lay with 10e-11
                cams_n2o_lay = np.where(cams_n2o_lay == 0.0, 10e-11, cams_n2o_lay)

                # Then assign values where the pressure condition is met
                n2o_iasi_lev[np.where(cams_pre_lay[-1] > iasi_pre)[0]] = cams_n2o_lay[-1]
                # if cams_n2o_lay[-1] == 0.0:
                #     n2o_iasi_lev[-1] = 10e-11
                # else:
                #     n2o_iasi_lev[np.where(cams_pre_lay[-1] > iasi_pre)[0]] = cams_n2o_lay[-1]

                # interpolate n2o ppm values for iasi lev. All flipped bc np.interp needs ascending values
                n2o_iasi_lev[np.where(np.isnan(n2o_iasi_lev))[0][::-1]] = np.interp(
                    iasi_pre[np.where(np.isnan(n2o_iasi_lev))[0]][::-1],  # new levels on which to interpolate
                    cams_pre_lay[::-1], cams_n2o_lay[::-1])  # cams data

                # how would iasi see the cams atmosphere
                x_sim = np.exp(np.matmul(avk, np.log(n2o_iasi_lev[::-1])) +
                                   np.matmul((np.eye(avk.shape[0]) - avk), np.log(apri[::-1])))[::-1]

                cams_avk.append(x_sim * 1000)   # back to ppb
                pressure.append(iasi_pre)

    cams_avk = pad_front_nan_29(cams_avk)
    pressure = pad_front_nan_29(pressure)

    p_min = np.nanmin(pressure)
    p_max = np.nanmax(pressure)

    # Define log-spaced pressure grid (e.g. 40 levels)
    pressure_grid = np.logspace(np.log10(p_max), np.log10(p_min), 40)

    interp_profiles = []

    print('Interpolating CAMS AVK:')
    for i in tqdm(range(cams_avk.shape[0])):
        p_profile = pressure[i, :]
        n2o_profile = cams_avk[i, :]

        mask = ~np.isnan(p_profile) & ~np.isnan(n2o_profile)
        if np.sum(mask) < 2:  # need at least 2 points to interpolate
            continue

        interp_n2o = np.interp(pressure_grid, p_profile[mask][::-1], n2o_profile[mask][::-1])
        interp_profiles.append(interp_n2o)

    interp_profiles = np.array(interp_profiles)
    mean_profile = np.nanmean(interp_profiles, axis=0)

    return mean_profile, pressure_grid


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


def read_iasi_day(month, day):
    #path = f'/data/Data/IASI_N2O/2020_{month:02d}_{day:02d}_qf3.nc'
    path = f'/misc/ghgcci7/fabian/iasi_data/2020_{month:02d}_{day:02d}_qf3.nc'
    ds = xr.open_dataset(path)
    return ds



def pad_front_nan_29(arr_list):
    """
    arr_list: iterable of 1D arrays/lists, each length <= 29
    returns: (N, 29) array, front-padded with NaN
    """
    target = 29
    out = np.full((len(arr_list), target), np.nan, dtype=float)
    for i, a in enumerate(arr_list):
        a = np.asarray(a, dtype=float)
        out[i, -a.size:] = a  # right-align -> pads front with NaN
    return out