import glob2
from os.path import join

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta

from matplotlib.pyplot import figure
from tqdm import tqdm
import matplotlib.pyplot as plt


def plot_mean_columns(month, coord):
    dir = '/misc/ghgcci3/data/CAMS_N2O_v21r1/'
    paths = glob2.glob(join(f'{dir}', 'cams73_v21r1_n2o_conc_surface_inst_2020*.nc'))
    paths.sort()

    dataset = nc.Dataset(paths[month-1], mode='r')

    N2O = dataset.variables['N2O'][:, :, :, :]  # Shape (time, lay, lat, lon)
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

    P_lay = layer_mid_pressure(P)  # Shape (time, lay, lat, lon)

    cams_prof = np.full((N2O.shape[0] * N2O.shape[2] * N2O.shape[3], N2O.shape[1] + 2), np.nan)
    cams_pressure = np.full((N2O.shape[0] * N2O.shape[2] * N2O.shape[3], N2O.shape[1] + 2), np.nan)

    print(f'Reshaping cams prof:')
    idx_count = 0
    for t in tqdm(list(range(P_lay.shape[0]))):  # loading bar
        for lat in range(P_lay.shape[2]):
            for lon in range(P_lay.shape[3]):
                cams_prof[idx_count, :-2] = N2O[t, :, lat, lon]
                cams_prof[idx_count, -2] = longitude[lon]
                cams_prof[idx_count, -1] = latitude[lat]
                cams_pressure[idx_count, :-2] = P_lay[t, :, lat, lon]
                cams_pressure[idx_count, -2] = longitude[lon]
                cams_pressure[idx_count, -1] = latitude[lat]
                idx_count += 1

    cams_data = pd.DataFrame(cams_prof)
    cams_pressure = pd.DataFrame(cams_pressure)

    del dataset
    del cams_prof
    del A
    del B
    del P
    del P_lay
    del Psurf

    cams_data = cams_data.rename(columns={cams_data.columns[-2]: "lon",
                                        cams_data.columns[-1]: "lat"})
    cams_pressure = cams_pressure.rename(columns={cams_pressure.columns[-2]: "lon",
                                        cams_pressure.columns[-1]: "lat"})

    filtered_cams = cams_data[
        (cams_data["lon"] >= coord[0][0]) & (cams_data["lon"] <= coord[0][1]) &
        (cams_data["lat"] >= coord[1][0]) & (cams_data["lat"] <= coord[1][1])
        ]

    filtered_pre = cams_pressure[
        (cams_pressure["lon"] >= coord[0][0]) & (cams_pressure["lon"] <= coord[0][1]) &
        (cams_pressure["lat"] >= coord[1][0]) & (cams_pressure["lat"] <= coord[1][1])
        ]

    p_min = np.nanmin(filtered_pre.iloc[:, :-2].to_numpy())
    p_max = np.nanmax(filtered_pre.iloc[:, :-2].to_numpy())

    # Define log-spaced pressure grid (e.g. 40 levels)
    pressure_grid_cams = np.logspace(np.log10(p_max), np.log10(p_min), 40)

    interp_profiles = []

    print('Interpolating CAMS:')
    for i in tqdm(range(filtered_cams.shape[0])):
        p_profile = filtered_pre.iloc[i, :-2].to_numpy()
        n2o_profile = filtered_cams.iloc[i, :-2].to_numpy()

        interp_n2o = np.interp(pressure_grid_cams, p_profile[::-1], n2o_profile[::-1])
        interp_profiles.append(interp_n2o)

    interp_profiles = np.array(interp_profiles)
    cams_data = np.nanmean(interp_profiles, axis=0)
    cams_pressure = pressure_grid_cams / 100

    iasi_data, iasi_pressure = read_iasi(month, coord)

    fig = figure()
    plt.plot(cams_data, np.log(cams_pressure), label='CAMS')
    plt.plot(iasi_data, np.log(iasi_pressure), label='IASI')

    plt.gca().invert_yaxis()  # invert the y-axis
    plt.xlabel('Data ppb')
    plt.ylabel('Pressure hPa')
    plt.legend()
    plt.show()


def read_iasi(month, coord):
    path = f'/misc/ghgcci7/fabian/iasi_datav1/2020_{month:02d}_*_qf3.nc'
    paths = glob2.glob(f'{path}')
    paths.sort()

    ds_iasi = xr.open_dataset(paths[0])[['n2o_lev_dry', 'pre_lev', 'lon', 'lat']]
    x_iasi = ds_iasi['n2o_lev_dry'].values
    p_iasi = ds_iasi['pre_lev'].values
    lon_iasi = ds_iasi['lon'].values
    lat_iasi = ds_iasi['lat'].values

    iasi_data = np.full((x_iasi.shape[0], 29 + 2), np.nan)
    iasi_data[:, 29-x_iasi.shape[1]:-2] = x_iasi
    iasi_data[:, -2] = lon_iasi
    iasi_data[:, -1] = lat_iasi

    iasi_data = pd.DataFrame(iasi_data)
    iasi_data = iasi_data.rename(columns={iasi_data.columns[-2]: "lon",
                                          iasi_data.columns[-1]: "lat"})

    pre_data = np.full((p_iasi.shape[0], 29 + 2), np.nan)
    pre_data[:, 29 - p_iasi.shape[1]:-2] = p_iasi
    pre_data[:, -2] = lon_iasi
    pre_data[:, -1] = lat_iasi

    pre_data = pd.DataFrame(pre_data)
    pre_data = pre_data.rename(columns={pre_data.columns[-2]: "lon",
                                          pre_data.columns[-1]: "lat"})

    filtered = iasi_data[
        (iasi_data["lon"] >= coord[0][0]) & (iasi_data["lon"] <= coord[0][1]) &
        (iasi_data["lat"] >= coord[1][0]) & (iasi_data["lat"] <= coord[1][1])
        ]

    pre_fil = pre_data[
        (pre_data["lon"] >= coord[0][0]) & (pre_data["lon"] <= coord[0][1]) &
        (pre_data["lat"] >= coord[1][0]) & (pre_data["lat"] <= coord[1][1])
        ]

    for pa in paths[1:]:
        ds_new = xr.open_dataset(pa)[['n2o_lev_dry', 'pre_lev', 'lon', 'lat']]
        new_iasi = ds_new['n2o_lev_dry'].values
        p_new = ds_new['pre_lev'].values
        lon_new = ds_new['lon'].values
        lat_new = ds_new['lat'].values

        new_data = np.full((new_iasi.shape[0], 29 + 2), np.nan)
        new_data[:, 29-new_iasi.shape[1]:-2] = new_iasi
        new_data[:, -2] = lon_new
        new_data[:, -1] = lat_new

        new_data = pd.DataFrame(new_data)
        new_data = new_data.rename(columns={new_data.columns[-2]: "lon",
                                              new_data.columns[-1]: "lat"})

        pre_new = np.full((p_new.shape[0], 29 + 2), np.nan)
        pre_new[:, 29 - p_new.shape[1]:-2] = p_new
        pre_new[:, -2] = lon_new
        pre_new[:, -1] = lat_new

        pre_new = pd.DataFrame(pre_new)
        pre_new = pre_new.rename(columns={pre_new.columns[-2]: "lon",
                                            pre_new.columns[-1]: "lat"})

        new_filtered = new_data[
            (new_data["lon"] >= coord[0][0]) & (new_data["lon"] <= coord[0][1]) &
            (new_data["lat"] >= coord[1][0]) & (new_data["lat"] <= coord[1][1])
            ]

        pre_fil_new = pre_new[
            (pre_new["lon"] >= coord[0][0]) & (pre_new["lon"] <= coord[0][1]) &
            (pre_new["lat"] >= coord[1][0]) & (pre_new["lat"] <= coord[1][1])
            ]

        filtered = pd.concat([filtered, new_filtered])
        pre_fil = pd.concat([pre_fil, pre_fil_new])

    filtered_plots = filtered.iloc[:, :-2].to_numpy()
    pre_fil_plots = pre_fil.iloc[:, :-2].to_numpy()

    fig = plt.figure()
    for i in range(filtered_plots.shape[0]):
        plt.plot(filtered_plots[i, :], np.log(pre_fil_plots[i, :]))

    plt.gca().invert_yaxis()  # invert the y-axis
    plt.xlabel('Data ppb')
    plt.ylabel('Pressure hPa')
    plt.show()

    p_min = np.nanmin(pre_fil.iloc[:, :-2].to_numpy())
    p_max = np.nanmax(pre_fil.iloc[:, :-2].to_numpy())

    # Define log-spaced pressure grid (e.g. 40 levels)
    pressure_grid = np.logspace(np.log10(p_max), np.log10(p_min), 40)

    interp_profiles = []

    print('Interpolating IASI:')
    for i in tqdm(range(filtered.shape[0])):
        p_profile = pre_fil.iloc[i, :-2].to_numpy()
        n2o_profile = filtered.iloc[i, :-2].to_numpy()

        mask = ~np.isnan(p_profile) & ~np.isnan(n2o_profile)
        if np.sum(mask) < 2:  # need at least 2 points to interpolate
            continue

        interp_n2o = np.interp(pressure_grid, p_profile[mask][::-1], n2o_profile[mask][::-1])
        interp_profiles.append(interp_n2o)

    interp_profiles = np.array(interp_profiles)
    mean_profile = np.nanmean(interp_profiles, axis=0) * 1000

    return mean_profile, pressure_grid


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
    P_top = P[:, :-1, :, :]  # Shape (time, lev, lat, lon)
    P_bottom = P[:, 1:, :, :]

    # Geometric/logarithmic mean (exponential profile assumption)
    P_mid = np.exp(0.5 * (np.log(P_top) + np.log(P_bottom)))

    return P_mid