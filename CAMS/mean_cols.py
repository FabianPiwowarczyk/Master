import glob2
from os.path import join
import os
import matplotlib as mpl

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta

from matplotlib.pyplot import figure
from tqdm import tqdm
import matplotlib.pyplot as plt
import calendar

from .cams_iasi_avk import cams_iasi_avk


# ---------- helpers for caching ----------
def _coord_tag(coord):
    (lon0, lon1), (lat0, lat1) = coord
    def s(v): return str(v).replace('.', 'p').replace('-', 'm')
    return f"lon_{s(lon0)}_{s(lon1)}__lat_{s(lat0)}_{s(lat1)}"


def _vec_path(varname, month, coord, outdir="derived/mean_columns"):
    os.makedirs(outdir, exist_ok=True)
    return os.path.join(outdir, f"{varname}__m{month:02d}__{_coord_tag(coord)}.txt")


def save_vec(varname, month, coord, vec):
    path = _vec_path(varname, month, coord)
    np.savetxt(path, np.asarray(vec), fmt="%.10g", header=varname)
    return path


def load_vec(varname, month, coord):
    path = _vec_path(varname, month, coord)
    if os.path.exists(path):
        return np.loadtxt(path, ndmin=1)
    return None
# ----------------------------------------


def plot_mean_columns_multi(months, coord, idx):
    """
    Plot CAMS / IASI / IASI a-priori / CAMS (IASI AVK) profiles
    for a list of months side by side. Caches each 1D array to TXT.
    """

    mpl.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.size": 10.0,  # <- from FONTSIZE in your log
        "axes.labelsize": 9.0,
        "axes.titlesize": 10.95,  # actually does matter
        "legend.fontsize": 8,  # ~\footnotesize
        "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage{lmodern}",
        # add packages you actually use in labels, e.g.:
        # r"\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{siunitx}\usepackage{mhchem}"
    })

    # Ensure months is an iterable of ints
    if isinstance(months, int):
        months = [months]

    print(f'Plot {idx}')

    ncols = len(months)
    FIGSIZE = set_size(width_fraction=0.33 * ncols, height_fraction=0.35)
    data_dir = '/misc/ghgcci3/data/CAMS_N2O_v21r1/'
    paths = glob2.glob(join(data_dir, 'cams73_v21r1_n2o_conc_surface_inst_2020*.nc'))
    paths.sort()

    fig, axes = plt.subplots(1, ncols, figsize=FIGSIZE, sharey=True)
    if ncols == 1:
        axes = np.array([axes])

    for ax, month in zip(axes, months):

        # -------- try cache first --------
        cams_data = load_vec('cams_data', month, coord)
        cams_pressure = load_vec('cams_pressure', month, coord)
        iasi_data = load_vec('iasi_data', month, coord)
        iasi_pressure = load_vec('iasi_pressure', month, coord)
        mean_apri = load_vec('mean_apri', month, coord)
        cams_data_avk = load_vec('cams_data_avk', month, coord)
        cams_avk_pressure_grid = load_vec('cams_avk_pressure_grid', month, coord)

        need_cams = cams_data is None or cams_pressure is None
        need_iasi = (iasi_data is None) or (iasi_pressure is None) or (mean_apri is None)
        need_cams_avk = (cams_data_avk is None) or (cams_avk_pressure_grid is None)

        # -------- compute what’s missing --------
        if need_cams or need_cams_avk:
            ds = nc.Dataset(paths[month-1], mode='r')

            N2O = ds.variables['N2O'][:]       # (time, lay, lat, lon)
            latitude = ds.variables['latitude'][:]
            longitude = ds.variables['longitude'][:]
            time = ds.variables['time'][:]
            A = ds.variables['ap'][:]
            B = ds.variables['bp'][:]
            Psurf = ds.variables['Psurf'][:]

            # Hybrid levels -> pressure
            A_exp = A[None, :, None, None]
            B_exp = B[None, :, None, None]
            Psurf_exp = Psurf[:, None, :, :]
            P = A_exp + B_exp * Psurf_exp              # (time, lev, lat, lon)
            P_lay = layer_mid_pressure(P)              # (time, lay, lat, lon)

            del A, B, Psurf, A_exp, B_exp, Psurf_exp

            # Reshape into rows with lon/lat appended
            nrows = N2O.shape[0] * N2O.shape[2] * N2O.shape[3]
            cams_prof = np.full((nrows, N2O.shape[1] + 2), np.nan)
            cams_pres = np.full((nrows, N2O.shape[1] + 2), np.nan)

            idx = 0
            for t in range(P_lay.shape[0]):
                for ilat in range(P_lay.shape[2]):
                    for ilon in range(P_lay.shape[3]):
                        cams_prof[idx, :-2] = N2O[t, :, ilat, ilon]
                        cams_prof[idx, -2] = longitude[ilon]
                        cams_prof[idx, -1] = latitude[ilat]
                        cams_pres[idx, :-2] = P_lay[t, :, ilat, ilon]
                        cams_pres[idx, -2] = longitude[ilon]
                        cams_pres[idx, -1] = latitude[ilat]
                        idx += 1

            df = pd.DataFrame(cams_prof)
            dfP = pd.DataFrame(cams_pres)
            df = df.rename(columns={df.columns[-2]: "lon", df.columns[-1]: "lat"})
            dfP = dfP.rename(columns={dfP.columns[-2]: "lon", dfP.columns[-1]: "lat"})

            in_box = (
                (df["lon"] >= coord[0][0]) & (df["lon"] <= coord[0][1]) &
                (df["lat"] >= coord[1][0]) & (df["lat"] <= coord[1][1])
            )
            f_df, f_dfP = df[in_box], dfP[in_box]

            # CAMS interp grid (Pa -> later hPa)
            p_min = np.nanmin(f_dfP.iloc[:, :-2].to_numpy())
            p_max = np.nanmax(f_dfP.iloc[:, :-2].to_numpy())
            pressure_grid_cams = np.logspace(np.log10(p_max), np.log10(p_min), 40)

            if need_cams:
                interp_profiles = []
                for i in range(f_df.shape[0]):
                    p_prof = f_dfP.iloc[i, :-2].to_numpy()
                    n2o_prof = f_df.iloc[i, :-2].to_numpy()
                    interp_n2o = np.interp(pressure_grid_cams, p_prof[::-1], n2o_prof[::-1])
                    interp_profiles.append(interp_n2o)
                interp_profiles = np.asarray(interp_profiles)
                cams_data = np.nanmean(interp_profiles, axis=0)
                cams_pressure = pressure_grid_cams / 100.0  # hPa
                save_vec('cams_data', month, coord, cams_data)
                save_vec('cams_pressure', month, coord, cams_pressure)

                del interp_profiles, f_df, f_dfP, P_lay, in_box, df, dfP

            if need_cams_avk:
                cams_data_avk, cams_avk_pressure_grid = cams_iasi_avk(
                    time, latitude, longitude, N2O, P, month, coord
                )
                save_vec('cams_data_avk', month, coord, cams_data_avk)
                save_vec('cams_avk_pressure_grid', month, coord, cams_avk_pressure_grid)

                del time, latitude, longitude, N2O, P
            ds.close()

        if need_iasi:
            iasi_data, mean_apri, iasi_pressure = read_iasi(month, coord)
            save_vec('iasi_data', month, coord, iasi_data)
            save_vec('mean_apri', month, coord, mean_apri)
            save_vec('iasi_pressure', month, coord, iasi_pressure)

        # reload any still-None (from cache)
        if cams_data is None: cams_data = load_vec('cams_data', month, coord)
        if cams_pressure is None: cams_pressure = load_vec('cams_pressure', month, coord)
        if cams_data_avk is None: cams_data_avk = load_vec('cams_data_avk', month, coord)
        if cams_avk_pressure_grid is None:
            cams_avk_pressure_grid = load_vec('cams_avk_pressure_grid', month, coord)
        if iasi_data is None: iasi_data = load_vec('iasi_data', month, coord)
        if mean_apri is None: mean_apri = load_vec('mean_apri', month, coord)
        if iasi_pressure is None: iasi_pressure = load_vec('iasi_pressure', month, coord)

        # -------- plot this month --------
        ax.plot(cams_data,         cams_pressure,          label='CAMS')
        ax.plot(iasi_data,         iasi_pressure,          label='IASI')
        ax.plot(mean_apri,         iasi_pressure,          label='IASI a-priori')
        ax.plot(cams_data_avk,     cams_avk_pressure_grid, label='CAMS (IASI AVK)')

        ax.invert_yaxis()
        ax.set_xlabel('N$_2$O (ppb)')
        ax.set_title(calendar.month_name[month])  # e.g., "February"

    # shared y label on the leftmost axis
    axes[0].set_ylabel('Pressure (hPa)')
    # put a single legend outside (right) or below
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.03), frameon=False)

    fig.tight_layout()
    outpath = f'pictures/mean_cols_{idx}.png'
    plt.savefig(outpath)

# def plot_mean_columns(month, coord):
#     dir = '/misc/ghgcci3/data/CAMS_N2O_v21r1/'
#     paths = glob2.glob(join(f'{dir}', 'cams73_v21r1_n2o_conc_surface_inst_2020*.nc'))
#     paths.sort()
#
#     dataset = nc.Dataset(paths[month-1], mode='r')
#
#     N2O = dataset.variables['N2O'][:, :, :, :]  # Shape (time, lay, lat, lon)
#     latitude = dataset.variables['latitude'][:]
#     longitude = dataset.variables['longitude'][:]
#     time = dataset.variables['time'][:]
#
#     A = dataset.variables['ap'][:]
#     B = dataset.variables['bp'][:]
#     Psurf = dataset.variables['Psurf'][:, :, :]
#
#     # Expand A and B to allow broadcasting
#     A_expanded = A[None, :, None, None]  # Shape (1, lev, 1, 1)
#     B_expanded = B[None, :, None, None]  # Shape (1, lev, 1, 1)
#
#     # Expand Psurf to match
#     Psurf_expanded = Psurf[:, None, :, :]  # Shape (time, 1, lat, lon)
#
#     P = A_expanded + B_expanded * Psurf_expanded  # Shape (time, lev, lat, lon)
#
#     P_lay = layer_mid_pressure(P)  # Shape (time, lay, lat, lon)
#
#     cams_prof = np.full((N2O.shape[0] * N2O.shape[2] * N2O.shape[3], N2O.shape[1] + 2), np.nan)
#     cams_pressure = np.full((N2O.shape[0] * N2O.shape[2] * N2O.shape[3], N2O.shape[1] + 2), np.nan)
#
#     print(f'Reshaping cams prof:')
#     idx_count = 0
#     for t in tqdm(list(range(P_lay.shape[0]))):  # loading bar
#         for lat in range(P_lay.shape[2]):
#             for lon in range(P_lay.shape[3]):
#                 cams_prof[idx_count, :-2] = N2O[t, :, lat, lon]
#                 cams_prof[idx_count, -2] = longitude[lon]
#                 cams_prof[idx_count, -1] = latitude[lat]
#                 cams_pressure[idx_count, :-2] = P_lay[t, :, lat, lon]
#                 cams_pressure[idx_count, -2] = longitude[lon]
#                 cams_pressure[idx_count, -1] = latitude[lat]
#                 idx_count += 1
#
#     cams_data = pd.DataFrame(cams_prof)
#     cams_pressure = pd.DataFrame(cams_pressure)
#
#     cams_data_avk, cams_avk_pressure_grid = cams_iasi_avk(time, latitude, longitude, N2O, P, month, coord)
#
#     del dataset
#     del cams_prof
#     del A
#     del B
#     del P
#     del P_lay
#     del Psurf
#
#     cams_data = cams_data.rename(columns={cams_data.columns[-2]: "lon",
#                                         cams_data.columns[-1]: "lat"})
#     cams_pressure = cams_pressure.rename(columns={cams_pressure.columns[-2]: "lon",
#                                         cams_pressure.columns[-1]: "lat"})
#
#     filtered_cams = cams_data[
#         (cams_data["lon"] >= coord[0][0]) & (cams_data["lon"] <= coord[0][1]) &
#         (cams_data["lat"] >= coord[1][0]) & (cams_data["lat"] <= coord[1][1])
#         ]
#
#     filtered_pre = cams_pressure[
#         (cams_pressure["lon"] >= coord[0][0]) & (cams_pressure["lon"] <= coord[0][1]) &
#         (cams_pressure["lat"] >= coord[1][0]) & (cams_pressure["lat"] <= coord[1][1])
#         ]
#
#     p_min = np.nanmin(filtered_pre.iloc[:, :-2].to_numpy())
#     p_max = np.nanmax(filtered_pre.iloc[:, :-2].to_numpy())
#
#     # Define log-spaced pressure grid (e.g. 40 levels)
#     pressure_grid_cams = np.logspace(np.log10(p_max), np.log10(p_min), 40)
#
#     interp_profiles = []
#
#     print('Interpolating CAMS:')
#     for i in tqdm(range(filtered_cams.shape[0])):
#         p_profile = filtered_pre.iloc[i, :-2].to_numpy()
#         n2o_profile = filtered_cams.iloc[i, :-2].to_numpy()
#
#         interp_n2o = np.interp(pressure_grid_cams, p_profile[::-1], n2o_profile[::-1])
#         interp_profiles.append(interp_n2o)
#
#     interp_profiles = np.array(interp_profiles)
#     cams_data = np.nanmean(interp_profiles, axis=0)
#     cams_pressure = pressure_grid_cams / 100
#
#     iasi_data, mean_apri, iasi_pressure = read_iasi(month, coord)
#
#     fig = figure()
#     plt.plot(cams_data, cams_pressure, label='CAMS')
#     plt.plot(iasi_data, iasi_pressure, label='IASI')
#     plt.plot(mean_apri, iasi_pressure, label='IASI a-priori')
#     plt.plot(cams_data_avk, cams_avk_pressure_grid, label='CAMS (IASI AVK)')
#
#     plt.gca().invert_yaxis()  # invert the y-axis
#     plt.xlabel('Data ppb')
#     plt.ylabel('Pressure hPa')
#     plt.legend()
#     plt.show()


def read_iasi(month, coord):
    path = f'/misc/ghgcci7/fabian/iasi_data/2020_{month:02d}_*_qf3.nc'
    paths = glob2.glob(f'{path}')
    paths.sort()

    ds_iasi = xr.open_dataset(paths[0])[['n2o_lev_dry', 'pre_lev', 'lon', 'lat', 'apri']]
    x_iasi = ds_iasi['n2o_lev_dry'].values
    x_apri = ds_iasi['apri'].values
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

    apri_data = np.full((x_apri.shape[0], 29 + 2), np.nan)
    apri_data[:, 29 - x_apri.shape[1]:-2] = x_apri
    apri_data[:, -2] = lon_iasi
    apri_data[:, -1] = lat_iasi

    apri_data = pd.DataFrame(apri_data)
    apri_data = apri_data.rename(columns={apri_data.columns[-2]: "lon",
                                          apri_data.columns[-1]: "lat"})

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

    apri_fil = apri_data[
        (apri_data["lon"] >= coord[0][0]) & (apri_data["lon"] <= coord[0][1]) &
        (apri_data["lat"] >= coord[1][0]) & (apri_data["lat"] <= coord[1][1])
        ]

    pre_fil = pre_data[
        (pre_data["lon"] >= coord[0][0]) & (pre_data["lon"] <= coord[0][1]) &
        (pre_data["lat"] >= coord[1][0]) & (pre_data["lat"] <= coord[1][1])
        ]

    del ds_iasi, x_iasi, x_apri, p_iasi, lon_iasi, lat_iasi, iasi_data, apri_data, pre_data

    for pa in paths[1:]:
        ds_new = xr.open_dataset(pa)[['n2o_lev_dry', 'pre_lev', 'lon', 'lat', 'apri']]
        new_iasi = ds_new['n2o_lev_dry'].values
        new_apri = ds_new['apri'].values
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

        new_apri_data = np.full((new_apri.shape[0], 29 + 2), np.nan)
        new_apri_data[:, 29 - new_apri.shape[1]:-2] = new_apri
        new_apri_data[:, -2] = lon_new
        new_apri_data[:, -1] = lat_new

        new_apri_data = pd.DataFrame(new_apri_data)
        new_apri_data = new_apri_data.rename(columns={new_apri_data.columns[-2]: "lon",
                                            new_apri_data.columns[-1]: "lat"})

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

        apri_fil_new = new_apri_data[
            (new_apri_data["lon"] >= coord[0][0]) & (new_apri_data["lon"] <= coord[0][1]) &
            (new_apri_data["lat"] >= coord[1][0]) & (new_apri_data["lat"] <= coord[1][1])
            ]

        pre_fil_new = pre_new[
            (pre_new["lon"] >= coord[0][0]) & (pre_new["lon"] <= coord[0][1]) &
            (pre_new["lat"] >= coord[1][0]) & (pre_new["lat"] <= coord[1][1])
            ]

        filtered = pd.concat([filtered, new_filtered])
        apri_fil = pd.concat([apri_fil, apri_fil_new])
        pre_fil = pd.concat([pre_fil, pre_fil_new])

        del ds_new, new_iasi, new_apri, p_new, new_filtered, apri_fil_new, pre_fil_new

    # filtered_plots = filtered.iloc[:, :-2].to_numpy()
    # pre_fil_plots = pre_fil.iloc[:, :-2].to_numpy()
    #
    # fig = plt.figure()
    # for i in range(filtered_plots.shape[0]):
    #     plt.plot(filtered_plots[i, :], np.log(pre_fil_plots[i, :]))
    #
    # plt.gca().invert_yaxis()  # invert the y-axis
    # plt.xlabel('Data ppb')
    # plt.ylabel('Pressure hPa')
    # plt.show()

    p_min = np.nanmin(pre_fil.iloc[:, :-2].to_numpy())
    p_max = np.nanmax(pre_fil.iloc[:, :-2].to_numpy())

    # Define log-spaced pressure grid (e.g. 40 levels)
    pressure_grid = np.logspace(np.log10(p_max), np.log10(p_min), 40)

    interp_profiles = []
    interp_apris = []

    print('Interpolating IASI:')
    for i in tqdm(range(filtered.shape[0])):
        p_profile = pre_fil.iloc[i, :-2].to_numpy()
        n2o_profile = filtered.iloc[i, :-2].to_numpy()
        apri_profile = apri_fil.iloc[i, :-2].to_numpy()

        mask = ~np.isnan(p_profile) & ~np.isnan(n2o_profile)
        if np.sum(mask) < 2:  # need at least 2 points to interpolate
            continue

        interp_n2o = np.interp(pressure_grid, p_profile[mask][::-1], n2o_profile[mask][::-1])
        interp_apri = np.interp(pressure_grid, p_profile[mask][::-1], apri_profile[mask][::-1])
        interp_profiles.append(interp_n2o)
        interp_apris.append(interp_apri)

    interp_profiles = np.array(interp_profiles)
    interp_apris = np.array(interp_apris)
    mean_profile = np.nanmean(interp_profiles, axis=0) * 1000
    mean_apri = np.nanmean(interp_apris, axis=0) * 1000

    del interp_profiles, interp_apris, filtered, apri_fil, pre_fil

    return mean_profile, mean_apri, pressure_grid


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


def set_size(textwidth_pt=466.62503, textheight_pt=674.33032, width_fraction=1.0, height_fraction=0.4):
    """
    Convert LaTeX text sizes (pt) to a Matplotlib figure size (inches).

    * \textwidth=466.62503pt * \textheight=674.33032pt

    Parameters
    ----------
    textwidth_pt : float   # \textwidth  in pt (from LaTeX log)
    textheight_pt: float   # \textheight in pt (from LaTeX log)
    width_fraction : float # e.g. 1.0 for full width, 0.5 for half width
    height_fraction: float # e.g. 0.4 for 40% of \textheight

    Returns
    -------
    (width_in, height_in) in inches
    TEXTWIDTH_PT  = 418.25555
    TEXTHEIGHT_PT = 595.80026
    """
    inches_per_pt = 1/72.27
    width_in  = textwidth_pt  * inches_per_pt * width_fraction
    height_in = textheight_pt * inches_per_pt * height_fraction
    return (width_in, height_in)