import numpy as np
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
import glob2
import netCDF4 as nc
from datetime import datetime, timedelta

from .total_column_conversion import dry_column, lev2lay, grav_geoid, grav_acc
from.chng_prior import layer_mid_pressure

from os.path import join

from .conv_constants import *


def _coord_tag(coord):
    # coord = [(lon_min, lon_max), (lat_min, lat_max)]
    (lon0, lon1), (lat0, lat1) = coord
    # make filenames safe (e.g., minus signs)
    def s(v):
        return str(v).replace('.', 'p').replace('-', 'm')
    return f"lon_{s(lon0)}_{s(lon1)}__lat_{s(lat0)}_{s(lat1)}"


def _outfile(month, coord, outdir="derived/apri_means"):
    os.makedirs(outdir, exist_ok=True)
    return os.path.join(outdir, f"apri_means_m{month:02d}__{_coord_tag(coord)}.txt")


def _save_txt(path, pressure_hpa, mean_iasi_ppb, mean_gosat_ppb):
    arr = np.column_stack([pressure_hpa, mean_iasi_ppb, mean_gosat_ppb])
    header = "pressure_hPa mean_iasi_ppb mean_gosat_ppb"
    np.savetxt(path, arr, fmt="%.6f", header=header)


def _load_txt(path):
    arr = np.loadtxt(path, ndmin=2)
    # handle case of single row
    pressure_hpa = arr[:, 0]
    mean_iasi_ppb = arr[:, 1]
    mean_gosat_ppb = arr[:, 2]
    return pressure_hpa, mean_iasi_ppb, mean_gosat_ppb


def mean_apris():
    data_path = '/misc/hypatia/data/IASI/L2_MUSICA/V3.30/MetopA/2020/{:02d}/'
    months = [2, 7, 11]  # columns
    coords = [
        [(-105, -100), (35, 40)],
        [(35, 40), (-2.5, 2.5)],
        [(135, 140), (-35, -30)]
    ]  # rows
    quality_flag = 3

    # Prepare figure
    fig, axes = plt.subplots(len(coords), len(months), figsize=(15, 12), sharex=True, sharey=True)

    for j, month in enumerate(months):  # loop over months (columns)
        month_path = data_path.format(month)
        days_paths = [f.path for f in os.scandir(month_path) if f.is_dir()]
        days_paths.sort()

        for i, coord in enumerate(coords):  # loop over coords (rows)
            # ---- Load from TXT if available ----
            out_txt = _outfile(month, coord)
            if os.path.exists(out_txt):
                print(out_txt, 'Already exists.')
                pressure_grid, mean_iasi, mean_gosat = _load_txt(out_txt)
            else:
                # ---- Collect and compute means ----
                iasi_profiles, gosat_profiles, pre_lev = iasi_day_profiles(days_paths[0], quality_flag, coord)
                for dir in days_paths[1:]:
                    new_iasi, new_gosat, new_pre = iasi_day_profiles(dir, quality_flag, coord)
                    iasi_profiles = np.vstack([iasi_profiles, new_iasi])
                    gosat_profiles = np.vstack([gosat_profiles, new_gosat])
                    pre_lev = np.vstack([pre_lev, new_pre])

                # Interpolate and average
                mean_iasi, pressure_grid = inter_profiles(iasi_profiles, pre_lev)
                mean_gosat, _ = inter_profiles(gosat_profiles, pre_lev)

                # ---- Save to TXT for future runs ----
                _save_txt(out_txt, pressure_grid, mean_iasi, mean_gosat)

            # ---- Plot into subplot ----
            ax = axes[i, j]
            ax.plot(mean_iasi, pressure_grid, label="IASI")
            ax.plot(mean_gosat, pressure_grid, label="GOSAT")
            ax.invert_yaxis()

            # ============================== ADDED ↓↓↓
            # Annotation: difference of the FIRST values (GOSAT − IASI).
            # Falls back to the first finite pair if index 0 is NaN.
            idx0 = 0
            diff_0 = mean_gosat[idx0] - mean_iasi[idx0]
            note = ""
            if not np.isfinite(diff_0):
                valid = np.where(np.isfinite(mean_gosat) & np.isfinite(mean_iasi))[0]
                if valid.size:
                    idx0 = valid[0]
                    diff_0 = mean_gosat[idx0] - mean_iasi[idx0]
                    note = f" (first finite idx={idx0})"
                else:
                    note = " (no finite overlap)"
            # bottom-left
            ax.text(
                0.03, 0.05,  # << moved near bottom-left
                f"Δ(GOSAT-2 minus IASI) ground = {diff_0:.1f} ppb{note}",
                transform=ax.transAxes, ha='left', va='bottom',
                fontsize=9, bbox=dict(boxstyle='round', fc='white', alpha=0.7, ec='none')
            )
            # (optional) for absolute difference instead, use:
            # diff_0 = np.abs(mean_gosat[idx0] - mean_iasi[idx0])
            # ============================== ADDED ↑↑↑

            # Label axes
            ax.set_xlabel("N₂O (ppb)")
            if j == 0:
                lon_rng, lat_rng = coord
                ax.set_ylabel(f"Lon: {lon_rng[0]}–{lon_rng[1]}\nLat: {lat_rng[0]}–{lat_rng[1]}\nPressure (hPa)")

            # Titles for months
            if i == 0:
                ax.set_title(f"Month {month}")

    # Legend only once (outside grid)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right")
    fig.tight_layout()

    outpath = 'pictures/apri_vergleich_iasi_gosat.png'
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.savefig(outpath)


def conv_iasi_time_readable(time):
    # Define the epoch start date
    epoch_start = datetime(2000, 1, 1, 0, 0, 0)

    # Vectorized calculation of datetime objects
    utc_times = np.array([epoch_start + timedelta(seconds=int(seconds)) for seconds in time])

    return utc_times


def read_iasi(path):
    data = dict()

    dataset = nc.Dataset(path, mode='r')

    data['n2o_lev_dry'] = dataset.variables['musica_ghg'][:, 0]
    data['time'] = dataset.variables['time'][:]
    data['lon'] = dataset.variables['lon'][:]
    data['lat'] = dataset.variables['lat'][:]
    data['pre_lev'] = dataset.variables['musica_pressure_levels'][:]
    data['h2o_lev_dry'] = dataset.variables['musica_wv'][:, 0]
    data['tem_lev'] = dataset.variables['musica_at'][:]
    data['quality_flag'] = dataset.variables['musica_fit_quality_flag'][:]
    data['observation_id'] = dataset.variables['observation_id'][:]
    data['alt_lev'] = dataset.variables['musica_altitude_levels'][:].filled(np.nan)
    data['num_lev'] = dataset.variables['musica_nol'][:]
    data['apri'] = dataset.variables['musica_ghg_apriori'][:, 0]

    # data['avk_rank'] = dataset.variables['musica_ghg_avk_rank'][:]
    #
    # # make that shit stackable
    # if dataset.variables['musica_ghg_avk_val'][:, :].shape[1] != 16:
    #     profile_count = dataset.variables['musica_ghg_avk_lvec'][:, 0, :, :].shape[0]
    #     nan_count = 16 - dataset.variables['musica_ghg_avk_lvec'][:, 0, :, :].shape[1]
    #     nol_count = dataset.variables['musica_ghg_avk_lvec'][:, 0, :, :].shape[2]
    #
    #     nan_slice_val = np.ma.masked_array(np.full((profile_count, nan_count), np.nan))
    #     nan_slice_vec = np.ma.masked_array(np.full((profile_count, nan_count, nol_count), np.nan))
    #
    #     data['avk_val'] = np.concatenate(
    #         [dataset.variables['musica_ghg_avk_val'][:, :], nan_slice_val], axis=1)
    #
    #     data['avk_lvec'] = np.concatenate(
    #         [dataset.variables['musica_ghg_avk_lvec'][:, 0, :, :], nan_slice_vec], axis=1)
    #
    #     data['avk_rvec'] = np.concatenate(
    #         [dataset.variables['musica_ghg_avk_rvec'][:, 0, :, :], nan_slice_vec], axis=1)
    #
    # else:
    #     data['avk_val'] = dataset.variables['musica_ghg_avk_val'][:, :]
    #     data['avk_lvec'] = dataset.variables['musica_ghg_avk_lvec'][:, 0, :, :]
    #     data['avk_rvec'] = dataset.variables['musica_ghg_avk_rvec'][:, 0, :, :]

    dataset.close()

    return data


def read_iasi_day(dir, quality_flag):
    date = dir.split('/')[-2:]

    paths = glob2.glob(join(dir, '*.nc'))
    paths.sort()

    combined_dict = read_iasi(paths[0])

    for path in paths[1:]:
        data = read_iasi(path)
        for key in combined_dict:
            combined_dict[key] = np.concatenate((combined_dict[key], data[key]), axis=0)

    combined_dict['readable_time'] = conv_iasi_time_readable(combined_dict['time'])
    low_mask = combined_dict['readable_time'] >= datetime(2020, int(date[0]), int(date[1]), 0, 0, 0)
    upp_mask = combined_dict['readable_time'] <= datetime(2020, int(date[0]), int(date[1]), 23, 59, 59)
    mask = low_mask & upp_mask

    for key in combined_dict.keys():
        combined_dict[key] = combined_dict[key][mask]

    # filtering quality_flag
    quality_mask = combined_dict['quality_flag'] >= quality_flag
    for key in combined_dict.keys():
        combined_dict[key] = combined_dict[key][quality_mask]

    flip_keys = ['n2o_lev_dry', 'pre_lev', 'h2o_lev_dry', 'tem_lev', 'alt_lev', 'apri']
    for key in flip_keys:
        combined_dict[key] = combined_dict[key][:, ::-1]

    return combined_dict


def iasi_day_profiles(dir, quality_flag, coord):
    iasi_data = read_iasi_day(dir=dir, quality_flag=quality_flag)

    mask = (
            (iasi_data['lon'] >= coord[0][0]) & (iasi_data['lon'] <= coord[0][1]) &
            (iasi_data['lat'] >= coord[1][0]) & (iasi_data['lat'] <= coord[1][1])
    )

    # Apply mask
    iasi_profiles = iasi_data['apri'][mask]
    lon = iasi_data['lon'][mask]
    lat = iasi_data['lon'][mask]
    pre_lev = iasi_data['pre_lev'][mask]
    h2o_lev = iasi_data['h2o_lev_dry'][mask]
    alt_lev = iasi_data['alt_lev'][mask]

    nan_counts = np.isnan(iasi_profiles).sum(axis=1)

    gosat_profiles = np.full(iasi_profiles.shape, np.nan)

    items = list(range(iasi_profiles.shape[0]))
    for row in tqdm(items):

        alt_lay = lev2lay(alt_lev[row, nan_counts[row]:])

        gra_geo = grav_geoid(lat=lat[row])

        gra_lay = grav_acc(gra_geo=gra_geo, gmh_lay=alt_lay)

        h2o_lay = lev2lay(h2o_lev[row, nan_counts[row]:])

        dry_col = dry_column(pre_lev=pre_lev[row, nan_counts[row]:], h2o_dry=h2o_lay, gra_lay=gra_lay)

        xn2o = xn2o_norm  # ppm
        n2o_prof = gosat_prior

        n2o_prof = xn2o * n2o_prof / np.mean(n2o_prof)  # GOSAT prior

        pre_lay = layer_mid_pressure(pre_lev[row, nan_counts[row]:])

        cum_sum_par = np.cumsum(dry_col)  # Particle count for interpolation

        # Particles per layer for GOSAT minus the last layer bc this would give the last iasi level by design
        n_lay = [i * np.sum(dry_col) / len(n2o_prof) for i in range(1, len(n2o_prof))]

        gosat_pre_lev = np.zeros(len(n2o_prof) + 1)  # array for gosat pre levels
        gosat_pre_lev[0] = pre_lev[row, nan_counts[row]:][0]
        gosat_pre_lev[-1] = 0

        # interpolation
        gosat_pre_lev[1:-1] = np.interp(n_lay, cum_sum_par, pre_lay)

        # cal weights for weighted average
        deltas = np.zeros((len(gosat_pre_lev[:-1]), len(pre_lay)))  # -1 element for layer instead of level

        for i in range(deltas.shape[0]):
            gpre_1 = gosat_pre_lev[i]
            gpre_2 = gosat_pre_lev[i + 1]
            for j in range(deltas.shape[1]):
                ipre_1 = pre_lev[row, nan_counts[row]:][j]
                ipre_2 = pre_lev[row, nan_counts[row]:][j + 1]

                if ipre_2 >= gpre_1:  # iasi layer is under gosat layer
                    deltas[i, j] = 0.
                elif (ipre_1 <= gpre_1) and (ipre_2 >= gpre_2):  # iasi layer completely  inside gosat layer
                    deltas[i, j] = 1.
                elif (ipre_1 > gpre_1) and (ipre_2 < gpre_1):  # cut from below
                    i_del = ipre_1 - ipre_2
                    pre_del = gpre_1 - ipre_2
                    if ipre_2 < gpre_2:  # edge case if 1 iasi layer touches 3 gosat layers
                        overshoot = gpre_2 - ipre_2
                    else:
                        overshoot = 0.
                    deltas[i, j] = (pre_del - overshoot) / i_del
                elif (ipre_1 > gpre_2) and (ipre_2 < gpre_2):  # cut from top
                    i_del = ipre_1 - ipre_2
                    pre_del = ipre_1 - gpre_2
                    if ipre_1 > gpre_1:  # edge case if 1 iasi layer touches 3 gosat layers
                        overshoot = ipre_1 - gpre_1
                    else:
                        overshoot = 0.
                    deltas[i, j] = (pre_del - overshoot) / i_del
                elif ipre_1 <= gpre_2:  # iasi layer over gosat layer
                    deltas[i, j] = 0.

        # gosat n2o on iasi layers
        n2o_gosat_lay = np.sum(n2o_prof[:, np.newaxis] * deltas, axis=0)  # calc weighted average

        gosat_n2o_lev = np.full(pre_lev[row, nan_counts[row]:].shape, np.nan)

        # interpolate n2o ppm values for iasi lev. All flipped bc np.interp needs ascending values
        gosat_n2o_lev[1:-1] = np.interp(pre_lev[row, nan_counts[row]:][1:-1][::-1], pre_lay[::-1], n2o_gosat_lay[::-1])[
                              ::-1]

        # calc first and last level n2o ppm value
        gosat_n2o_lev[0] = 2 * n2o_gosat_lay[0] - gosat_n2o_lev[1]
        gosat_n2o_lev[-1] = 2 * n2o_gosat_lay[-1] - gosat_n2o_lev[-2]

        gosat_profiles[row, nan_counts[row]:] = gosat_n2o_lev

    return iasi_profiles, gosat_profiles, pre_lev


def inter_profiles(profiles, pre_lev):
    p_min = np.nanmin(pre_lev[:, :])
    p_max = np.nanmax(pre_lev[:, :])

    # Define log-spaced pressure grid (e.g. 40 levels)
    pressure_grid = np.logspace(np.log10(p_max), np.log10(p_min), 40)

    interp_profiles = []

    print('Interpolating:')
    for i in tqdm(range(profiles.shape[0])):
        p_profile = pre_lev[i, :]
        n2o_profile = profiles[i, :]

        mask = ~np.isnan(p_profile) & ~np.isnan(n2o_profile)
        if np.sum(mask) < 2:  # need at least 2 points to interpolate
            continue

        interp_n2o = np.interp(pressure_grid, p_profile[mask][::-1], n2o_profile[mask][::-1])
        interp_profiles.append(interp_n2o)

    interp_profiles = np.array(interp_profiles)
    mean_profile = np.nanmean(interp_profiles, axis=0)

    return mean_profile * 1000, pressure_grid / 100



