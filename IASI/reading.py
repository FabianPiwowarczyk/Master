import glob2
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta

from os.path import join


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

    data['avk_rank'] = dataset.variables['musica_ghg_avk_rank'][:]

    # make that shit stackable
    if dataset.variables['musica_ghg_avk_val'][:, :].shape[1] != 16:
        profile_count = dataset.variables['musica_ghg_avk_lvec'][:, 0, :, :].shape[0]
        nan_count = 16 - dataset.variables['musica_ghg_avk_lvec'][:, 0, :, :].shape[1]
        nol_count = dataset.variables['musica_ghg_avk_lvec'][:, 0, :, :].shape[2]

        nan_slice_val = np.ma.masked_array(np.full((profile_count, nan_count), np.nan))
        nan_slice_vec = np.ma.masked_array(np.full((profile_count, nan_count, nol_count), np.nan))

        data['avk_val'] = np.concatenate(
            [dataset.variables['musica_ghg_avk_val'][:, :], nan_slice_val], axis=1)

        data['avk_lvec'] = np.concatenate(
            [dataset.variables['musica_ghg_avk_lvec'][:, 0, :, :], nan_slice_vec], axis=1)

        data['avk_rvec'] = np.concatenate(
            [dataset.variables['musica_ghg_avk_rvec'][:, 0, :, :], nan_slice_vec], axis=1)

    else:
        data['avk_val'] = dataset.variables['musica_ghg_avk_val'][:, :]
        data['avk_lvec'] = dataset.variables['musica_ghg_avk_lvec'][:, 0, :, :]
        data['avk_rvec'] = dataset.variables['musica_ghg_avk_rvec'][:, 0, :, :]

    dataset.close()

    return data


def _file_paths(dir, i, date_tuples, org_path):
    if i != 0:
        paths_day_before = glob2.glob(join(org_path, date_tuples[i-1, 0], date_tuples[i-1, 1], '*.nc'))
        paths_day_before.sort()

        paths = glob2.glob(join(f'{dir}', '*.nc'))
        paths.sort()
        paths.insert(0, paths_day_before[-1])
    else:
        paths = glob2.glob(join(f'{dir}', '*.nc'))
        paths.sort()
    return paths


def read_all_iasi(dir, date, i, date_tuples, org_path, quality_flag):
    paths = _file_paths(dir, i, date_tuples, org_path)
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
