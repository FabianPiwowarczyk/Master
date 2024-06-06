import glob2
import netCDF4 as nc
import numpy as np

from os.path import join


def read_iasi(path):
    data = dict()

    dataset = nc.Dataset(path, mode='r')

    data['n2o'] = dataset.variables['musica_ghg'][:, 0]
    data['time'] = dataset.variables['time'][:]
    data['lon'] = dataset.variables['lon'][:]
    data['lat'] = dataset.variables['lat'][:]
    data['pressure_levels'] = dataset.variables['musica_pressure_levels'][:]
    data['h2o'] = dataset.variables['musica_wv'][:, 0]

    dataset.close()

    return data


def _file_paths(dir):
    paths = glob2.glob(join(f'{dir}', '*.nc'))
    return paths


def read_all_iasi(dir):
    paths = _file_paths(dir)
    combined_dict = read_iasi(paths[0])

    for path in paths[1:]:
        data = read_iasi(path)
        for key in combined_dict:
            combined_dict[key] = np.concatenate((combined_dict[key], data[key]), axis=0)

    return combined_dict
