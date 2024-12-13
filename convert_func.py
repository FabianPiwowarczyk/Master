import netCDF4 as nc
import numpy as np
from os.path import join
import os

from IASI.total_column_conversion import data2total_col


def convert():
    data_path = '/misc/hypatia/data/IASI/L2_MUSICA/V3.30/MetopA/2020'

    quality_flag = 3

    months_paths = [f.path for f in os.scandir(data_path) if f.is_dir()]
    months_paths.sort()
    months = [m.replace(data_path, '').replace('/', '') for m in months_paths]

    date_array = np.array([0, 0])
    for m in months:
        days_paths = [f.path for f in os.scandir(join(data_path, m)) if f.is_dir()]
        days_paths.sort()
        days = [d.replace(join(data_path, m), '').replace('/', '') for d in days_paths]
        for day in days:
            date_array = np.vstack((date_array, np.array([m, day])))

    date_tuples = date_array[1:, :]

    for i, date in enumerate(date_tuples):
        out_path = f'new_finished_data/2020_{date[0]}_{date[1]}.nc'
        if os.path.exists(out_path):
            print('File: ', out_path, ' already exists.')
        else:
            path = join(data_path, str(date[0]), str(date[1]))
            print(path)

            data = data2total_col(path=path, date=date, i=i, date_tuples=date_tuples,
                                  org_path=data_path, quality_flag=quality_flag)

            descriptions = {
                'n2o_lev_dry': "N2O leveled data, dry air mole fraction (ppmv).",
                'time': "UTC time in seconds after 2000-1-1 00:00:00.",
                'lon': "longitude.",
                'lat': "latitude.",
                'pre_lev': "pressure levels.",
                'h2o_lev_dry': "water vapor, dry air mole fraction (ppmv).",
                'tem_lev': "atmospheric temperatures.",
                'total_column': "total column in ppm.",
                'alt_lev': "Altitude levels in meters above sea level."
            }

            for key in [el for el in data.keys()]:
                if key not in descriptions.keys():
                    del data[key]

            print('Creating output file:', out_path)
            ds = nc.Dataset(out_path, 'w', format='NETCDF4')

            num_entries = data['n2o_lev_dry'].shape[0]
            ds.createDimension('index', num_entries)
            level_entries = data['n2o_lev_dry'].shape[1]
            ds.createDimension('level', level_entries)

            for key, values in data.items():

                twodim_keys = ['n2o_lev_dry', 'pre_lev', 'h2o_lev_dry', 'tem_lev', 'alt_lev']

                if key in twodim_keys:
                    var = ds.createVariable(key, values.dtype, ('index', 'level'))
                    var[:, :] = values
                    var.description = descriptions[key]
                else:
                    var = ds.createVariable(key, values.dtype, ('index',))
                    var[:] = values
                    var.description = descriptions[key]

            ds.close()