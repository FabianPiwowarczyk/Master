import netCDF4 as nc
import pandas as pd
import numpy as np
import glob2
from os.path import join
import xarray as xr


def grid_centers(x_res, y_res):
    lat_bnd = (-90, 90)
    lon_bnd = (-180, 180)

    lat = np.array([lat_bnd[0] + i + y_res/2 for i in range(0, lat_bnd[1] + abs(lat_bnd[0]), y_res)])
    lon = np.array([lon_bnd[0] + i + x_res / 2 for i in range(0, lon_bnd[1] + abs(lon_bnd[0]), x_res)])

    return lon, lat


def find_cen(value, centers, res):
    for center in centers:
        if center - (res / 2) <= value < center + (res / 2):
            return center


def read_df(dir_path, m, tot_var, lon_var, lat_var, sat, qf=None):

    if sat == 'iasi':
        paths = glob2.glob(dir_path + m + f'*qf{qf}.nc')
    else:
        paths = glob2.glob(dir_path + m + '*.nc')
    paths.sort()

    tot_col_list = []
    lon_list = []
    lat_list = []

    # reading and combining all data files
    for path in paths:
        dataset = nc.Dataset(path, mode='r')

        tot_col_list.append(dataset.variables[tot_var][:])
        lon_list.append(dataset.variables[lon_var][:])
        lat_list.append(dataset.variables[lat_var][:])

        dataset.close()

    tot_col_combined = np.concatenate(tot_col_list)
    lon_combined = np.concatenate(lon_list)
    lat_combined = np.concatenate(lat_list)

    if sat == 'iasi':
        tot_col_combined = tot_col_combined * 1000

    df = pd.DataFrame({
        'tot_col': tot_col_combined,
        'lon': lon_combined,
        'lat': lat_combined
    })

    return df


def grid(df, lon_cen, x_res, lat_cen, y_res, th):
    # assigning grid box centers
    df['lon_cen'] = df['lon'].apply(lambda x: find_cen(x, lon_cen, x_res))
    df['lat_cen'] = df['lat'].apply(lambda x: find_cen(x, lat_cen, y_res))

    # averaging
    df_grouped = df.groupby(['lon_cen', 'lat_cen'], as_index=False).agg(
        mean_tot=('tot_col', 'mean'),
        std_tot=('tot_col', 'std'),
        count=('tot_col', 'count'))

    # filter low count boxes
    df_filtered = df_grouped[df_grouped['count'] >= th]

    # filling every box without value with np.nan
    full_grid = pd.MultiIndex.from_product([lon_cen, lat_cen], names=['lon_cen', 'lat_cen'])
    full_grid_df = pd.DataFrame(index=full_grid).reset_index()

    df_filled = pd.merge(full_grid_df, df_filtered, on=['lon_cen', 'lat_cen'], how='left')

    return df_filled


def monthly_mean(x_res, y_res, th, dir_path, tot_var, lon_var, lat_var, sat, qf=None):

    months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

    lon_cen, lat_cen = grid_centers(x_res=x_res, y_res=y_res)

    for m in months:
        if sat == 'iasi':
            df = read_df(dir_path, m, tot_var, lon_var, lat_var, sat, qf)
        else:
            df = read_df(dir_path, m, tot_var, lon_var, lat_var, sat)
        df_filled = grid(df, lon_cen, x_res, lat_cen, y_res, th)

        # saving the data
        pivot_mean_tot = df_filled.pivot(index='lat_cen', columns='lon_cen', values='mean_tot')
        mean_tot = pivot_mean_tot.values

        pivot_std_tot = df_filled.pivot(index='lat_cen', columns='lon_cen', values='std_tot')
        std_tot = pivot_std_tot.values

        pivot_count = df_filled.pivot(index='lat_cen', columns='lon_cen', values='count')
        count = pivot_count.values

        lon_centers = df_filled['lon_cen'].unique()
        lat_centers = df_filled['lat_cen'].unique()

        ds = xr.Dataset(
            {
                'mean_tot': (['lat_cen', 'lon_cen'], mean_tot),
                'std_tot': (['lat_cen', 'lon_cen'], std_tot),
                'count': (['lat_cen', 'lon_cen'], count)
                },
                coords={
                    'lon_cen': lon_centers,
                    'lat_cen': lat_centers
                }
            )

        print(f'Saving {m}')
        if sat == 'iasi':
            if tot_var in ['tc_cor_met0', 'tc_cor_met1', 'tc_cor_met2']:
                output_path = f'monthly_means/{sat}_met{tot_var[-1]}_{m}_{x_res}x{y_res}_th{th}_qf{qf}.nc'
            else:
                output_path = f'monthly_means/{sat}_{m}_{x_res}x{y_res}_th{th}_qf{qf}.nc'
        else:
            output_path = f'monthly_means/{sat}_{m}_{x_res}x{y_res}_th{th}.nc'
        ds.to_netcdf(output_path)
