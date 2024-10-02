import netCDF4 as nc
import pandas as pd
import numpy as np


def monthly_mean(x_res, y_res):

    path = 'finished_data/2020_01_01.nc'

    dataset = nc.Dataset(path, mode='r')

    # for var_name in dataset.variables:
    #     print(var_name)
    #     print(dataset.variables[var_name].description)

    df = pd.DataFrame({
        'tot_col': dataset.variables['total_column'][:],
        'lon': dataset.variables['lon'][:],
        'lat': dataset.variables['lat'][:]
    })
    dataset.close()

    lon_cen, lat_cen = grid_centers(x_res=x_res, y_res=y_res)

    def find_cen(value, centers, res):
        for center in centers:
            if center - (res / 2) <= value < center + (res / 2):
                return center

    # Apply the function to each longitude value in the DataFrame
    df['lon_cen'] = df['lon'].apply(lambda x: find_cen(x, lon_cen, x_res))
    df['lat_cen'] = df['lat'].apply(lambda x: find_cen(x, lat_cen, y_res))

    df_grouped = df.groupby(['lon_cen', 'lat_cen'], as_index=False).agg(
    mean_tot=('tot_col', 'mean'),
    std_tot=('tot_col', 'std'),
    count=('tot_col', 'count')
)

    print(df_grouped.head(10))


def grid_centers(x_res, y_res):
    lat_bnd = (-90, 90)
    lon_bnd = (-180, 180)

    lat = np.array([lat_bnd[0] + i + y_res/2 for i in range(0, lat_bnd[1] + abs(lat_bnd[0]), y_res)])
    lon = np.array([lon_bnd[0] + i + x_res / 2 for i in range(0, lon_bnd[1] + abs(lon_bnd[0]), x_res)])

    return lon, lat
