import xarray as xr

from gridding.grid_func import grid, grid_centers
from .cams_tc import cams_tc


def cams_monthly_means(x_res, y_res, months):

    months[1] = months[1] + 1
    lon_cen, lat_cen = grid_centers(x_res=x_res, y_res=y_res)
    for m in range(*months):

        print(f'Griding Month {m}')
        df_cams, df_iasi = cams_tc(m)
        df_cams_filled = grid(df_cams, lon_cen, x_res, lat_cen, y_res, th=0)
        df_iasi_filled = grid(df_iasi, lon_cen, x_res, lat_cen, y_res, th=0)

        save_data(df_cams_filled, x_res, y_res, m, 'cams')
        save_data(df_iasi_filled, x_res, y_res, m, 'cams_iasi')


def save_data(df_filled, x_res, y_res, m, path):
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
    output_path = f'monthly_means/{path}_{m:02d}_{x_res}x{y_res}_th{0}.nc'
    ds.to_netcdf(output_path)
