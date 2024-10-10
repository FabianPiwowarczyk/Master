import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np


def bnds(cen):
    cell_bnds = np.zeros(cen.size + 1)
    cell_bnds[1:-1] = (cen[:-1] + cen[1:]) / 2
    cell_bnds[0] = cen[0] - (cen[1] - cen[0]) / 2
    cell_bnds[-1] = cen[-1] + (cen[-1] - cen[-2]) / 2
    return cell_bnds


def plot_lv3_data(x_res, y_res, th, dir_path, sat, vmin=None, vmax=None):

    months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

    for m in months:
        ds = xr.open_dataset(dir_path.format(m, x_res, y_res, th))

        lon = ds['lon_cen'].values
        lat = ds['lat_cen'].values
        tot_col = ds['mean_tot'].values
        count = ds['count'].values
        std = ds['std_tot'].values

        lon_bnds = bnds(lon)
        lat_bnds = bnds(lat)

        plt.figure(figsize=(12, 6))
        ax = plt.axes(projection=ccrs.PlateCarree())

        ax.coastlines()
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.LAND, facecolor='white')
        ax.add_feature(cfeature.OCEAN, facecolor='white')

        # Add gridlines and set latitude/longitude labels
        gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray')
        gl.top_labels = gl.right_labels = False
        gl.xlabel_style = {'size': 10, 'color': 'black'}
        gl.ylabel_style = {'size': 10, 'color': 'black'}

        plt.pcolormesh(lon_bnds, lat_bnds, tot_col, transform=ccrs.PlateCarree(), cmap='jet'
                       , vmin=vmin, vmax=vmax)

        plt.colorbar(label=r'$\mathrm{xN}_2\mathrm{O}$ in ppm')
        plt.title(f'{sat} {x_res}x{y_res} gridded data 2020-{m}')

        outpath = f'pictures/{sat}_{x_res}x{y_res}_{m}.png'
        plt.savefig(outpath)


def combined_plot(x_res, y_res, th, dir_path, var, vmin=None, vmax=None):

    months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

    for m in months:
        ds_iasi = xr.open_dataset(dir_path.format('iasi', m, x_res, y_res, th))
        ds_gosat = xr.open_dataset(dir_path.format('gosat', m, x_res, y_res, th))

        lon_iasi = ds_iasi['lon_cen'].values
        lat_iasi = ds_iasi['lat_cen'].values
        variable_iasi = ds_iasi[var].values

        variable_gosat = ds_gosat[var].values

        lon_bnds = bnds(lon_iasi)
        lat_bnds = bnds(lat_iasi)

        # Create a figure with two subplots side by side
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6), subplot_kw={'projection': ccrs.PlateCarree()})

        # ---------------- Plot 1: iasi ----------------
        ax1.coastlines()
        ax1.add_feature(cfeature.BORDERS)
        ax1.add_feature(cfeature.LAND, facecolor='white')
        ax1.add_feature(cfeature.OCEAN, facecolor='white')

        # Add gridlines and set latitude/longitude labels
        gl1 = ax1.gridlines(draw_labels=True, linestyle='--', color='gray')
        gl1.top_labels = gl1.right_labels = False
        gl1.xlabel_style = {'size': 10, 'color': 'black'}
        gl1.ylabel_style = {'size': 10, 'color': 'black'}

        # Plot the tot_col variable for IASI
        iasi_plot = ax1.pcolormesh(lon_bnds, lat_bnds, variable_iasi, transform=ccrs.PlateCarree(), cmap='jet',
                                   vmin=vmin, vmax=vmax)
        ax1.set_title('IASI')

        # ---------------- Plot 2: gosat ----------------
        ax2.coastlines()
        ax2.add_feature(cfeature.BORDERS)
        ax2.add_feature(cfeature.LAND, facecolor='white')
        ax2.add_feature(cfeature.OCEAN, facecolor='white')

        # Add gridlines and set latitude/longitude labels
        gl2 = ax2.gridlines(draw_labels=True, linestyle='--', color='gray')
        gl2.top_labels = gl2.right_labels = False
        gl2.xlabel_style = {'size': 10, 'color': 'black'}
        gl2.ylabel_style = {'size': 10, 'color': 'black'}

        # Plot the count variable for GOSAT
        gosat_plot = ax2.pcolormesh(lon_bnds, lat_bnds, variable_gosat, transform=ccrs.PlateCarree(), cmap='jet',
                                    vmin=vmin, vmax=vmax)
        ax2.set_title('GOSAT')

        # ---------------- Shared Colorbar ----------------
        # Create a shared colorbar for both subplots
        cbar = fig.colorbar(iasi_plot, ax=[ax1, ax2], orientation='vertical', fraction=0.02, pad=0.04)
        cbar.set_label(r'$\mathrm{xN}_2\mathrm{O}$ in ppb')

        outpath = f'pictures/combined_{x_res}x{y_res}_{m}.png'
        plt.savefig(outpath)


def plot_dif(x_res, y_res, th, dir_path, vmin=None, vmax=None):

    months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

    # mini = []
    # maxi = []

    for m in months:
        ds_iasi = xr.open_dataset(dir_path.format('iasi', m, x_res, y_res, th))
        ds_gosat = xr.open_dataset(dir_path.format('gosat', m, x_res, y_res, th))

        lon = ds_iasi['lon_cen'].values
        lat = ds_iasi['lat_cen'].values
        tot_col_iasi = ds_iasi['mean_tot'].values
        tot_col_gosat = ds_gosat['mean_tot'].values

        com_col = tot_col_gosat - tot_col_iasi

        # mini.append(np.nanmin(com_col))
        # maxi.append(np.nanmax(com_col))

        lon_bnds = bnds(lon)
        lat_bnds = bnds(lat)

        plt.figure(figsize=(12, 6))
        ax = plt.axes(projection=ccrs.PlateCarree())

        ax.coastlines()
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.LAND, facecolor='white')
        ax.add_feature(cfeature.OCEAN, facecolor='white')

        # Add gridlines and set latitude/longitude labels
        gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray')
        gl.top_labels = gl.right_labels = False
        gl.xlabel_style = {'size': 10, 'color': 'black'}
        gl.ylabel_style = {'size': 10, 'color': 'black'}

        plt.pcolormesh(lon_bnds, lat_bnds, com_col, transform=ccrs.PlateCarree(), cmap='jet'
                       , vmin=vmin, vmax=vmax)

        plt.colorbar(label=r'$\mathrm{xN}_2\mathrm{O}$ in ppm')
        plt.title(f'gosat minus iasi {x_res}x{y_res} gridded data 2020-{m}')

        outpath = f'pictures/dif_{x_res}x{y_res}_{m}.png'
        plt.savefig(outpath)

    # print(np.mean(mini), np.mean(maxi))