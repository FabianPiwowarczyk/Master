import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib.colors as mcolors


def bnds(cen):
    cell_bnds = np.zeros(cen.size + 1)
    cell_bnds[1:-1] = (cen[:-1] + cen[1:]) / 2
    cell_bnds[0] = cen[0] - (cen[1] - cen[0]) / 2
    cell_bnds[-1] = cen[-1] + (cen[-1] - cen[-2]) / 2
    return cell_bnds


def plot_lv3_data(x_res, y_res, th, dir_path, sat, var='mean_tot', vmin=None, vmax=None, qf=None, met=None):

    months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

    for m in months:
        if sat == 'iasi':
            if met in [0, 1, 2]:
                ds = xr.open_dataset(dir_path.format(sat, met, m, x_res, y_res, th, qf))
            else:
                ds = xr.open_dataset(dir_path.format(sat, m, x_res, y_res, th, qf))
        else:
            ds = xr.open_dataset(dir_path.format(sat, m, x_res, y_res, th))

        lon = ds['lon_cen'].values
        lat = ds['lat_cen'].values
        plot_var = ds[var].values

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

        plt.pcolormesh(lon_bnds, lat_bnds, plot_var, transform=ccrs.PlateCarree(), cmap='jet'
                       , vmin=vmin, vmax=vmax)

        plt.colorbar(label=r'$\mathrm{xN}_2\mathrm{O}$ in ppb')
        if met in [0, 1, 2]:
            plt.title(f'{sat} {x_res}x{y_res} {var} met{met} 2020-{m}')
        else:
            plt.title(f'{sat} {x_res}x{y_res} {var} 2020-{m}')

        if sat == 'iasi':
            if met in [0, 1, 2]:
                outpath = f'pictures/{sat}_{var}_met{met}_{x_res}x{y_res}_qf{qf}_{m}.png'
            else:
                outpath = f'pictures/{sat}_{var}_{x_res}x{y_res}_qf{qf}_{m}.png'
        else:
            outpath = f'pictures/{sat}_{var}_{x_res}x{y_res}_{m}.png'
        plt.savefig(outpath)


def combined_plot(x_res, y_res, th, dir_path, var, qf, vmin=None, vmax=None,
                  met_path='monthly_means/{}_met{}_{}_{}x{}_th{}_qf{}.nc', met=0):

    months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)  # norm for colorbar

    for m in months:
        ds_iasi = xr.open_dataset(dir_path.format('iasi', m, x_res, y_res, th, qf))
        ds_met = xr.open_dataset(met_path.format('iasi', met, m, x_res, y_res, th, qf))

        gosat_path = 'monthly_means/{}_{}_{}x{}_th{}.nc'
        ds_gosat = xr.open_dataset(gosat_path.format('gosat', m, x_res, y_res, th))

        lon_iasi = ds_iasi['lon_cen'].values
        lat_iasi = ds_iasi['lat_cen'].values
        variable_iasi = ds_iasi[var].values
        variable_met = ds_met[var].values
        variable_gosat = ds_gosat[var].values

        var1 = variable_iasi - variable_gosat
        var2 = variable_met - variable_gosat

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

        # Plot the tot_col variable for IASI - GOSAT
        iasi_plot = ax1.pcolormesh(lon_bnds, lat_bnds, var1, transform=ccrs.PlateCarree(), cmap='seismic',
                                   vmin=vmin, vmax=vmax, norm=norm)
        ax1.set_title(f'IASI - GOSAT {x_res}x{y_res} grid 2020-{m}')

        # ---------------- Plot 2: methode --------------
        ax2.coastlines()
        ax2.add_feature(cfeature.BORDERS)
        ax2.add_feature(cfeature.LAND, facecolor='white')
        ax2.add_feature(cfeature.OCEAN, facecolor='white')

        # Add gridlines and set latitude/longitude labels
        gl2 = ax2.gridlines(draw_labels=True, linestyle='--', color='gray')
        gl2.top_labels = gl2.right_labels = False
        gl2.xlabel_style = {'size': 10, 'color': 'black'}
        gl2.ylabel_style = {'size': 10, 'color': 'black'}

        # Plot the count variable for IASI MET0 - GOSAT
        gosat_plot = ax2.pcolormesh(lon_bnds, lat_bnds, var2, transform=ccrs.PlateCarree(), cmap='seismic',
                                    vmin=vmin, vmax=vmax, norm=norm)
        ax2.set_title(f'Cor met{met} IASI - GOSAT {x_res}x{y_res} grid 2020-{m}')

        # ---------------- Shared Colorbar ----------------
        # Create a shared colorbar for both subplots
        cbar = fig.colorbar(iasi_plot, ax=[ax1, ax2], orientation='vertical', fraction=0.02, pad=0.04)
        cbar.set_label(r'$\mathrm{xN}_2\mathrm{O}$ in ppb')

        outpath = f'pictures/combined_met{met}_{x_res}x{y_res}_{m}.png'
        plt.savefig(outpath)


def iasi_dif_plot(x_res, y_res, th, dir_path, qf, vmin=None, vmax=None,
                  met_path='monthly_means/{}_met{}_{}_{}x{}_th{}_qf{}.nc'):

    months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

    for m in months:
        ds_iasi = xr.open_dataset(dir_path.format('iasi', m, x_res, y_res, th, qf))
        ds_met0 = xr.open_dataset(met_path.format('iasi', '0', m, x_res, y_res, th, qf))
        ds_met1 = xr.open_dataset(met_path.format('iasi', '1', m, x_res, y_res, th, qf))
        ds_met2 = xr.open_dataset(met_path.format('iasi', '2', m, x_res, y_res, th, qf))

        lon_iasi = ds_iasi['lon_cen'].values
        lat_iasi = ds_iasi['lat_cen'].values
        variable_iasi = ds_iasi['mean_tot'].values
        variable_met0 = ds_met0['mean_tot'].values
        variable_met1 = ds_met1['mean_tot'].values
        variable_met2 = ds_met2['mean_tot'].values

        var0 = variable_iasi - variable_met0
        var1 = variable_iasi - variable_met1
        var2 = variable_iasi - variable_met2

        lon_bnds = bnds(lon_iasi)
        lat_bnds = bnds(lat_iasi)

        # Create a figure with two subplots side by side
        fig, axs = plt.subplots(2, 2, figsize=(18, 12), subplot_kw={'projection': ccrs.PlateCarree()})

        # ---------------- Plot 0: Met 0 ----------------
        axs[0, 0].coastlines()
        axs[0, 0].add_feature(cfeature.BORDERS)
        axs[0, 0].add_feature(cfeature.LAND, facecolor='white')
        axs[0, 0].add_feature(cfeature.OCEAN, facecolor='white')

        # Add gridlines and set latitude/longitude labels
        gl0 = axs[0, 0].gridlines(draw_labels=True, linestyle='--', color='gray')
        gl0.top_labels = gl0.right_labels = False
        gl0.xlabel_style = {'size': 10, 'color': 'black'}
        gl0.ylabel_style = {'size': 10, 'color': 'black'}

        # Plot
        plot0 = axs[0, 0].pcolormesh(lon_bnds, lat_bnds, var0, transform=ccrs.PlateCarree(), cmap='jet',
                                   vmin=vmin, vmax=vmax)
        axs[0, 0].set_title(f'IASI - MET 0 {x_res}x{y_res} grid 2020-{m}')

        # ---------------- Plot 1: Met 1 ----------------
        axs[0, 1].coastlines()
        axs[0, 1].add_feature(cfeature.BORDERS)
        axs[0, 1].add_feature(cfeature.LAND, facecolor='white')
        axs[0, 1].add_feature(cfeature.OCEAN, facecolor='white')

        # Add gridlines and set latitude/longitude labels
        gl1 = axs[0, 1].gridlines(draw_labels=True, linestyle='--', color='gray')
        gl1.top_labels = gl1.right_labels = False
        gl1.xlabel_style = {'size': 10, 'color': 'black'}
        gl1.ylabel_style = {'size': 10, 'color': 'black'}

        # Plot
        plot1 = axs[0, 1].pcolormesh(lon_bnds, lat_bnds, var1, transform=ccrs.PlateCarree(), cmap='jet',
                                   vmin=vmin, vmax=vmax)
        axs[0, 1].set_title(f'IASI - MET 1 {x_res}x{y_res} grid 2020-{m}')

        # ---------------- Plot 2: Met 2 ----------------
        axs[1, 0].coastlines()
        axs[1, 0].add_feature(cfeature.BORDERS)
        axs[1, 0].add_feature(cfeature.LAND, facecolor='white')
        axs[1, 0].add_feature(cfeature.OCEAN, facecolor='white')

        # Add gridlines and set latitude/longitude labels
        gl2 = axs[1, 0].gridlines(draw_labels=True, linestyle='--', color='gray')
        gl2.top_labels = gl2.right_labels = False
        gl2.xlabel_style = {'size': 10, 'color': 'black'}
        gl2.ylabel_style = {'size': 10, 'color': 'black'}

        # Plot
        plot2 = axs[1, 0].pcolormesh(lon_bnds, lat_bnds, var2, transform=ccrs.PlateCarree(), cmap='jet',
                                   vmin=vmin, vmax=vmax)
        axs[1, 0].set_title(f'IASI - MET 2 {x_res}x{y_res} grid 2020-{m}')

        # ---------------- Shared Colorbar ----------------
        # Create a shared colorbar for both subplots
        cbar = fig.colorbar(plot0, ax=axs, orientation='vertical', fraction=0.02, pad=0.04)
        cbar.set_label(r'$\mathrm{xN}_2\mathrm{O}$ in ppb')

        outpath = f'pictures/iasi_dif_{x_res}x{y_res}_{m}.png'
        plt.savefig(outpath)


def plot_dif(x_res, y_res, th, dir_path, vmin=None, vmax=None):
    """Not yet updated for newer data format"""

    months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

    mini = []
    maxi = []

    for m in months:
        ds_iasi = xr.open_dataset(dir_path.format('iasi', m, x_res, y_res, th))
        ds_gosat = xr.open_dataset(dir_path.format('gosat', m, x_res, y_res, th))

        lon = ds_iasi['lon_cen'].values
        lat = ds_iasi['lat_cen'].values
        tot_col_iasi = ds_iasi['mean_tot'].values
        tot_col_gosat = ds_gosat['mean_tot'].values

        com_col = tot_col_gosat - tot_col_iasi

        mini.append(np.nanmin(com_col))
        maxi.append(np.nanmax(com_col))

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

        plt.colorbar(label=r'$\mathrm{xN}_2\mathrm{O}$ in ppb')
        plt.title(f'gosat minus iasi {x_res}x{y_res} gridded data 2020-{m}')

        outpath = f'pictures/dif_{x_res}x{y_res}_{m}.png'
        plt.savefig(outpath)

    print(np.mean(mini), np.mean(maxi))