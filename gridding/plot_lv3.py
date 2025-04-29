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

def zonal_mean(var):
    var_zonal = np.zeros((var.shape[0], 1))
    for i in range(var.shape[0]):
        if np.isnan(var[i, :]).all():
            var_zonal[i] = np.nan
        else:
            var_zonal[i] = np.nanmean(var[i, :])
    return var_zonal


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
        ax.add_feature(cfeature.LAND, facecolor='gray')
        ax.add_feature(cfeature.OCEAN, facecolor='lightgray')

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

        mean_var1 = np.nanmean(var1)
        std_var1 = np.nanstd(var1)
        mean_var2 = np.nanmean(var2)
        std_var2 = np.nanstd(var2)

        lon_bnds = bnds(lon_iasi)
        lat_bnds = bnds(lat_iasi)

        # Create a figure with two subplots side by side
        fig, axs = plt.subplots(2, 2, figsize=(18, 10), subplot_kw={'projection': ccrs.PlateCarree()})

        # ---------------- Plot 1: iasi ----------------
        axs[0, 0].coastlines()
        axs[0, 0].add_feature(cfeature.BORDERS)
        axs[0, 0].add_feature(cfeature.LAND, facecolor='gray')
        axs[0, 0].add_feature(cfeature.OCEAN, facecolor='lightgray')

        # Add gridlines and set latitude/longitude labels
        gl00 = axs[0, 0].gridlines(draw_labels=True, linestyle='--', color='black')
        gl00.top_labels = gl00.right_labels = False
        gl00.xlabel_style = {'size': 10, 'color': 'black'}
        gl00.ylabel_style = {'size': 10, 'color': 'black'}

        # Plot the tot_col variable for IASI - GOSAT
        iasi_plot = axs[0, 0].pcolormesh(lon_bnds, lat_bnds, var1, transform=ccrs.PlateCarree(), cmap='seismic',
                                   vmin=vmin, vmax=vmax, norm=norm)
        axs[0, 0].set_title(f'IASI - GOSAT {x_res}x{y_res} grid 2020-{m} \n {mean_var1:.2f} \u00B1 {std_var1:.2f}')

        # ---------------- Plot 2: methode --------------
        axs[0, 1].coastlines()
        axs[0, 1].add_feature(cfeature.BORDERS)
        axs[0, 1].add_feature(cfeature.LAND, facecolor='gray')
        axs[0, 1].add_feature(cfeature.OCEAN, facecolor='lightgray')

        # Add gridlines and set latitude/longitude labels
        gl01 = axs[0, 1].gridlines(draw_labels=True, linestyle='--', color='black')
        gl01.top_labels = gl01.right_labels = False
        gl01.xlabel_style = {'size': 10, 'color': 'black'}
        gl01.ylabel_style = {'size': 10, 'color': 'black'}

        # Plot the count variable for IASI MET0 - GOSAT
        gosat_plot = axs[0, 1].pcolormesh(lon_bnds, lat_bnds, var2, transform=ccrs.PlateCarree(), cmap='seismic',
                                    vmin=vmin, vmax=vmax, norm=norm)
        axs[0, 1].set_title(f'Cor met{met} IASI - GOSAT {x_res}x{y_res} grid 2020-{m} \n {mean_var2:.2f} \u00B1 {std_var2:.2f}')

        # ---------------- Shared Colorbar ----------------
        # Create a shared colorbar for both subplots
        cbar = fig.colorbar(iasi_plot, ax=[axs[0, 0], axs[0, 1]], orientation='vertical', fraction=0.02, pad=0.04)
        cbar.set_label(r'$\mathrm{xN}_2\mathrm{O}$ in ppb')

        # ---------------- Plot 3: iasi count --------------
        iasi_count = ds_iasi['count'].values
        gosat_count = ds_gosat['count'].values

        imin, imax = [5, 6000]
        gmin, gmax = [5, 60]

        axs[1, 0].coastlines()
        axs[1, 0].add_feature(cfeature.BORDERS)
        axs[1, 0].add_feature(cfeature.LAND, facecolor='gray')
        axs[1, 0].add_feature(cfeature.OCEAN, facecolor='lightgray')

        # Add gridlines and set latitude/longitude labels
        gl10 = axs[1, 0].gridlines(draw_labels=True, linestyle='--', color='black')
        gl10.top_labels = gl10.right_labels = False
        gl10.xlabel_style = {'size': 10, 'color': 'black'}
        gl10.ylabel_style = {'size': 10, 'color': 'black'}

        # Plot the count of iasi data per grid cell
        icount_plot = axs[1, 0].pcolormesh(lon_bnds, lat_bnds, iasi_count, transform=ccrs.PlateCarree(), cmap='jet',
                                           vmin=imin, vmax=imax)
        axs[1, 0].set_title(f'IASI count {x_res}x{y_res} grid 2020-{m}')

        cbar = fig.colorbar(icount_plot, ax=axs[1, 0], orientation='vertical', fraction=0.02, pad=0.04)
        cbar.set_label(r'N')

        # ---------------- Plot 4: gosat count --------------
        axs[1, 1].coastlines()
        axs[1, 1].add_feature(cfeature.BORDERS)
        axs[1, 1].add_feature(cfeature.LAND, facecolor='gray')
        axs[1, 1].add_feature(cfeature.OCEAN, facecolor='lightgray')

        # Add gridlines and set latitude/longitude labels
        gl11 = axs[1, 1].gridlines(draw_labels=True, linestyle='--', color='black')
        gl11.top_labels = gl11.right_labels = False
        gl11.xlabel_style = {'size': 10, 'color': 'black'}
        gl11.ylabel_style = {'size': 10, 'color': 'black'}

        # Plot the count of gosat data per grid cell
        gcount_plot = axs[1, 1].pcolormesh(lon_bnds, lat_bnds, gosat_count, transform=ccrs.PlateCarree(), cmap='jet',
                                           vmin=gmin, vmax=gmax)
        axs[1, 1].set_title(f'IASI count {x_res}x{y_res} grid 2020-{m}')

        cbar = fig.colorbar(gcount_plot, ax=axs[1, 1], orientation='vertical', fraction=0.02, pad=0.04)
        cbar.set_label(r'N')

        outpath = f'pictures/combined_met{met}_{x_res}x{y_res}_{m}.png'
        plt.savefig(outpath)


def zonal_plot(x_res, y_res, th, dir_path, var, qf, vmin=None, vmax=None,
                  met_path='monthly_means/{}_met{}_{}_{}x{}_th{}_qf{}.nc', met=0):

    months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)  # norm for colorbar

    for m in months:
        ds_iasi = xr.open_dataset(dir_path.format('iasi', m, x_res, y_res, th, qf))
        ds_met = xr.open_dataset(met_path.format('iasi', met, m, x_res, y_res, th, qf))

        gosat_path = 'monthly_means/{}_{}_{}x{}_th{}.nc'
        ds_gosat = xr.open_dataset(gosat_path.format('gosat', m, x_res, y_res, th))

        lat_iasi = ds_iasi['lat_cen'].values

        variable_iasi = ds_iasi[var].values
        iasi_zonal = zonal_mean(variable_iasi)

        variable_met = ds_met[var].values
        met_zonal = zonal_mean(variable_met)

        variable_gosat = ds_gosat[var].values
        gosat_zonal = zonal_mean(variable_gosat)

        var1 = variable_iasi - variable_gosat
        var1_zonal = zonal_mean(var1)

        var2 = variable_met - variable_gosat
        var2_zonal = zonal_mean(var2)

        lon_bnds = [-180.0, 180.0]
        lat_bnds = bnds(lat_iasi)

        # # Create a figure with two subplots side by side
        # fig, axs = plt.subplots(2, 2, figsize=(18, 6), subplot_kw={'projection': ccrs.PlateCarree()})
        #
        # # ---------------- Plot 1: iasi ----------------
        # axs[0, 0].coastlines()
        # axs[0, 0].add_feature(cfeature.BORDERS)
        # axs[0, 0].add_feature(cfeature.LAND, facecolor='white')
        # axs[0, 0].add_feature(cfeature.OCEAN, facecolor='white')
        #
        # # Add gridlines and set latitude/longitude labels
        # gl1 = axs[0, 0].gridlines(draw_labels=True, linestyle='--', color='gray')
        # gl1.top_labels = gl1.right_labels = False
        # gl1.xlabel_style = {'size': 10, 'color': 'black'}
        # gl1.ylabel_style = {'size': 10, 'color': 'black'}
        #
        # # Plot the zonal tot_col from IASI
        # iasi_plot = axs[0, 0].pcolormesh(lon_bnds, lat_bnds, iasi_zonal, transform=ccrs.PlateCarree(), cmap='jet',
        #                            vmin=vmin, vmax=vmax)
        # axs[0,0].set_title(f'IASI zonal {x_res}x{y_res} grid 2020-{m}')
        #
        # # ---------------- Plot 2: methode --------------
        # axs[0, 1].coastlines()
        # axs[0, 1].add_feature(cfeature.BORDERS)
        # axs[0, 1].add_feature(cfeature.LAND, facecolor='white')
        # axs[0, 1].add_feature(cfeature.OCEAN, facecolor='white')
        #
        # # Add gridlines and set latitude/longitude labels
        # gl2 = axs[0, 1].gridlines(draw_labels=True, linestyle='--', color='gray')
        # gl2.top_labels = gl2.right_labels = False
        # gl2.xlabel_style = {'size': 10, 'color': 'black'}
        # gl2.ylabel_style = {'size': 10, 'color': 'black'}
        #
        # # Plot the zonal tot_col from IASI met
        # met_plot = axs[0, 1].pcolormesh(lon_bnds, lat_bnds, met_zonal, transform=ccrs.PlateCarree(), cmap='jet',
        #                             vmin=vmin, vmax=vmax)
        # axs[0, 1].set_title(f'Met{met} IASI zonal {x_res}x{y_res} grid 2020-{m}')
        #
        # # ---------------- Shared Colorbar ----------------
        # # Create a shared colorbar for both subplots
        # cbar = fig.colorbar(iasi_plot, ax=[axs[0, 0], axs[0, 1]], orientation='vertical', fraction=0.02, pad=0.04)
        # cbar.set_label(r'$\mathrm{xN}_2\mathrm{O}$ in ppb')
        #
        # # ---------------- Plot 3: iasi - gosat ----------------
        # axs[1, 0].coastlines()
        # axs[1, 0].add_feature(cfeature.BORDERS)
        # axs[1, 0].add_feature(cfeature.LAND, facecolor='white')
        # axs[1, 0].add_feature(cfeature.OCEAN, facecolor='white')
        #
        # # Add gridlines and set latitude/longitude labels
        # gl3 = axs[1, 0].gridlines(draw_labels=True, linestyle='--', color='gray')
        # gl3.top_labels = gl3.right_labels = False
        # gl3.xlabel_style = {'size': 10, 'color': 'black'}
        # gl3.ylabel_style = {'size': 10, 'color': 'black'}
        #
        # # Plot the zonal tot_col from IASI - GOSAT
        # dif_iasi_plot = axs[1, 0].pcolormesh(lon_bnds, lat_bnds, var1_zonal, transform=ccrs.PlateCarree(),
        #                                      cmap='seismic', vmin=vmin, vmax=vmax, norm=norm)
        # axs[1, 0].set_title(f'IASI - GOSAT zonal {x_res}x{y_res} grid 2020-{m}')
        #
        # # ---------------- Plot 2: methode --------------
        # axs[1, 1].coastlines()
        # axs[1, 1].add_feature(cfeature.BORDERS)
        # axs[1, 1].add_feature(cfeature.LAND, facecolor='white')
        # axs[1, 1].add_feature(cfeature.OCEAN, facecolor='white')
        #
        # # Add gridlines and set latitude/longitude labels
        # gl4 = axs[1, 1].gridlines(draw_labels=True, linestyle='--', color='gray')
        # gl4.top_labels = gl4.right_labels = False
        # gl4.xlabel_style = {'size': 10, 'color': 'black'}
        # gl4.ylabel_style = {'size': 10, 'color': 'black'}
        #
        # # Plot the zonal tot_col from IASI met - GOSAT
        # dif_met_plot = axs[1, 1].pcolormesh(lon_bnds, lat_bnds, var2_zonal, transform=ccrs.PlateCarree(),
        #                                     cmap='seismic', vmin=vmin, vmax=vmax, norm=norm)
        # axs[1, 1].set_title(f'Met{met} IASI - GOSAT zonal {x_res}x{y_res} grid 2020-{m}')
        #
        # # ---------------- Shared Colorbar ----------------
        # # Create a shared colorbar for both subplots
        # cbar = fig.colorbar(dif_iasi_plot, ax=[axs[1, 0], axs[1, 1]], orientation='vertical', fraction=0.02, pad=0.04)
        # cbar.set_label(r'$\mathrm{xN}_2\mathrm{O}$ in ppb')
        #
        # outpath = f'pictures/zonal_met{met}_{x_res}x{y_res}_{m}.png'
        # plt.savefig(outpath)

        # Graph plot

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))

        ax1.plot(iasi_zonal, lat_iasi, color='red')
        ax1.scatter(iasi_zonal, lat_iasi, color='red', label='True IASI')

        ax1.plot(met_zonal, lat_iasi, color='cyan')
        ax1.scatter(met_zonal, lat_iasi, color='cyan', label='Met0 IASI')

        ax1.plot(gosat_zonal, lat_iasi, color='purple')
        ax1.scatter(gosat_zonal, lat_iasi, color='purple', label='True GOSAT2')

        ax1.axhline(0, color='black', linewidth=1, linestyle='--')

        ax1.set_title(f'Zonal {x_res}x{y_res} grid 2020-{m}')
        ax1.set_ylabel('Latitude')
        ax1.set_ylim(-90, 90)
        ax1.set_xlabel('Zonal mean N2O in ppb')
        ax1.legend()

        ax2.plot(var1_zonal, lat_iasi, color='red')
        ax2.scatter(var1_zonal, lat_iasi, color='red', label='True IASI - GOSAT2')

        ax2.plot(var2_zonal, lat_iasi, color='cyan')
        ax2.scatter(var2_zonal, lat_iasi, color='cyan', label='Met0 IASI - GOSAT2')

        ax2.axhline(0, color='black', linewidth=1, linestyle='--')
        ax2.axvline(0, color='black', linewidth=1, linestyle='--')

        ax2.set_title(f'Dif Zonal {x_res}x{y_res} grid 2020-{m}')
        ax2.set_ylabel('Latitude')
        ax2.set_ylim(-90, 90)
        ax2.set_xlabel('Zonal mean N2O in ppb')
        ax2.legend()

        outpath = f'pictures/zonal_met{met}_{x_res}x{y_res}_{m}.png'
        plt.savefig(outpath)


# def iasi_dif_plot(x_res, y_res, th, dir_path, qf, vmin=None, vmax=None,
#                   met_path='monthly_means/{}_met{}_{}_{}x{}_th{}_qf{}.nc'):
#
#     months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
#
#     for m in months:
#         ds_iasi = xr.open_dataset(dir_path.format('iasi', m, x_res, y_res, th, qf))
#         ds_met0 = xr.open_dataset(met_path.format('iasi', '0', m, x_res, y_res, th, qf))
#         ds_met1 = xr.open_dataset(met_path.format('iasi', '1', m, x_res, y_res, th, qf))
#         ds_met2 = xr.open_dataset(met_path.format('iasi', '2', m, x_res, y_res, th, qf))
#
#         lon_iasi = ds_iasi['lon_cen'].values
#         lat_iasi = ds_iasi['lat_cen'].values
#         variable_iasi = ds_iasi['mean_tot'].values
#         variable_met0 = ds_met0['mean_tot'].values
#         variable_met1 = ds_met1['mean_tot'].values
#         variable_met2 = ds_met2['mean_tot'].values
#
#         var0 = variable_iasi - variable_met0
#         var1 = variable_iasi - variable_met1
#         var2 = variable_iasi - variable_met2
#
#         lon_bnds = bnds(lon_iasi)
#         lat_bnds = bnds(lat_iasi)
#
#         # Create a figure with two subplots side by side
#         fig, axs = plt.subplots(2, 2, figsize=(18, 12), subplot_kw={'projection': ccrs.PlateCarree()})
#
#         # ---------------- Plot 0: Met 0 ----------------
#         axs[0, 0].coastlines()
#         axs[0, 0].add_feature(cfeature.BORDERS)
#         axs[0, 0].add_feature(cfeature.LAND, facecolor='white')
#         axs[0, 0].add_feature(cfeature.OCEAN, facecolor='white')
#
#         # Add gridlines and set latitude/longitude labels
#         gl0 = axs[0, 0].gridlines(draw_labels=True, linestyle='--', color='gray')
#         gl0.top_labels = gl0.right_labels = False
#         gl0.xlabel_style = {'size': 10, 'color': 'black'}
#         gl0.ylabel_style = {'size': 10, 'color': 'black'}
#
#         # Plot
#         plot0 = axs[0, 0].pcolormesh(lon_bnds, lat_bnds, var0, transform=ccrs.PlateCarree(), cmap='jet',
#                                    vmin=vmin, vmax=vmax)
#         axs[0, 0].set_title(f'IASI - MET 0 {x_res}x{y_res} grid 2020-{m}')
#
#         # ---------------- Plot 1: Met 1 ----------------
#         axs[0, 1].coastlines()
#         axs[0, 1].add_feature(cfeature.BORDERS)
#         axs[0, 1].add_feature(cfeature.LAND, facecolor='white')
#         axs[0, 1].add_feature(cfeature.OCEAN, facecolor='white')
#
#         # Add gridlines and set latitude/longitude labels
#         gl1 = axs[0, 1].gridlines(draw_labels=True, linestyle='--', color='gray')
#         gl1.top_labels = gl1.right_labels = False
#         gl1.xlabel_style = {'size': 10, 'color': 'black'}
#         gl1.ylabel_style = {'size': 10, 'color': 'black'}
#
#         # Plot
#         plot1 = axs[0, 1].pcolormesh(lon_bnds, lat_bnds, var1, transform=ccrs.PlateCarree(), cmap='jet',
#                                    vmin=vmin, vmax=vmax)
#         axs[0, 1].set_title(f'IASI - MET 1 {x_res}x{y_res} grid 2020-{m}')
#
#         # ---------------- Plot 2: Met 2 ----------------
#         axs[1, 0].coastlines()
#         axs[1, 0].add_feature(cfeature.BORDERS)
#         axs[1, 0].add_feature(cfeature.LAND, facecolor='white')
#         axs[1, 0].add_feature(cfeature.OCEAN, facecolor='white')
#
#         # Add gridlines and set latitude/longitude labels
#         gl2 = axs[1, 0].gridlines(draw_labels=True, linestyle='--', color='gray')
#         gl2.top_labels = gl2.right_labels = False
#         gl2.xlabel_style = {'size': 10, 'color': 'black'}
#         gl2.ylabel_style = {'size': 10, 'color': 'black'}
#
#         # Plot
#         plot2 = axs[1, 0].pcolormesh(lon_bnds, lat_bnds, var2, transform=ccrs.PlateCarree(), cmap='jet',
#                                    vmin=vmin, vmax=vmax)
#         axs[1, 0].set_title(f'IASI - MET 2 {x_res}x{y_res} grid 2020-{m}')
#
#         # ---------------- Shared Colorbar ----------------
#         # Create a shared colorbar for both subplots
#         cbar = fig.colorbar(plot0, ax=axs, orientation='vertical', fraction=0.02, pad=0.04)
#         cbar.set_label(r'$\mathrm{xN}_2\mathrm{O}$ in ppb')
#
#         outpath = f'pictures/iasi_dif_{x_res}x{y_res}_{m}.png'
#         plt.savefig(outpath)


# def plot_dif(x_res, y_res, th, dir_path, vmin=None, vmax=None):
#     """Not yet updated for newer data format"""
#
#     months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
#
#     mini = []
#     maxi = []
#
#     for m in months:
#         ds_iasi = xr.open_dataset(dir_path.format('iasi', m, x_res, y_res, th))
#         ds_gosat = xr.open_dataset(dir_path.format('gosat', m, x_res, y_res, th))
#
#         lon = ds_iasi['lon_cen'].values
#         lat = ds_iasi['lat_cen'].values
#         tot_col_iasi = ds_iasi['mean_tot'].values
#         tot_col_gosat = ds_gosat['mean_tot'].values
#
#         com_col = tot_col_gosat - tot_col_iasi
#
#         mini.append(np.nanmin(com_col))
#         maxi.append(np.nanmax(com_col))
#
#         lon_bnds = bnds(lon)
#         lat_bnds = bnds(lat)
#
#         plt.figure(figsize=(12, 6))
#         ax = plt.axes(projection=ccrs.PlateCarree())
#
#         ax.coastlines()
#         ax.add_feature(cfeature.BORDERS)
#         ax.add_feature(cfeature.LAND, facecolor='white')
#         ax.add_feature(cfeature.OCEAN, facecolor='white')
#
#         # Add gridlines and set latitude/longitude labels
#         gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray')
#         gl.top_labels = gl.right_labels = False
#         gl.xlabel_style = {'size': 10, 'color': 'black'}
#         gl.ylabel_style = {'size': 10, 'color': 'black'}
#
#         plt.pcolormesh(lon_bnds, lat_bnds, com_col, transform=ccrs.PlateCarree(), cmap='jet'
#                        , vmin=vmin, vmax=vmax)
#
#         plt.colorbar(label=r'$\mathrm{xN}_2\mathrm{O}$ in ppb')
#         plt.title(f'gosat minus iasi {x_res}x{y_res} gridded data 2020-{m}')
#
#         outpath = f'pictures/dif_{x_res}x{y_res}_{m}.png'
#         plt.savefig(outpath)
#
#     print(np.mean(mini), np.mean(maxi))