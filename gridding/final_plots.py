from calendar import month

import matplotlib as mpl
import matplotlib.pyplot as plt

import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib.colors as mcolors

from matplotlib.ticker import FuncFormatter


def correction_example(month):

    mpl.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.size": 10.0,  # <- from FONTSIZE in your log
        "axes.labelsize": 9.0,
        "axes.titlesize": 10.95,  # actually does matter
        "legend.fontsize": 8,  # ~\footnotesize
        "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage{lmodern}",
        # add packages you actually use in labels, e.g.:
        # r"\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{siunitx}\usepackage{mhchem}"
    })

    FIGSIZE = set_size(width_fraction=1.3, height_fraction=0.66)
    cbar_fontsize = 9

    dir_path = 'monthly_means/{}_{:02d}_{}x{}_th{}_qf{}.nc'
    met_path = 'monthly_means/{}_met{}_{:02d}_{}x{}_th{}_qf{}.nc'

    #month = 2
    x_res = 5
    y_res = 5
    th_iasi = 10
    th_gosat = 3
    qf = 3
    var = 'mean_tot'
    met = 0

    ds_iasi = xr.open_dataset(dir_path.format('iasi', month, x_res, y_res, th_iasi, qf))
    ds_met = xr.open_dataset(met_path.format('iasi', met, month, x_res, y_res, th_iasi, qf))

    gosat_path = 'monthly_means/{}_{:02d}_{}x{}_th{}.nc'
    ds_gosat = xr.open_dataset(gosat_path.format('gosat', month, x_res, y_res, th_gosat))

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

    # keep copies for the absolute-difference row
    var1_abs = var1.copy()   # IASI - GOSAT
    var2_abs = var2.copy()   # Cor met{met} IASI - GOSAT

    # prepare zonal-mean-removed fields for the second row
    var1_zm  = var1_abs - zonal_mean(var1_abs)
    var2_zm  = var2_abs - zonal_mean(var2_abs)

    # norms per row
    norm_row1 = mcolors.TwoSlopeNorm(vmin=-30, vcenter=0, vmax=30)
    norm_row2 = mcolors.TwoSlopeNorm(vmin=-10, vcenter=0, vmax=10)

    # new combined 2x2 figure
    fig, axs = plt.subplots(2, 2, figsize=FIGSIZE, subplot_kw={'projection': ccrs.PlateCarree()},
    constrained_layout=True)

    def style_ax(ax, pos):
        ax.coastlines()
        #ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.LAND, facecolor='gray')
        ax.add_feature(cfeature.OCEAN, facecolor='lightgray')
        gl = ax.gridlines(draw_labels=True, linestyle='--', color='black')
        if pos == 'left':
            gl.top_labels = gl.right_labels = False
        else:
            gl.top_labels = gl.right_labels = gl.left_labels = False
        gl.xlabel_style = {'size': 8, 'color': 'black'}
        gl.ylabel_style = {'size': 8, 'color': 'black'}

    # ---- Row 1: absolute differences, [-30, 30] ----
    style_ax(axs[0, 0], pos='left')
    pcm00 = axs[0, 0].pcolormesh(lon_bnds, lat_bnds, var1_abs, transform=ccrs.PlateCarree(),
                                 cmap='seismic', norm=norm_row1)
    axs[0, 0].set_title(rf'IASI - GOSAT-2 {x_res}$^\circ$$\times${y_res}$^\circ$ grid 2020-{month}'
                        + '\n'
                        + f'global mean {mean_var1:.2f} \u00B1 {std_var1:.2f}')

    style_ax(axs[0, 1], pos='right')
    pcm01 = axs[0, 1].pcolormesh(lon_bnds, lat_bnds, var2_abs, transform=ccrs.PlateCarree(),
                                 cmap='seismic', norm=norm_row1)
    axs[0, 1].set_title(rf'IASI (a-priori) - GOSAT-2 {x_res}$^\circ$$\times${y_res}$^\circ$ grid 2020-{month}'
                        + '\n'
                        f' global mean {mean_var2:.2f} \u00B1 {std_var2:.2f}')

    # shared colorbar for top row
    cbar1 = fig.colorbar(pcm00, ax=axs[0, :], orientation='horizontal', fraction=0.02, pad=0.04, aspect=80)
    cbar1.set_label(r'$\mathrm{xN}_2\mathrm{O}$ in ppb', fontsize=cbar_fontsize)

    # ---- Row 2: zonal-mean-removed differences, [-10, 10] ----
    style_ax(axs[1, 0], pos='left')
    pcm10 = axs[1, 0].pcolormesh(lon_bnds, lat_bnds, var1_zm, transform=ccrs.PlateCarree(),
                                 cmap='seismic', norm=norm_row2)
    axs[1, 0].set_title(rf'IASI - GOSAT-2 - ZONAL {x_res}$^\circ$$\times${y_res}$^\circ$ grid 2020-{month}')

    style_ax(axs[1, 1], pos='right')
    pcm11 = axs[1, 1].pcolormesh(lon_bnds, lat_bnds, var2_zm, transform=ccrs.PlateCarree(),
                                 cmap='seismic', norm=norm_row2)
    axs[1, 1].set_title(rf'IASI (a-priori) - GOSAT-2 - ZONAL {x_res}$^\circ$$\times${y_res}$^\circ$ grid 2020-{month}')

    # shared colorbar for bottom row
    cbar2 = fig.colorbar(pcm10, ax=axs[1, :], orientation='horizontal', fraction=0.02, pad=0.04, aspect=80)
    cbar2.set_label(r'$\mathrm{xN}_2\mathrm{O}$ in ppb', fontsize=cbar_fontsize)

    fig.set_constrained_layout_pads(
        w_pad=1 / 72,  # padding around axes (inches, ~1pt)
        h_pad=1 / 72,
        wspace=0.02,  # space between columns (fraction of axes)
        hspace=0.00  # space between rows
    )

    #outpath_combined = f'pictures/differences_combined_{x_res}x{y_res}_{month}.png'
    #plt.savefig(outpath_combined, dpi=300)

    base = f"pictures/final_plots/iasi_correction_{x_res}x{y_res}_{month}"
    fig.savefig(base + ".pdf")  # vector PDF (great for LaTeX)


def tot_mean_example(month):

    mpl.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.size": 10.0,  # <- from FONTSIZE in your log
        "axes.labelsize": 9.0,
        "axes.titlesize": 10.95,  # actually does matter
        "legend.fontsize": 8,  # ~\footnotesize
        "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage{lmodern}",
        # add packages you actually use in labels, e.g.:
        # r"\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{siunitx}\usepackage{mhchem}"
    })

    FIGSIZE = set_size(width_fraction=1.3, height_fraction=0.63)
    cbar_fontsize = 9

    dir_path = 'monthly_means/{}_{:02d}_{}x{}_th{}_qf{}.nc'

    #month = 3
    x_res = 5
    y_res = 5
    th_iasi = 10
    th_gosat = 3
    th_cams = 0
    qf = 3
    var = 'mean_tot'
    met = 0
    vmin = 290
    vmax = 340

    ds_iasi = xr.open_dataset(dir_path.format('iasi', month, x_res, y_res, th_iasi, qf))

    gosat_path = 'monthly_means/{}_{:02d}_{}x{}_th{}.nc'
    ds_gosat = xr.open_dataset(gosat_path.format('gosat', month, x_res, y_res, th_gosat))

    ds_cams = xr.open_dataset(gosat_path.format('cams', month, x_res, y_res, th_cams))
    ds_cams_iasi = xr.open_dataset(gosat_path.format('cams_iasi', month, x_res, y_res, th_cams))

    lon_iasi = ds_iasi['lon_cen'].values
    lat_iasi = ds_iasi['lat_cen'].values
    variable_iasi = ds_iasi[var].values
    variable_gosat = ds_gosat[var].values
    variable_cams = ds_cams[var].values
    variable_cams_i = ds_cams_iasi[var].values

    if month == 3:
        print(np.nanmin(variable_iasi), np.nanmax(variable_iasi))
        print(np.nanmin(variable_gosat), np.nanmax(variable_gosat))
        print(np.nanmin(variable_cams), np.nanmax(variable_cams))
        print(np.nanmin(variable_cams_i), np.nanmax(variable_cams_i))

    iasi_zonal = zonal_mean(variable_iasi)
    gosat_zonal = zonal_mean(variable_gosat)
    cams_zonal = zonal_mean(variable_cams)
    cams_i_zonal = zonal_mean(variable_cams_i)

    # iasi_grad = gradiant(iasi_zonal)
    # gosat_grad = gradiant(gosat_zonal)
    # cams_grad = gradiant(cams_zonal)
    # cams_igrad = gradiant(cams_i_zonal)
    #
    # print(iasi_grad)
    # print(gosat_grad)
    # print(cams_grad)
    # print(cams_igrad)

    lon_bnds = bnds(lon_iasi)
    lat_bnds = bnds(lat_iasi)

    ######################## Plot all products #######################################
    # new combined 2x2 figure
    fig, axs = plt.subplots(2, 2, figsize=FIGSIZE, subplot_kw={'projection': ccrs.PlateCarree()},
    constrained_layout=True)

    def style_ax(ax, pos):
        ax.coastlines()
        #ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.LAND, facecolor='gray')
        ax.add_feature(cfeature.OCEAN, facecolor='lightgray')
        gl = ax.gridlines(draw_labels=True, linestyle='--', color='black')
        if pos == 'left':
            gl.top_labels = gl.right_labels = False
        else:
            gl.top_labels = gl.right_labels = gl.left_labels = False
        gl.xlabel_style = {'size': 8, 'color': 'black'}
        gl.ylabel_style = {'size': 8, 'color': 'black'}

    # ---- Row 1: IASI GOSAT ----
    style_ax(axs[0, 0], pos='left')
    pcm00 = axs[0, 0].pcolormesh(lon_bnds, lat_bnds, variable_iasi, transform=ccrs.PlateCarree(), cmap='jet'
                       , vmin=vmin, vmax=vmax)
    axs[0, 0].set_title(rf'IASI total columns {x_res}$^\circ$$\times${y_res}$^\circ$ grid 2020-{month:02d}')

    style_ax(axs[0, 1], pos='right')
    pcm01 = axs[0, 1].pcolormesh(lon_bnds, lat_bnds, variable_gosat, transform=ccrs.PlateCarree(), cmap='jet'
                       , vmin=vmin, vmax=vmax)
    axs[0, 1].set_title(rf'GOSAT-2 total columns {x_res}$^\circ$$\times${y_res}$^\circ$ grid 2020-{month:02d}')

    # ---- Row 2: CAMS ----
    style_ax(axs[1, 0], pos='left')
    pcm10 = axs[1, 0].pcolormesh(lon_bnds, lat_bnds, variable_cams_i, transform=ccrs.PlateCarree(), cmap='jet'
                       , vmin=vmin, vmax=vmax)
    axs[1, 0].set_title(rf'CAMS (IASI AVK) total columns {x_res}$^\circ$$\times${y_res}$^\circ$ grid 2020-{month:02d}')

    style_ax(axs[1, 1], pos='right')
    pcm11 = axs[1, 1].pcolormesh(lon_bnds, lat_bnds, variable_cams, transform=ccrs.PlateCarree(), cmap='jet'
                       , vmin=vmin, vmax=vmax)
    axs[1, 1].set_title(rf'CAMS total columns {x_res}$^\circ$$\times${y_res}$^\circ$ grid 2020-{month:02d}')

    # shared colorbar
    cbar = fig.colorbar(pcm10, ax=axs[:, :], orientation='horizontal', fraction=0.02, pad=0.04, aspect=50)
    cbar.set_label(r'$\mathrm{xN}_2\mathrm{O}$ in ppb', fontsize=cbar_fontsize)

    fig.set_constrained_layout_pads(
        w_pad=1 / 72,  # padding around axes (inches, ~1pt)
        h_pad=1 / 72,
        wspace=0.02,  # space between columns (fraction of axes)
        hspace=0.00  # space between rows
    )

    #plt.show()

    base = f"pictures/final_plots/total_columns_{x_res}x{y_res}_{month}"
    fig.savefig(base + ".pdf")  # vector PDF (great for LaTeX)
    #fig.savefig(f"pictures/total_columns_{x_res}x{y_res}_{month}.png", dpi=300)

    ######################## Plot IASI and CAMS - ZONAL #######################################

    variable_iasi_zonal = variable_iasi - iasi_zonal
    variable_cams_i_zonal = variable_cams_i - cams_i_zonal

    norm = mcolors.TwoSlopeNorm(vmin=-10, vcenter=0, vmax=10)  # norm for colorbar

    FIGSIZE_zonal = set_size(width_fraction=1.3, height_fraction=0.32)

    fig, axs = plt.subplots(1, 2, figsize=FIGSIZE_zonal, subplot_kw={'projection': ccrs.PlateCarree()},
                            constrained_layout=True)

    def style_ax(ax, pos):
        ax.coastlines()
        # ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.LAND, facecolor='gray')
        ax.add_feature(cfeature.OCEAN, facecolor='lightgray')
        gl = ax.gridlines(draw_labels=True, linestyle='--', color='black')
        if pos == 'left':
            gl.top_labels = gl.right_labels = False
        else:
            gl.top_labels = gl.right_labels = gl.left_labels = False
        gl.xlabel_style = {'size': 8, 'color': 'black'}
        gl.ylabel_style = {'size': 8, 'color': 'black'}

    # ---- Row 1: IASI GOSAT ----
    style_ax(axs[0], pos='left')
    pcm00 = axs[0].pcolormesh(lon_bnds, lat_bnds, variable_iasi_zonal, transform=ccrs.PlateCarree(), cmap='seismic'
                                 , norm=norm)
    axs[0].set_title(rf'IASI minus zonal mean {x_res}$^\circ$$\times${y_res}$^\circ$ grid 2020-{month:02d}')

    style_ax(axs[1], pos='right')
    pcm01 = axs[1].pcolormesh(lon_bnds, lat_bnds, variable_cams_i_zonal, transform=ccrs.PlateCarree(), cmap='seismic'
                                 , norm=norm)
    axs[1].set_title(rf'CAMS (IASI AVK) minus zonal mean {x_res}$^\circ$$\times${y_res}$^\circ$ grid 2020-{month:02d}')

    # shared colorbar
    cbar = fig.colorbar(pcm00, ax=axs[:], orientation='horizontal', fraction=0.02, pad=0.04, aspect=50)
    cbar.set_label(r'$\mathrm{xN}_2\mathrm{O}$ in ppb', fontsize=cbar_fontsize)

    fig.set_constrained_layout_pads(
        w_pad=1 / 72,  # padding around axes (inches, ~1pt)
        h_pad=1 / 72,
        wspace=0.02,  # space between columns (fraction of axes)
        hspace=0.00  # space between rows
    )

    base = f"pictures/final_plots/total_columns-zonal_{x_res}x{y_res}_{month}"
    # fig.savefig(base + ".pdf")  # vector PDF (great for LaTeX)
    fig.savefig(f"pictures/total_columns-zonal_{x_res}x{y_res}_{month}.png", dpi=300)


def seasonal_plots():

    mpl.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.size": 8.0,  # <- from FONTSIZE in your log
        "axes.labelsize": 9.0,
        "axes.titlesize": 10.95,  # actually does matter
        "legend.fontsize": 8,  # ~\footnotesize
        "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage{lmodern}",
        # add packages you actually use in labels, e.g.:
        # r"\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{siunitx}\usepackage{mhchem}"
    })

    FIGSIZE = set_size(width_fraction=1.3, height_fraction=0.63)
    cbar_fontsize = 9

    var = 'mean_tot'
    x_res = 5
    y_res = 5
    path = 'monthly_means/{}_{:02}_{}x{}_th{}.nc'
    path_iasi = 'monthly_means/{}_{:02}_{}x{}_th{}_qf{}.nc'
    th_iasi = 10
    th_gosat = 3
    th_cams = 0

    # Load first month just to get lat array length
    ds_iasi_first = xr.open_dataset(path_iasi.format('iasi', 1, x_res, y_res, th_iasi, 3))
    lat_cen = ds_iasi_first['lat_cen'].values
    n_lat = len(lat_cen)

    # Create array: rows = lat, cols = 12 months
    all_iasi_zonal = np.full((n_lat, 12), np.nan)
    all_gosat_zonal = np.full((n_lat, 12), np.nan)
    all_cams_zonal = np.full((n_lat, 12), np.nan)
    all_cams_iasi_zonal = np.full((n_lat, 12), np.nan)

    for m in range(1, 13):

        ds_iasi = xr.open_dataset(path_iasi.format('iasi', m, x_res, y_res, th_iasi, 3))
        ds_gosat = xr.open_dataset(path.format('gosat', m, x_res, y_res, th_gosat))
        ds_cams = xr.open_dataset(path.format('cams', m, x_res, y_res, th_cams))
        ds_cams_iasi = xr.open_dataset(path.format('cams_iasi', m, x_res, y_res, th_cams))

        variable_iasi = ds_iasi[var].values
        iasi_zonal = zonal_mean(variable_iasi)
        all_iasi_zonal[:, m - 1] = iasi_zonal[:,0]

        variable_gosat = ds_gosat[var].values
        gosat_zonal = zonal_mean(variable_gosat)
        all_gosat_zonal[:, m - 1] = gosat_zonal[:,0]

        variable_cams = ds_cams[var].values
        cams_zonal = zonal_mean(variable_cams)
        all_cams_zonal[:, m - 1] = cams_zonal[:,0]

        variable_cams_iasi = ds_cams_iasi[var].values
        cams_iasi_zonal = zonal_mean(variable_cams_iasi)
        all_cams_iasi_zonal[:, m - 1] = cams_iasi_zonal[:,0]

    ################# Plot Seasonal: Plot original data #################
    months = np.arange(1, 13)  # 1 to 12
    M, L = np.meshgrid(months, lat_cen)

    all_data = [
        (all_iasi_zonal, "IASI Zonal Mean"),
        (all_gosat_zonal, "GOSAT-2 Zonal Mean"),
        (all_cams_iasi_zonal, "CAMS (IASI AVK) Zonal Mean"),
        (all_cams_zonal, "CAMS Zonal Mean")
    ]

    vmin = min(np.nanmin(d[0]) for d in all_data)
    vmax = max(np.nanmax(d[0]) for d in all_data)

    fig, axs = plt.subplots(2, 2, figsize=FIGSIZE, sharex=True, sharey=True,
                            constrained_layout=True)
    cs = []

    yticks = np.arange(-60, 61, 30)  # -90, -60, -30, 0, 30, 60, 90

    for ax, (data, title) in zip(axs.ravel(), all_data):
        c = ax.pcolormesh(M, L, data, cmap='jet', vmin=vmin, vmax=vmax, shading='auto')
        ax.set_title(title)
        if ax in axs[:, 0]:
            ax.set_ylim(-90, 90)
            ax.set_yticks(yticks)
            ax.yaxis.set_major_formatter(FuncFormatter(lat_fmt))
            #ax.set_ylabel("Latitude")
        cs.append(c)
    for i in range(2):
        axs[1, i].set_xlabel("Month")

    xticks = [1, 3, 6, 9, 12]

    for ax in axs.ravel():  # all 4 axes
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks)
        ax.tick_params(axis='x', labelbottom=True)  # force labels to show

    # Add overall title
    #fig.suptitle("Seasonal Zonal Means for 2020", fontsize=16)

    # Add colorbar to the right of all subplots
    cbar = fig.colorbar(cs[0], ax=axs, orientation="horizontal", fraction=0.05, pad=0.04)
    cbar.set_label(r"mean $\mathrm{xN}_2\mathrm{O}$ in ppb")

    fig.set_constrained_layout_pads(
        w_pad=1 / 72,  # padding around axes (inches, ~1pt)
        h_pad=1 / 72,
        wspace=0.02,  # space between columns (fraction of axes)
        hspace=0.0  # space between rows
    )

    #plt.show()
    base = f'pictures/seasonal_{x_res}x{y_res}'
    #fig.savefig(base + ".pdf")  # vector PDF (great for LaTeX)
    fig.savefig(base + ".png", dpi=300)

    ################# Plot Seasonal: Plot diffs #################

    FIGSIZE_row = set_size(width_fraction=1.3, height_fraction=0.4)

    diff_data = [
        (all_iasi_zonal - all_cams_iasi_zonal, "IASI - CAMS (IASI AVK) Zonal Mean"),
        (all_gosat_zonal - all_cams_zonal, "GOSAT-2 - CAMS Zonal Mean")
    ]

    vmin = -10
    vmax = 10

    # Create diverging colormap: blue-white-red
    cmap = plt.cm.bwr.copy()
    cmap.set_bad("gray")  # gray background for NaNs

    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)  # norm for colorbar

    fig, axs = plt.subplots(1, 2, figsize=FIGSIZE_row, sharex=True, sharey=True,
                            constrained_layout=True)
    cs = []

    yticks = np.arange(-60, 61, 30)  # -90, -60, -30, 0, 30, 60, 90

    for ax, (data, title) in zip(axs.ravel(), diff_data):
        c = ax.pcolormesh(M, L, data, cmap=cmap, norm=norm, shading='auto')
        ax.set_title(title)
        if ax == axs[0]:
            ax.set_ylim(-90, 90)
            ax.set_yticks(yticks)
            ax.yaxis.set_major_formatter(FuncFormatter(lat_fmt))
            #ax.set_ylabel("Latitude")
        cs.append(c)
        ax.set_xlabel("Month")

    xticks = [1, 3, 6, 9, 12]

    for ax in axs.ravel():  # all 4 axes
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks)
        ax.tick_params(axis='x', labelbottom=True)  # force labels to show

    # Add overall title
    #fig.suptitle("Seasonal Zonal Means for 2020", fontsize=16)

    # Add colorbar to the right of all subplots
    cbar = fig.colorbar(cs[0], ax=axs, orientation="horizontal", fraction=0.05, pad=0.04)
    cbar.set_label(r"mean $\mathrm{xN}_2\mathrm{O}$ in ppb")

    fig.set_constrained_layout_pads(
        w_pad=1 / 72,  # padding around axes (inches, ~1pt)
        h_pad=1 / 72,
        wspace=0.02,  # space between columns (fraction of axes)
        hspace=0.0  # space between rows
    )

    #plt.show()
    base = f'pictures/seasonal_diff_{x_res}x{y_res}'
    # fig.savefig(base + ".pdf")  # vector PDF (great for LaTeX)
    fig.savefig(base + ".png", dpi=300)

    ################# Plot Seasonal: Plot diffs (IASI - GOSAT) #################

    FIGSIZE_one = set_size(width_fraction=0.75, height_fraction=0.4)

    # use only the first diff
    diff0_data, diff0_title = (all_iasi_zonal - all_gosat_zonal, "IASI - GOSAT-2 Zonal Mean")

    vmin, vmax = -10, 10

    # diverging colormap centered at 0
    cmap = plt.cm.bwr.copy()
    cmap.set_bad("gray")  # NaNs as gray
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    # single subplot
    fig, ax = plt.subplots(1, 1, figsize=FIGSIZE_one, constrained_layout=True)

    # ticks/labels
    yticks = np.arange(-60, 61, 30)

    # plot
    c = ax.pcolormesh(M, L, diff0_data, cmap=cmap, norm=norm, shading='auto')
    ax.set_title(diff0_title)
    ax.set_xlabel("Month")
    ax.set_ylim(-90, 90)
    ax.set_yticks(yticks)
    ax.yaxis.set_major_formatter(FuncFormatter(lat_fmt))

    # x ticks
    xticks = [1, 3, 6, 9, 12]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks)
    ax.tick_params(axis='x', labelbottom=True)

    # bottom colorbar
    cbar = fig.colorbar(c, ax=ax, orientation="horizontal", fraction=0.05, pad=0.04)
    cbar.set_label(r"mean $\mathrm{xN}_2\mathrm{O}$ in ppb")

    # save
    base = f'pictures/seasonal_diff_{x_res}x{y_res}'
    fig.savefig(base + "_iasi_gosat.png", dpi=300)


def proxy_plots(month):

    mpl.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.size": 10.0,  # <- from FONTSIZE in your log
        "axes.labelsize": 9.0,
        "axes.titlesize": 10.95,  # actually does matter
        "legend.fontsize": 8,  # ~\footnotesize
        "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage{lmodern}",
        # add packages you actually use in labels, e.g.:
        # r"\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{siunitx}\usepackage{mhchem}"
    })

    FIGSIZE = set_size(width_fraction=1.3, height_fraction=0.63)
    cbar_fontsize = 9

    dir_path = 'monthly_means/{}_{:02d}_{}x{}_th{}_qf{}.nc'

    #month = 3
    x_res = 5
    y_res = 5
    th_iasi = 10
    th_gosat = 3
    th_cams = 0
    qf = 3
    var = 'mean_tot'
    met = 0
    vmin = 280
    vmax = 340

    ds_iasi = xr.open_dataset(dir_path.format('iasi', month, x_res, y_res, th_iasi, qf))

    gosat_path = 'monthly_means/{}_{:02d}_{}x{}_th{}.nc'
    ds_gosat = xr.open_dataset(gosat_path.format('gosat', month, x_res, y_res, th_gosat))

    ds_cams = xr.open_dataset(gosat_path.format('cams', month, x_res, y_res, th_cams))
    ds_cams_iasi = xr.open_dataset(gosat_path.format('cams_iasi', month, x_res, y_res, th_cams))

    lon_iasi = ds_iasi['lon_cen'].values
    lat_iasi = ds_iasi['lat_cen'].values
    variable_iasi = ds_iasi[var].values
    variable_gosat = ds_gosat[var].values
    variable_cams = ds_cams[var].values
    variable_cams_i = ds_cams_iasi[var].values

    iasi_zonal = zonal_mean(variable_iasi)
    gosat_zonal = zonal_mean(variable_gosat)
    cams_zonal = zonal_mean(variable_cams)
    cams_i_zonal = zonal_mean(variable_cams_i)

    lon_bnds = bnds(lon_iasi)
    lat_bnds = bnds(lat_iasi)

    iasi_proxy = variable_iasi - iasi_zonal
    gosat_proxy = variable_gosat - gosat_zonal
    cams_proxy = variable_cams - cams_zonal
    camsI_proxy = variable_cams_i - cams_i_zonal

    ######################## Proxy diffs #######################################
    proxy_dif1 = iasi_proxy - camsI_proxy
    proxy_dif2 = gosat_proxy - cams_proxy

    norm = mcolors.TwoSlopeNorm(vmin=-10, vcenter=0, vmax=10)  # norm for colorbar

    FIGSIZE_zonal = set_size(width_fraction=1.3, height_fraction=0.32)

    fig, axs = plt.subplots(1, 2, figsize=FIGSIZE_zonal, subplot_kw={'projection': ccrs.PlateCarree()},
                            constrained_layout=True)

    def style_ax(ax, pos):
        ax.coastlines()
        # ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.LAND, facecolor='gray')
        ax.add_feature(cfeature.OCEAN, facecolor='lightgray')
        gl = ax.gridlines(draw_labels=True, linestyle='--', color='black')
        if pos == 'left':
            gl.top_labels = gl.right_labels = False
        else:
            gl.top_labels = gl.right_labels = gl.left_labels = False
        gl.xlabel_style = {'size': 8, 'color': 'black'}
        gl.ylabel_style = {'size': 8, 'color': 'black'}

    # ---- Row 1: IASI GOSAT ----
    style_ax(axs[0], pos='left')
    pcm00 = axs[0].pcolormesh(lon_bnds, lat_bnds, proxy_dif1, transform=ccrs.PlateCarree(), cmap='seismic'
                              , norm=norm)
    axs[0].set_title(rf'IASI - CAMS (IASI AVK) proxys {x_res}$^\circ$$\times${y_res}$^\circ$ grid 2020-{month:02d}')

    style_ax(axs[1], pos='right')
    pcm01 = axs[1].pcolormesh(lon_bnds, lat_bnds, proxy_dif2, transform=ccrs.PlateCarree(), cmap='seismic'
                              , norm=norm)
    axs[1].set_title(rf'GOSAT-2 - CAMS proxys {x_res}$^\circ$$\times${y_res}$^\circ$ grid 2020-{month:02d}')

    # shared colorbar
    cbar = fig.colorbar(pcm00, ax=axs[:], orientation='horizontal', fraction=0.02, pad=0.04, aspect=50)
    cbar.set_label(r'$\mathrm{xN}_2\mathrm{O}$ in ppb', fontsize=cbar_fontsize)

    fig.set_constrained_layout_pads(
        w_pad=1 / 72,  # padding around axes (inches, ~1pt)
        h_pad=1 / 72,
        wspace=0.02,  # space between columns (fraction of axes)
        hspace=0.00  # space between rows
    )

    base = f"pictures/proxy_{x_res}x{y_res}_{month}"
    # fig.savefig(base + ".pdf")  # vector PDF (great for LaTeX)
    fig.savefig(f"pictures/proxy_{x_res}x{y_res}_{month}.png", dpi=300)


    ######################## ABS diffs #######################################
    ABS_dif1 = variable_iasi - variable_cams_i
    ABS_dif2 = variable_gosat - variable_cams

    norm = mcolors.TwoSlopeNorm(vmin=-10, vcenter=0, vmax=10)  # norm for colorbar

    FIGSIZE_zonal = set_size(width_fraction=1.3, height_fraction=0.32)

    fig, axs = plt.subplots(1, 2, figsize=FIGSIZE_zonal, subplot_kw={'projection': ccrs.PlateCarree()},
                            constrained_layout=True)

    def style_ax(ax, pos):
        ax.coastlines()
        # ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.LAND, facecolor='gray')
        ax.add_feature(cfeature.OCEAN, facecolor='lightgray')
        gl = ax.gridlines(draw_labels=True, linestyle='--', color='black')
        if pos == 'left':
            gl.top_labels = gl.right_labels = False
        else:
            gl.top_labels = gl.right_labels = gl.left_labels = False
        gl.xlabel_style = {'size': 8, 'color': 'black'}
        gl.ylabel_style = {'size': 8, 'color': 'black'}

    # ---- Row 1: IASI GOSAT ----
    style_ax(axs[0], pos='left')
    pcm00 = axs[0].pcolormesh(lon_bnds, lat_bnds, ABS_dif1, transform=ccrs.PlateCarree(), cmap='seismic'
                              , norm=norm)
    axs[0].set_title(rf'IASI - CAMS (IASI AVK) {x_res}$^\circ$$\times${y_res}$^\circ$ grid 2020-{month:02d}')

    style_ax(axs[1], pos='right')
    pcm01 = axs[1].pcolormesh(lon_bnds, lat_bnds, ABS_dif2, transform=ccrs.PlateCarree(), cmap='seismic'
                              , norm=norm)
    axs[1].set_title(rf'GOSAT-2 - CAMS {x_res}$^\circ$$\times${y_res}$^\circ$ grid 2020-{month:02d}')

    # shared colorbar
    cbar = fig.colorbar(pcm00, ax=axs[:], orientation='horizontal', fraction=0.02, pad=0.04, aspect=50)
    cbar.set_label(r'$\mathrm{xN}_2\mathrm{O}$ in ppb', fontsize=cbar_fontsize)

    fig.set_constrained_layout_pads(
        w_pad=1 / 72,  # padding around axes (inches, ~1pt)
        h_pad=1 / 72,
        wspace=0.02,  # space between columns (fraction of axes)
        hspace=0.00  # space between rows
    )

    base = f"pictures/abs_difs_{x_res}x{y_res}_{month}"
    # fig.savefig(base + ".pdf")  # vector PDF (great for LaTeX)
    fig.savefig(f"pictures/abs_difs_{x_res}x{y_res}_{month}.png", dpi=300)


def IASI_GOSAT_plots(month):

    mpl.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.size": 10.0,  # <- from FONTSIZE in your log
        "axes.labelsize": 9.0,
        "axes.titlesize": 10.95,  # actually does matter
        "legend.fontsize": 8,  # ~\footnotesize
        "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage{lmodern}",
        # add packages you actually use in labels, e.g.:
        # r"\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{siunitx}\usepackage{mhchem}"
    })

    FIGSIZE = set_size(width_fraction=1.3, height_fraction=0.63)
    cbar_fontsize = 9

    dir_path = 'monthly_means/{}_{:02d}_{}x{}_th{}_qf{}.nc'

    #month = 3
    x_res = 5
    y_res = 5
    th_iasi = 10
    th_gosat = 3
    th_cams = 0
    qf = 3
    var = 'mean_tot'
    met = 0
    vmin = 280
    vmax = 340

    ds_iasi = xr.open_dataset(dir_path.format('iasi', month, x_res, y_res, th_iasi, qf))

    gosat_path = 'monthly_means/{}_{:02d}_{}x{}_th{}.nc'
    ds_gosat = xr.open_dataset(gosat_path.format('gosat', month, x_res, y_res, th_gosat))

    lon_iasi = ds_iasi['lon_cen'].values
    lat_iasi = ds_iasi['lat_cen'].values
    variable_iasi = ds_iasi[var].values
    variable_gosat = ds_gosat[var].values

    iasi_zonal = zonal_mean(variable_iasi)
    gosat_zonal = zonal_mean(variable_gosat)

    lon_bnds = bnds(lon_iasi)
    lat_bnds = bnds(lat_iasi)

    iasi_proxy = variable_iasi - iasi_zonal
    gosat_proxy = variable_gosat - gosat_zonal

    ######################## GOSAT IASI diffs #######################################
    abs_dif = variable_iasi - variable_gosat
    proxy_dif = iasi_proxy - gosat_proxy

    norm = mcolors.TwoSlopeNorm(vmin=-10, vcenter=0, vmax=10)  # norm for colorbar

    FIGSIZE_zonal = set_size(width_fraction=1.3, height_fraction=0.32)

    fig, axs = plt.subplots(1, 2, figsize=FIGSIZE_zonal, subplot_kw={'projection': ccrs.PlateCarree()},
                            constrained_layout=True)

    def style_ax(ax, pos):
        ax.coastlines()
        # ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.LAND, facecolor='gray')
        ax.add_feature(cfeature.OCEAN, facecolor='lightgray')
        gl = ax.gridlines(draw_labels=True, linestyle='--', color='black')
        if pos == 'left':
            gl.top_labels = gl.right_labels = False
        else:
            gl.top_labels = gl.right_labels = gl.left_labels = False
        gl.xlabel_style = {'size': 8, 'color': 'black'}
        gl.ylabel_style = {'size': 8, 'color': 'black'}

    # ---- Row 1: IASI GOSAT ----
    style_ax(axs[0], pos='left')
    pcm00 = axs[0].pcolormesh(lon_bnds, lat_bnds, abs_dif, transform=ccrs.PlateCarree(), cmap='seismic'
                              , norm=norm)
    axs[0].set_title(rf'IASI - GOSAT-2 {x_res}$^\circ$$\times${y_res}$^\circ$ grid 2020-{month:02d}')

    style_ax(axs[1], pos='right')
    pcm01 = axs[1].pcolormesh(lon_bnds, lat_bnds, proxy_dif, transform=ccrs.PlateCarree(), cmap='seismic'
                              , norm=norm)
    axs[1].set_title(rf'IASI - GOSAT-2 proxys {x_res}$^\circ$$\times${y_res}$^\circ$ grid 2020-{month:02d}')

    # shared colorbar
    cbar = fig.colorbar(pcm00, ax=axs[:], orientation='horizontal', fraction=0.02, pad=0.04, aspect=50)
    cbar.set_label(r'$\mathrm{xN}_2\mathrm{O}$ in ppb', fontsize=cbar_fontsize)

    fig.set_constrained_layout_pads(
        w_pad=1 / 72,  # padding around axes (inches, ~1pt)
        h_pad=1 / 72,
        wspace=0.02,  # space between columns (fraction of axes)
        hspace=0.00  # space between rows
    )

    base = f"pictures/iasi_gosat_diffs_{x_res}x{y_res}_{month}"
    # fig.savefig(base + ".pdf")  # vector PDF (great for LaTeX)
    fig.savefig(f"pictures/iasi_gosat_diffs_{x_res}x{y_res}_{month}.png", dpi=300)


def set_size(textwidth_pt=466.62503, textheight_pt=674.33032, width_fraction=1.0, height_fraction=0.4):
    """
    Convert LaTeX text sizes (pt) to a Matplotlib figure size (inches).

    * \textwidth=466.62503pt * \textheight=674.33032pt

    Parameters
    ----------
    textwidth_pt : float   # \textwidth  in pt (from LaTeX log)
    textheight_pt: float   # \textheight in pt (from LaTeX log)
    width_fraction : float # e.g. 1.0 for full width, 0.5 for half width
    height_fraction: float # e.g. 0.4 for 40% of \textheight

    Returns
    -------
    (width_in, height_in) in inches
    TEXTWIDTH_PT  = 418.25555
    TEXTHEIGHT_PT = 595.80026
    """
    inches_per_pt = 1/72.27
    width_in  = textwidth_pt  * inches_per_pt * width_fraction
    height_in = textheight_pt * inches_per_pt * height_fraction
    return (width_in, height_in)


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


def grad_difs(prof):
    return prof[:-1] - prof[1:]


def gradiant(zonal):
    max_idx = int(np.nanargmax(zonal))
    south = zonal[:max_idx+1][::-1]
    north = zonal[max_idx:]

    south_grads = grad_difs(south)
    north_grads = grad_difs(north)

    mean_south = np.nanmean(south_grads)
    mean_north = np.nanmean(north_grads)
    mean_global_grad = np.nanmean([mean_north, mean_south])

    return mean_global_grad


def lat_fmt(y, _pos):
    if y < 0:
        return f"{int(abs(y))}°S"
    elif y > 0:
        return f"{int(y)}°N"
    else:
        return "0°"

