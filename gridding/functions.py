from .plot_lv3 import plot_lv3_data, combined_plot, zonal_plot, dif_plot
from .grid_func import monthly_mean
from .final_plots import correction_example, tot_mean_example, seasonal_plots, proxy_plots, IASI_GOSAT_plots


def replot_all():

    th_gosat = 3
    th_iasi = 10
    th_cams = 0

    # plot_lv3_data(5, 5, th_gosat,
    #               'monthly_means/{}_{}_{}x{}_th{}.nc', 'gosat',
    #               'mean_tot', 280, 340)
    #
    # plot_lv3_data(5, 5, th_iasi,
    #               'monthly_means/{}_{}_{}x{}_th{}_qf{}.nc', 'iasi',
    #               'mean_tot', 280, 340, 3)
    # plot_lv3_data(5, 5, th_iasi,
    #               'monthly_means/{}_met{}_{}_{}x{}_th{}_qf{}.nc', 'iasi',
    #               'mean_tot', 280, 340, 3, 0)
    #
    # combined_plot(5, 5, th_iasi, th_gosat, 'monthly_means/{}_{}_{}x{}_th{}_qf{}.nc',
    #               'mean_tot', 3, vmin=-30, vmax=30,
    #               met_path='monthly_means/{}_met{}_{}_{}x{}_th{}_qf{}.nc', met=0)
    #
    # zonal_plot(5, 5, th_iasi, th_gosat, th_cams, 'monthly_means/{}_{}_{}x{}_th{}_qf{}.nc',
    #               'mean_tot', 3, vmin=None, vmax=None,
    #               met_path='monthly_means/{}_met{}_{}_{}x{}_th{}_qf{}.nc', met=0)
    #
    # # CAMS Plots
    # plot_lv3_data(5, 5, th_cams,
    #               'monthly_means/{}_{}_{}x{}_th{}.nc', 'cams',
    #               'mean_tot', 280, 340)
    #
    # plot_lv3_data(5, 5, th_cams,
    #               'monthly_means/{}_{}_{}x{}_th{}.nc', 'cams_iasi',
    #               'mean_tot', 280, 340)
    #
    # dif_plot(5, 5, th_iasi, th_gosat, th_cams, 'monthly_means/{}_{}_{}x{}_th{}_qf{}.nc',
    #               'mean_tot', 3, vmin=-30, vmax=30,
    #               met_path='monthly_means/{}_met{}_{}_{}x{}_th{}_qf{}.nc', met=0)
    #
    #
    for m in range(1, 13):
        correction_example(m)

    for m in range(1, 13):
       tot_mean_example(m)
    seasonal_plots()

    for m in range(1, 13):
       proxy_plots(m)

    for m in range(1, 13):
        IASI_GOSAT_plots(m)



def regrid_iasi():

    print('Grid gosat data th=0')
    monthly_mean(x_res=5, y_res=5, th=0,
                 dir_path='gosat_data/IUP-GHG-L2-N2O-GOSAT2-FOCAL-2020',
                 tot_var='xn2o', lon_var='longitude', lat_var='latitude', sat='gosat')

    print('Grid gosat data th=3')
    monthly_mean(x_res=5, y_res=5, th=3,
                 dir_path='gosat_data/IUP-GHG-L2-N2O-GOSAT2-FOCAL-2020',
                 tot_var='xn2o', lon_var='longitude', lat_var='latitude', sat='gosat')

    # print('Regrid Iasi in 1x1 th=0')
    # print('Normal tot_col')
    path = f'/misc/ghgcci7/fabian/iasi_data/2020_'
    # monthly_mean(x_res=1, y_res=1, th=0,
    #              dir_path=path,
    #              tot_var='total_column', lon_var='lon', lat_var='lat', sat='iasi', qf=3)
    #
    # print('Method 0 correction')
    # monthly_mean(x_res=1, y_res=1, th=0,
    #              dir_path=path,
    #              tot_var='tc_cor_met0', lon_var='lon', lat_var='lat', sat='iasi', qf=3)

    print('Regrid Iasi in 5x5 th=0')
    print('Normal tot_col')
    monthly_mean(x_res=5, y_res=5, th=0,
                 dir_path=path,
                 tot_var='total_column', lon_var='lon', lat_var='lat', sat='iasi', qf=3)

    print('Method 0 correction')
    monthly_mean(x_res=5, y_res=5, th=0,
                 dir_path=path,
                 tot_var='tc_cor_met0', lon_var='lon', lat_var='lat', sat='iasi', qf=3)

    print('Regrid Iasi in 5x5 th=10')
    print('Normal tot_col')
    monthly_mean(x_res=5, y_res=5, th=10,
                 dir_path=path,
                 tot_var='total_column', lon_var='lon', lat_var='lat', sat='iasi', qf=3)

    print('Method 0 correction')
    monthly_mean(x_res=5, y_res=5, th=10,
                 dir_path=path,
                 tot_var='tc_cor_met0', lon_var='lon', lat_var='lat', sat='iasi', qf=3)