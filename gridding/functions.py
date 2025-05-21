from .plot_lv3 import plot_lv3_data, combined_plot, zonal_plot, dif_plot
from .grid_func import monthly_mean


def replot_all():

    plot_lv3_data(5, 5, 0,
                  'monthly_means/{}_{}_{}x{}_th{}.nc', 'gosat',
                  'mean_tot', 280, 340)
    plot_lv3_data(1, 1, 0,
                  'monthly_means/{}_{}_{}x{}_th{}_qf{}.nc', 'iasi',
                  'mean_tot', 280, 340, 3)
    plot_lv3_data(1, 1, 0,
                  'monthly_means/{}_met{}_{}_{}x{}_th{}_qf{}.nc', 'iasi',
                  'mean_tot', 280, 340, 3, 0)

    plot_lv3_data(5, 5, 0,
                  'monthly_means/{}_{}_{}x{}_th{}_qf{}.nc', 'iasi',
                  'mean_tot', 280, 340, 3)
    plot_lv3_data(5, 5, 0,
                  'monthly_means/{}_met{}_{}_{}x{}_th{}_qf{}.nc', 'iasi',
                  'mean_tot', 280, 340, 3, 0)

    combined_plot(5, 5, 0, 'monthly_means/{}_{}_{}x{}_th{}_qf{}.nc',
                  'mean_tot', 3, vmin=-30, vmax=30,
                  met_path='monthly_means/{}_met{}_{}_{}x{}_th{}_qf{}.nc', met=0)

    zonal_plot(5, 5, 0, 'monthly_means/{}_{}_{}x{}_th{}_qf{}.nc',
                  'mean_tot', 3, vmin=None, vmax=None,
                  met_path='monthly_means/{}_met{}_{}_{}x{}_th{}_qf{}.nc', met=0)

    # CAMS Plots
    plot_lv3_data(5, 5, 0,
                  'monthly_means/{}_{}_{}x{}_th{}.nc', 'cams',
                  'mean_tot', 280, 340)

    dif_plot(5, 5, 0, 'monthly_means/{}_{}_{}x{}_th{}_qf{}.nc',
                  'mean_tot', 3, vmin=-30, vmax=30,
                  met_path='monthly_means/{}_met{}_{}_{}x{}_th{}_qf{}.nc', met=0)


def regrid_iasi():

    # monthly_mean(5, 5, 0,
    #              'gosat_data/IUP-GHG-L2-N2O-GOSAT2-FOCAL-2020',
    #              'xn2o', 'longitude', 'latitude', 'gosat')

    print('Regrid Iasi in 1x1')
    print('Normal tot_col')
    monthly_mean(1, 1, 0,
                 'finished_data/2020_',
                 'total_column', 'lon', 'lat', 'iasi', 3)

    print('Method 0 correction')
    monthly_mean(1, 1, 0,
                 'finished_data/2020_',
                 'tc_cor_met0', 'lon', 'lat', 'iasi', 3)

    print('Method 1 correction')
    monthly_mean(1, 1, 0,
                 'finished_data/2020_',
                 'tc_cor_met1', 'lon', 'lat', 'iasi', 3)

    print('Method 2 correction')
    monthly_mean(1, 1, 0,
                 'finished_data/2020_',
                 'tc_cor_met2', 'lon', 'lat', 'iasi', 3)

    print('Regrid Iasi in 5x5')
    print('Normal tot_col')
    monthly_mean(5, 5, 0,
                 'finished_data/2020_',
                 'total_column', 'lon', 'lat', 'iasi', 3)

    print('Method 0 correction')
    monthly_mean(5, 5, 0,
                 'finished_data/2020_',
                 'tc_cor_met0', 'lon', 'lat', 'iasi', 3)

    print('Method 1 correction')
    monthly_mean(5, 5, 0,
                 'finished_data/2020_',
                 'tc_cor_met1', 'lon', 'lat', 'iasi', 3)

    print('Method 2 correction')
    monthly_mean(5, 5, 0,
                 'finished_data/2020_',
                 'tc_cor_met2', 'lon', 'lat', 'iasi', 3)