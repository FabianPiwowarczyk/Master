from .plot_lv3 import plot_lv3_data, plot_dif, combined_plot
from .grid_func import monthly_mean


def replot_all():
    # monthly_mean(5, 5, 0,
    #              'gosat_data/IUP-GHG-L2-N2O-GOSAT2-FOCAL-2020',
    #              'xn2o', 'longitude', 'latitude', 'gosat')

    plot_lv3_data(5, 5, 0,
                  'monthly_means/gosat_{}_{}x{}_th{}.nc', 'gosat', 280, 340)

    plot_lv3_data(5, 5, 0,
                  'monthly_means/iasi_{}_{}x{}_th{}.nc', 'iasi', 280, 340)

    plot_dif(5, 5, 0,
             'monthly_means/{}_{}_{}x{}_th{}.nc', -26, 13)

    monthly_mean(1, 1, 0,
                 'finished_data/2020_',
                 'total_column', 'lon', 'lat', 'iasi')

    plot_lv3_data(1, 1, 0,
                  'monthly_means/iasi_{}_{}x{}_th{}.nc', 'iasi', 280, 340)

    combined_plot(5, 5, 0, 'monthly_means/{}_{}_{}x{}_th{}.nc',
                  'mean_tot', 280, 340)


def regrid_iasi():
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