#!/usr/bin/env python3


def main():
    # from gridding.grid_func import monthly_mean
    #
    # path = f'/misc/ghgcci7/fabian/iasi_data/2020_'
    #
    # print('Monthly means IASI org.')
    # monthly_mean(x_res=5, y_res=5, th=4,
    #              dir_path=path,
    #              tot_var='total_column', lon_var='lon', lat_var='lat', sat='iasi', qf=3)
    #
    # print('Method 0 correction')
    # monthly_mean(x_res=5, y_res=5, th=4,
    #              dir_path=path,
    #              tot_var='tc_cor_met0', lon_var='lon', lat_var='lat', sat='iasi', qf=3)

    import CAMS
    CAMS.cams_monthly_means(x_res=5, y_res=5, months=[1, 12])

if __name__ == "__main__":
    main()