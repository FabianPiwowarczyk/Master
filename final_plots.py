import numpy as np
import os
from os.path import join
from IASI.reading import read_all_iasi


def date_tuple(data_path):
    months_paths = [f.path for f in os.scandir(data_path) if f.is_dir()]
    months_paths.sort()
    months = [m.replace(data_path, '').replace('/', '') for m in months_paths]

    date_array = np.array([0, 0])
    for m in months:
        days_paths = [f.path for f in os.scandir(join(data_path, m)) if f.is_dir()]
        days_paths.sort()
        days = [d.replace(join(data_path, m), '').replace('/', '') for d in days_paths]
        for day in days:
            date_array = np.vstack((date_array, np.array([m, day])))

    date_tuples = date_array[1:, :]
    return date_tuples


def is_point_in_square(lon, lat, square=None):
    """
    Check if a point (lon, lat) is inside a predefined square inside the Sahara Desert.

    :param square: target square
    :param lon: Longitude of the point
    :param lat: Latitude of the point
    :return: True if the point is inside the square, False otherwise
    """
    if square is None:
        square = {'min_lon': 10.0,
                  'max_lon': 25.0,
                  'min_lat': 15.0,
                  'max_lat': 30.0}

    return (square['min_lon'] <= lon <= square['max_lon']) and (square['min_lat'] <= lat <= square['max_lat'])


def find_one_specific_profile():
    data_path = '/misc/hypatia/data/IASI/L2_MUSICA/V3.30/MetopA/2020'
    quality_flag = 3
    date_tuples = date_tuple(data_path)
    square = {
        'min_lon': 5.5,
        'max_lon': 15.5,
        'min_lat': 47.0,
        'max_lat': 55.5
                }

    for i, date in enumerate(date_tuples):

        if int(date[0]) != 7:
            continue
        else:
            path = join(data_path, str(date[0]), str(date[1]))
            print(path)

            iasi_data = read_all_iasi(dir=path, date=date, i=i, date_tuples=date_tuples,
                                      org_path=data_path, quality_flag=quality_flag)

            for i, time in enumerate(iasi_data['readable_time']):
                lat = iasi_data['lat'][i]
                lon = iasi_data['lon'][i]

                if is_point_in_square(lon=lon, lat=lat, square=square):
                    if 9 <= time.hour < 14:
                        nan_count = len(iasi_data['n2o_lev_dry'][i, :]) - iasi_data['num_lev'][i]
                        print(time)
                        print(lon, lat)
                        print(iasi_data['n2o_lev_dry'][i, nan_count:])
                        print(iasi_data['alt_lev'][i, nan_count:])

        # 2020 - 07 - 01
        # 09: 16:49
        # lon 6.2223215
        # lat 47.194534
        # n2o [0.33287    0.3329     0.33332002 0.33408    0.33514002 0.33636
        #  0.33769    0.3389     0.33977    0.34007    0.33978    0.33897
        #  0.33777    0.33630002 0.33436    0.33020997 0.32274    0.31509
        #  0.31426    0.29402    0.20883001 0.16556    0.10083999 0.023201
        #  0.0047898  0.0015097  0.00080287]
        # alt_lev [601.    833.2  1298.   1816.   2365.   2944.   3562.   4211.   4890.
        #  5609.   6367.   7166.   7995.   8864.   9763.  10900.  11980.  13640.
        #  15940.  18290.  22080.  26200.  30670.  36330.  42360.  48540.  55590.]