import numpy as np
import ECMWF_Read_V3_5 as ecmwf
from .reading import read_all_iasi
from tqdm import tqdm
from .chng_prior import change_prior

from .conv_constants import *


def total_column(gas_lay, dry_col):
    """
    compute total column (average dry air mole fraction)
    input:
        gas_lay : gas mixing ration in layers
        dry_col : dry column
    """

    tc = np.sum(gas_lay * dry_col) / np.sum(dry_col)  # last dim = altitude

    return tc


def dry_column(pre_lev, h2o_dry, gra_lay):
    """
    calculate dry column
    input:
        pre_lev : pressure on levels
        h2o_dry : dry air mole fraction of h2o in ppm
        gra_lay : gravitational acceleration at layers
    """

    # calculate delta pressure in Pa divided by 10^4 in order to calculate number of particles over 1cm^2 (not 1m^2)
    pre_del = (pre_lev[0:-1] - pre_lev[1:])
    pre_del /= 100.
    # convert dry mole fraction of h2o to wet mole fraction (not in ppm)
    hdry = h2o_dry * 1e-6  # ppm -> 1
    h2o_wet = hdry / (1. + hdry)
    # calculate molar mass of wet air

    mwet = mdry * (1. - h2o_wet) + mh2o * h2o_wet
    # calculate number of wet air particles per layer
    wet_col = pre_del * avogadro / mwet / gra_lay
    # calculate number of dry air particles per layer
    dry_col = wet_col * (1. - h2o_wet)

    return dry_col


def h2o_to_hum(h2o_lay):
    """
    convert dry air mole fraction (n_h2o / n_dryair) in ppm to specific humidity (mass h2o / mass wet air)
    """
    hum_lay = ((h2o_lay * (mh2o / mdry) * 1.e-6)**-1 + 1)**-1
    return hum_lay


def grav_acc(gra_geo, gmh_lay):
    """
    output:
        gra_lay : gravitational acceleration at layers
    """
    gra3d = np.zeros_like(gmh_lay)
    gra3d[:] = gra_geo
    gra_lay = gra3d + fac * gmh_lay

    return gra_lay


def geometric_height(phi_lay, gra_geo):
    """
    output:
        gmh_lay : geometric height at layers
    """

    gra3d = np.zeros_like(phi_lay)
    gra3d[:] = gra_geo

    gmh_lay = -gra3d / fac - np.sqrt(np.clip((gra3d ** 2.) / (fac ** 2.) + (2. / fac) * phi_lay, 0, None))

    return gmh_lay


def grav_geoid(lat):
    """
    input:
        lat : latitude
    output:
        gra_geo : gravity constant at geoid (NN)
    """

    gra_geo = 9.780327 * (1. + 0.0053024 * (np.sin(lat / 180. * np.pi)) ** 2. - 0.0000058 * (np.sin(2. * lat / 180. * np.pi)) ** 2.)

    return gra_geo


def geopot_layers(phi_lev, tem_lay, hum_lay, pre_lev):
    """
    output : geopotential at layers
    """

    phi_lay = np.zeros_like(tem_lay)
    phi_lay[0:-1] = phi_lev[0:-2] + gasc_air * tem_lay[0:-1] * (
                1. + (gasc_h2o / gasc_air - 1.) * hum_lay[0:-1]) * (
                                      1. - pre_lev[1:-1] / (pre_lev[0:-2] - pre_lev[1:-1]) * (
                                          np.log(pre_lev[0:-2]) - np.log(pre_lev[1:-1])))

    phi_lay[-1] = phi_lev[-2] + np.log(2.) * gasc_air * tem_lay[-1] * (
                1. + (gasc_h2o / gasc_air - 1.) * hum_lay[-1])

    return phi_lay


def geopot_levels(phi_sur, tem_lay, hum_lay, pre_lev):
    """
    output:
        phi_lev : geopotential at levels
    """

    phi_lev = np.zeros_like(pre_lev)
    phi_lev[0] = phi_sur
    phi_lev[1:-1] = gasc_air * tem_lay[0:-1] * (1. + (gasc_h2o / gasc_air - 1.) * hum_lay[0:-1]) * (np.log(pre_lev[0:-2] / pre_lev[1:-1]))
    phi_lev = phi_lev.cumsum()
    phi_lev[-1] = np.inf

    return phi_lev


def conv_IASI_time_to_ecmwf(IASI_time):
    # Julian date for the epoch (2000-01-01 00:00:00 UTC)
    julian_date_epoch = 2451545.0

    # Calculate the difference in days for each element in the array
    difference_in_days = IASI_time / 86400.0

    # Calculate the Julian dates
    julian_dates = julian_date_epoch + difference_in_days

    return julian_dates


def read_phi_sur(data):
    path = '/misc/hypatia/data/IUP_MetData/ECMWF_era5/'
    time = conv_IASI_time_to_ecmwf(data['time'])

    phi_sur = []

    print('Reading Surface Potential:')

    items = list(range(data['time'].shape[0]))

    for i in tqdm(items):
        lon = data['lon'][i]
        lat = data['lat'][i]

        dict = ecmwf.ecmwf_read(time[i], lon=lon, lat=lat, phi_sur=True, path=path)

        phi_sur.append(dict['phi_sur'])

    return np.array(phi_sur)


def lev2lay(x_lev):
    return (x_lev[1:] + x_lev[:-1]) / 2


def altitude_lev2lay(alt_lev, pre_lev):
    z0 = alt_lev[:-1]
    z1 = alt_lev[1:]
    p0 = pre_lev[:-1]
    p1 = pre_lev[1:]
    p = (p0 + p1) / 2

    alt_lay = z0 + ((z1 - z0) / (np.log(p1) - np.log(p0))) * (np.log(p) - np.log(p0))

    return alt_lay


def data2total_col(path, date, i, date_tuples, org_path, quality_flag):

    iasi_data = read_all_iasi(dir=path, date=date, i=i, date_tuples=date_tuples,
                              org_path=org_path, quality_flag=quality_flag)

    tot_col = np.zeros_like(iasi_data['lon'])
    met0_tc = np.zeros_like(iasi_data['lon'])
    apr_col = np.zeros_like(iasi_data['lon'])
    apr_gosat = np.full((iasi_data['lon'].shape[0], 29), np.nan)
    iasi_avk = np.full((iasi_data['lon'].shape[0], 29, 29), np.nan)
    iasi_drycol = np.full((iasi_data['lon'].shape[0], 28), np.nan)  # 28 layers

    print('Converting data to total columns:')
    items = list(range(iasi_data['time'].shape[0]))

    for row in tqdm(items):

        col_avk_dict = {
        'avk_rank': iasi_data['avk_rank'][row],
        'avk_val': iasi_data['avk_val'][row],
        'avk_lvec': iasi_data['avk_lvec'][row],
        'avk_rvec': iasi_data['avk_rvec'][row],
        'num_lev': iasi_data['num_lev'][row]
                        }

        nan_count = len(iasi_data['n2o_lev_dry'][row, :]) - iasi_data['num_lev'][row]  # number of nan in the profile

        col_dic = {
        'n2o_lev': iasi_data['n2o_lev_dry'][row, nan_count:],
        'h2o_lev': iasi_data['h2o_lev_dry'][row, nan_count:],
        'pre_lev': iasi_data['pre_lev'][row, nan_count:],
        'alt_lev': iasi_data['alt_lev'][row, nan_count:],
        'apri': iasi_data['apri'][row, nan_count:]
                    }

        row_lat = iasi_data['lat'][row]

        col_dic['pre_lev'] /= 100  # from Pa to hPa

        n2o_lay = lev2lay(x_lev=col_dic['n2o_lev'])
        h2o_lay = lev2lay(x_lev=col_dic['h2o_lev'])

        alt_lay = altitude_lev2lay(alt_lev=col_dic['alt_lev'], pre_lev=col_dic['pre_lev'])

        gra_geo = grav_geoid(lat=row_lat)

        gra_lay = grav_acc(gra_geo=gra_geo, gmh_lay=alt_lay)

        dry_col = dry_column(pre_lev=col_dic['pre_lev'], h2o_dry=h2o_lay, gra_lay=gra_lay)

        tc = total_column(gas_lay=n2o_lay, dry_col=dry_col)

        tot_col[row] = tc

        # prior correction and avk calculation
        met0_tc[row], apr_gosat[row, nan_count:29], iasi_avk[row, nan_count:29, nan_count:29] = change_prior(dry_col,
                                    col_dic['pre_lev'], col_avk_dict, col_dic['alt_lev'], col_dic['apri'], n2o_lay)

        # additional apriori total column
        apri_lay = lev2lay(col_dic['apri'])
        tc_apri = total_column(gas_lay=apri_lay, dry_col=dry_col)

        apr_col[row] = tc_apri
        iasi_drycol[row, nan_count:28] = dry_col[:]

        if 0.2 >= tc >= 0.5:
            print('tc is out of bounds: ', tc)

    iasi_data['total_column'] = tot_col
    iasi_data['tc_cor_met0'] = met0_tc
    iasi_data['tc_apri'] = apr_col
    iasi_data['gosat_apri'] = apr_gosat
    iasi_data['avk'] = iasi_avk
    iasi_data['dry_col'] = iasi_drycol
    print('Returning total columns.')

    return iasi_data


def is_point_in_square(lon, lat):
    """
    Check if a point (lon, lat) is inside a predefined square inside the Sahara Desert.

    :param lon: Longitude of the point
    :param lat: Latitude of the point
    :return: True if the point is inside the square, False otherwise
    """
    square = {
        'min_lon': 10.0,
        'max_lon': 25.0,
        'min_lat': 15.0,
        'max_lat': 30.0
    }

    return (square['min_lon'] <= lon <= square['max_lon']) and (square['min_lat'] <= lat <= square['max_lat'])