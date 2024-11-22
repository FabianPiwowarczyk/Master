import numpy as np
import ECMWF_Read_V3_5 as ecmwf
from .reading import read_all_iasi
from tqdm import tqdm

from .conv_constants import *


def total_column(gas_lay, dry_col, row):
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


def data2total_col(path, date, i, date_tuples, org_path):

    iasi_data = read_all_iasi(dir=path, date=date, i=i, date_tuples=date_tuples, org_path=org_path)

    # print(iasi_data['n2o_lev_dry'][4559, :])
    # print(iasi_data['h2o_lev_dry'][4559, :])
    # print(iasi_data['pre_lev'][4559, :])
    # print(iasi_data['tem_lev'][4559, :])

    # print(iasi_data['lon'][3788], iasi_data['lat'][3788])

    # all_phi_sur = read_phi_sur(data=iasi_data)

    # all_phi_sur = np.loadtxt('surface_pot.txt', delimiter=' ')

    # np.savetxt('surface_pot.txt', X=all_phi_sur, delimiter=' ')

    tot_col = np.zeros_like(iasi_data['lon'])
    nan_flags = np.zeros_like(iasi_data['lon'])

    tc_nan_count = 0

    print('Converting data to total columns:')
    items = list(range(iasi_data['time'].shape[0]))

    for row in tqdm(items):

        col_dic = {
        'n2o_lev': iasi_data['n2o_lev_dry'][row, :],
        'h2o_lev': iasi_data['h2o_lev_dry'][row, :],
        'pre_lev': iasi_data['pre_lev'][row, :],
        'alt_lev': iasi_data['alt_lev'][row, :]
                    }

        if (nan_count := np.isnan(col_dic['n2o_lev']).sum()) > 0:
            nan_flags[row] = nan_count
            if nan_count == len(col_dic['n2o_lev']):
                raise ValueError(f"In row {row} are only nan values.")

            for key in col_dic.keys():
                if np.isnan(col_dic[key]).sum() == nan_count:
                    col_dic[key] = col_dic[key][nan_count:]
                else:
                    raise ValueError(f"In row {row} the nan value counts are not equal!")

                if np.isnan(col_dic[key]).any():
                    raise ValueError(f"In row {row} are still nan values left after removal.")

        row_lat = iasi_data['lat'][row]

        # phi_sur = all_phi_sur[row]

        col_dic['pre_lev'] /= 100  # from Pa to hPa

        n2o_lay = lev2lay(x_lev=col_dic['n2o_lev'])
        # tem_lay = lev2lay(x_lev=col_dic['tem_lev'])
        h2o_lay = lev2lay(x_lev=col_dic['h2o_lev'])

        alt_lay = altitude_lev2lay(alt_lev=col_dic['alt_lev'], pre_lev=col_dic['pre_lev'])

        # hum_lay = h2o_to_hum(h2o_lay=h2o_lay)
        #
        # phi_lev = geopot_levels(phi_sur=phi_sur, tem_lay=tem_lay, hum_lay=hum_lay, pre_lev=col_dic['pre_lev'])
        #
        # phi_lay = geopot_layers(phi_lev=phi_lev, tem_lay=tem_lay, hum_lay=hum_lay, pre_lev=col_dic['pre_lev'])

        gra_geo = grav_geoid(lat=row_lat)

        # gmh_lay = geometric_height(phi_lay=phi_lay, gra_geo=gra_geo)

        gra_lay = grav_acc(gra_geo=gra_geo, gmh_lay=alt_lay)

        # if float(0) in gra_lay:
        #     zero_count = len(gra_lay) - np.count_nonzero(gra_lay)
        #
        #     if zero_count == len(gra_lay):
        #         print(iasi_data['observation_id'][row])
        #         break
        #         fucking_hell_count += 1
        #         continue
        #
        #     col_dic['pre_lev'] = col_dic['pre_lev'][:-zero_count]
        #     h2o_lay = h2o_lay[:-zero_count]
        #     gra_lay = gra_lay[:-zero_count]
        #     n2o_lay = n2o_lay[:-zero_count]

        dry_col = dry_column(pre_lev=col_dic['pre_lev'], h2o_dry=h2o_lay, gra_lay=gra_lay)

        tc = total_column(gas_lay=n2o_lay, dry_col=dry_col, row=row)

        from .chng_prior import change_prior
        change_prior(tc, dry_col)

        if np.isnan(tc):
            tc_nan_count += 1
        else:
            tot_col[row] = tc

        if 0.2 >= tc >= 0.5:
            print('tc is out of bounds: ', tc)

    iasi_data['total_column'] = tot_col
    iasi_data['nan_flg'] = nan_flags
    print('Returning total columns with nan count: ', tc_nan_count)

    return iasi_data