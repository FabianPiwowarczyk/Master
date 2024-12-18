import matplotlib.pyplot as plt
import numpy as np
from dask.dot import label


def change_prior(tc, dry_col, pre_lvl, avk_dict, row, lat, lon, alt_lev):

    avk = calc_avk(avk_dict=avk_dict)

    fig1 = plt.figure()

    # for i in range(avk_dict['num_lev']):
    #     plt.plot(avk[:, i], pre_lvl[::-1])

    fig2 = plt.figure()
    cols = list(range(0, avk_dict['num_lev'] + 1, 5))

    lev = list(range(0, avk_dict['num_lev']))
    for i in cols:
        plt.plot(avk[:, i], alt_lev[::-1], label=f'{alt_lev[::-1][i]}')
        plt.legend()
        plt.title(f'avk_row{row}, lat:{lat}, lon:{lon}')

    outpath = f'avk_pics/avk_row{row}.png'
    plt.savefig(outpath)

    #breakpoint()


    xn2o = 330.e-3  # ppm
    n2o_prof = np.array([0.31999999,
                         0.31999999,
                         0.31999999,
                         0.31999999,
                         0.31999999,
                         0.31999999,
                         0.31999999,
                         0.31999999,
                         0.31999999,
                         0.31999999,
                         0.31999999,
                         0.31999999,
                         0.31999999,
                         0.31967527,
                         0.31827208,
                         0.31372729,
                         0.30633333,
                         0.29544127,
                         0.26624754,
                         0.17457867])

    n2o_prof = xn2o * n2o_prof / np.mean(n2o_prof)  # GOSAT prior

    # pressure layers IASI
    pre_lay = np.zeros_like(pre_lvl[:-1])
    for i in range(len(pre_lay)):
        pre_lay[i] = (pre_lvl[i] + pre_lvl[i+1]) / 2

    cum_sum_par = np.cumsum(dry_col)  # Particle count for interpolation

    # Particles per layer for GOSAT minus the last layer bc this would give the last iasi level by design
    n_lay = [i * np.sum(dry_col) / len(n2o_prof) for i in range(1, len(n2o_prof))]

    gosat_pre_lev = np.zeros(len(n2o_prof)+1)  # array for gosat pre levels
    gosat_pre_lev[0] = pre_lvl[0]
    gosat_pre_lev[-1] = pre_lvl[-1]  # the first and last level are identical by design

    # interpolation
    gosat_pre_lev[1:-1] = np.interp(n_lay, cum_sum_par, pre_lay)

    gosat_pre_lay = np.zeros_like(gosat_pre_lev[:-1])
    for i in range(len(gosat_pre_lay)):
        gosat_pre_lay[i] = (gosat_pre_lev[i] + gosat_pre_lev[i+1]) / 2

    plt.plot(n_lay, gosat_pre_lay[:-1])
    plt.scatter(cum_sum_par, pre_lay)


def calc_avk(avk_dict):

    num_lev = avk_dict['num_lev']
    rank = int(avk_dict['avk_rank'])
    eig_val = avk_dict['avk_val'][:rank]
    lvec = avk_dict['avk_lvec'][:rank, :num_lev].T
    rvec = avk_dict['avk_rvec'][:rank, :num_lev].T

    rxr_val = np.diag(eig_val)
    avk = np.dot(np.dot(lvec, rxr_val), rvec.T)

    return avk