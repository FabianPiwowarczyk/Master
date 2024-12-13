import matplotlib.pyplot as plt
import numpy as np


def change_prior(tc, dry_col, pre_lvl):

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
    plt.show()

    breakpoint()