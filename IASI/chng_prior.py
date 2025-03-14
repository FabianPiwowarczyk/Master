import matplotlib.pyplot as plt
import numpy as np
from dask.dot import label
from scipy.optimize import curve_fit


def change_prior(dry_col, pre_lvl, avk_dict, alt_lev, old_prior, gas_lay):

    avk = calc_avk(avk_dict=avk_dict)

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

    # cal weights for weighted average
    deltas = np.zeros((len(gosat_pre_lay), len(pre_lay)))

    for i in range(deltas.shape[0]):
        gpre_1 = gosat_pre_lev[i]
        gpre_2 = gosat_pre_lev[i+1]
        for j in range(deltas.shape[1]):
            ipre_1 = pre_lvl[j]
            ipre_2 = pre_lvl[j+1]

            if ipre_2 >= gpre_1:  # iasi layer is under gosat layer
                deltas[i, j] = 0.
            elif (ipre_1 <= gpre_1) and (ipre_2 >= gpre_2):  # iasi layer completely  inside gosat layer
                deltas[i, j] = 1.
            elif (ipre_1 > gpre_1) and (ipre_2 < gpre_1):  # cut from below
                i_del = ipre_1 - ipre_2
                pre_del = gpre_1 - ipre_2
                if ipre_2 < gpre_2:  # edge case if 1 iasi layer touches 3 gosat layers
                    overshoot = gpre_2 - ipre_2
                else:
                    overshoot = 0.
                deltas[i, j] = (pre_del - overshoot) / i_del
            elif (ipre_1 > gpre_2) and (ipre_2 < gpre_2):  # cut from top
                i_del = ipre_1 - ipre_2
                pre_del = ipre_1 - gpre_2
                if ipre_1 > gpre_1:  # edge case if 1 iasi layer touches 3 gosat layers
                    overshoot = ipre_1 - gpre_1
                else:
                    overshoot = 0.
                deltas[i, j] = (pre_del - overshoot) / i_del
            elif ipre_1 <= gpre_2:  # iasi layer over gosat layer
                deltas[i, j] = 0.

    n2o_gosat_lay = np.sum(n2o_prof[:, np.newaxis] * deltas, axis=0)  # calc weighted average

    gosat_n2o_lev = np.zeros_like(pre_lvl)

    # interpolate n2o ppm values for iasi lev. All flipped bc np.interp needs ascending values
    gosat_n2o_lev[1:-1] = np.interp(pre_lvl[1:-1][::-1], pre_lay[::-1], n2o_gosat_lay[::-1])[::-1]

    # calc first and last level n2o ppm value
    gosat_n2o_lev[0] = 2 * n2o_gosat_lay[0] - gosat_n2o_lev[1]
    gosat_n2o_lev[-1] = 2 * n2o_gosat_lay[-1] - gosat_n2o_lev[-2]

    # correcting interpolation error. Normalize to xn2o
    cor_fac = xn2o / total_column(n2o_gosat_lay, dry_col)
    new_prior = gosat_n2o_lev * cor_fac

    met0_new_prior = new_prior.copy()
    met0_cor_lev, met0_cor = calc_correction(avk, met0_new_prior, old_prior, dry_col, gas_lay)

    met1_new_prior = linear_methode(deltas, new_prior.copy(), pre_lvl)
    met1_cor_lev, met1_cor = calc_correction(avk, met1_new_prior, old_prior, dry_col, gas_lay)

    met2_new_prior = ignor_top_layer(deltas, new_prior.copy(), old_prior)
    met2_cor_lev, met2_cor = calc_correction(avk, met2_new_prior, old_prior, dry_col, gas_lay)

    # plot_prior_methodes(met0=met0_new_prior, met0_cor_lev=met0_cor_lev, met0_cor=met0_cor,
    #                     met1=met1_new_prior, met1_cor_lev=met1_cor_lev, met1_cor=met1_cor,
    #                     met2=met2_new_prior, met2_cor_lev=met2_cor_lev, met2_cor=met2_cor,
    #                     pre_lev=pre_lvl, old_prior=old_prior)

    return met0_cor, met1_cor, met2_cor


def calc_avk(avk_dict):

    num_lev = avk_dict['num_lev']
    rank = int(avk_dict['avk_rank'])
    eig_val = avk_dict['avk_val'][:rank]
    lvec = avk_dict['avk_lvec'][:rank, :num_lev].T
    rvec = avk_dict['avk_rvec'][:rank, :num_lev].T

    rxr_val = np.diag(eig_val)
    avk = np.dot(np.dot(lvec, rxr_val), rvec.T)

    return avk


def total_column(gas_lay, dry_col):
    """
    compute total column (average dry air mole fraction)
    input:
        gas_lay : gas mixing ration in layers
        dry_col : dry column
    """

    tc = np.sum(gas_lay * dry_col) / np.sum(dry_col)  # last dim = altitude

    return tc


def lev2lay(x_lev):
    return (x_lev[1:] + x_lev[:-1]) / 2


def ignor_top_layer(deltas, gosat_prior, old_prior):
    """
    swap all gosat prior data against the iasi prior data in the top layer
    """

    # plus 1 to account for the change from layers to levels
    gosat_prior[np.where(deltas[19, :] == 1)[0][:]+1] = old_prior[np.where(deltas[19, :] == 1)[0][:]+1]

    return gosat_prior


def linear_methode(deltas, gosat_prior, pre_lev):
    """
    linear fit between 0 pressure and last valid gosat point
    """
    gosat_last_idx = np.where(deltas[19, :] == 1)[0][0]
    slope = gosat_prior[gosat_last_idx] / pre_lev[gosat_last_idx]
    gosat_prior[gosat_last_idx+1:] = slope * pre_lev[gosat_last_idx+1:]

    return gosat_prior


def calc_correction(avk, new_prior, old_prior, dry_col, gas_lay):

    tc_cor_lev = np.matmul((np.eye(avk.shape[0]) - avk), (np.log(new_prior[::-1]) - np.log(old_prior[::-1])).T)
    tc_cor_lay = np.exp(lev2lay(tc_cor_lev[::-1]))

    tc_cor = total_column_correction(gas_lay, dry_col, tc_cor_lay)

    if 0.2 >= tc_cor >= 0.5:
        print('tc is out of bounds: ', tc_cor)

    return tc_cor_lev[::-1], tc_cor


def total_column_correction(gas_lay, dry_col, tc_cor_lay):
    tc = np.sum(gas_lay * tc_cor_lay * dry_col) / np.sum(dry_col)  # last dim = altitude
    return tc


def plot_prior_methodes(met0, met0_cor_lev, met0_cor,
                        met1, met1_cor_lev, met1_cor,
                        met2, met2_cor_lev, met2_cor,
                        pre_lev, old_prior):

    # plotting priors
    fig, axes = plt.subplots(3, 3, figsize=(10, 10), sharex=False, sharey=True)
    y = [i for i in range(len(pre_lev))]
    #y = pre_lev[:]

    #axes[0, 0].invert_yaxis()

    axes[0, 0].scatter(met0, y)
    axes[0, 0].set_title("GOSAT prior pure")
    axes[0, 0].set_xlabel("ppm")
    axes[0, 0].set_ylabel("Height level")
    axes[0, 0].set_xlim(0, 0.35)

    axes[0, 1].scatter(met0-old_prior, y)
    axes[0, 1].set_title("New - Old")
    axes[0, 1].set_xlabel("ppm")

    axes[0, 2].scatter(met0_cor_lev, y, label=f'{met0_cor}')
    axes[0, 2].set_title("Methode 0 cor")
    axes[0, 2].set_xlabel("factor")
    axes[0, 2].legend()

    axes[1, 0].scatter(met1, y)
    axes[1, 0].set_title("GOSAT prior linear")
    axes[1, 0].set_xlabel("ppm")
    axes[1, 0].set_ylabel("Height level")
    axes[1, 0].set_xlim(0, 0.35)

    axes[1, 1].scatter(met1 - old_prior, y)
    axes[1, 1].set_title("New - Old")
    axes[1, 1].set_xlabel("ppm")

    axes[1, 2].scatter(met1_cor_lev, y, label=f'{met1_cor}')
    axes[1, 2].set_title("Methode 1 cor")
    axes[1, 2].set_xlabel("factor")
    axes[1, 2].legend()

    axes[2, 0].scatter(met2, y)
    axes[2, 0].set_title("GOSAT prior ignore")
    axes[2, 0].set_xlabel("ppm")
    axes[2, 0].set_ylabel("Height level")
    axes[2, 0].set_xlim(0, 0.35)

    axes[2, 1].scatter(met2 - old_prior, y)
    axes[2, 1].set_title("New - Old")
    axes[2, 1].set_xlabel("ppm")

    axes[2, 2].scatter(met2_cor_lev, y, label=f'{met2_cor}')
    axes[2, 2].set_title("Methode 0 cor")
    axes[2, 2].set_xlabel("factor")
    axes[2, 2].legend()
    plt.tight_layout()
    plt.show()


def plot_avk(avk, alt_lev):
    fig, axes = plt.subplots(2, 3, figsize=(10, 10), sharex=True, sharey=False)

    y_lev = [i for i in range(avk.shape[0])]
    y_alt = alt_lev[:]

    max_row = np.unravel_index(np.argmax(avk), avk.shape)[0]
    max_col = np.unravel_index(np.argmax(avk), avk.shape)[1]

    for i in range(avk.shape[0]):
        x0_0 = avk[:, i][::-1]

        if i == max_col:
            axes[0, 0].plot(x0_0, y_lev, label=f'height level: {avk.shape[0] - 1 - i} {alt_lev[::-1][i]}')
        else:
            axes[0, 0].plot(x0_0, y_lev)
    axes[0, 0].set_title(f"AVK cols height levels {avk.shape[0]}")
    axes[0, 0].set_xlabel("AVK")
    axes[0, 0].legend()

    for i in range(avk.shape[0]):
        x0_1 = avk[:, i][::-1]

        if i == max_col:
            axes[0, 1].plot(x0_1, y_alt, label=f'height level: {avk.shape[0] - 1 - i} {alt_lev[::-1][i]}')
        else:
            axes[0, 1].plot(x0_1, y_alt)
    axes[0, 1].set_title(f"AVK cols altitude levels {avk.shape[0]}")
    axes[0, 1].set_xlabel("AVK")
    axes[0, 1].legend()

    for i in range(avk.shape[0] - 8, avk.shape[0] - 3):
        x0_2 = avk[:, i][::-1]

        axes[0, 2].plot(x0_2, y_alt, label=f'height level: {avk.shape[0] - 1 - i}')

    axes[0, 2].set_title(f"AVK cols height levels {avk.shape[0]}")
    axes[0, 2].set_xlabel("AVK")
    axes[0, 2].legend()

    for i in range(avk.shape[0]):
        x1_0 = avk[i, :][::-1]

        if i == max_row:
            axes[1, 0].plot(x1_0, y_lev, label=f'height level: {avk.shape[0] - 1 - i} {alt_lev[::-1][i]}')
        else:
            axes[1, 0].plot(x1_0, y_lev)
        axes[1, 0].set_title(f"AVK rows height levels {avk.shape[0]}")
        axes[1, 0].set_xlabel("AVK")
        axes[1, 0].legend()

    for i in range(avk.shape[0]):
        x1_1 = avk[i, :][::-1]

        if i == max_row:
            axes[1, 1].plot(x1_1, y_alt, label=f'height level: {avk.shape[0] - 1 - i} {alt_lev[::-1][i]}')
        else:
            axes[1, 1].plot(x1_1, y_alt)
        axes[1, 1].set_title(f"AVK rows height levels {avk.shape[0]}")
        axes[1, 1].set_xlabel("AVK")
        axes[1, 1].legend()

    for i in range(avk.shape[0] - 8, avk.shape[0] - 3):
        x1_2 = avk[i, :][::-1]

        axes[1, 2].plot(x1_2, y_alt, label=f'height level: {avk.shape[0] - 1 - i}')

    axes[1, 2].set_title(f"AVK cols height levels {avk.shape[0]}")
    axes[1, 2].set_xlabel("AVK")
    axes[1, 2].legend()

    plt.tight_layout()
    plt.show()