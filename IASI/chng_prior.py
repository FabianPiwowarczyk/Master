import matplotlib.pyplot as plt
import numpy as np
from dask.dot import label


def change_prior(tc, dry_col, pre_lvl, avk_dict, row, lat, lon, alt_lev, old_prior):

    avk = calc_avk(avk_dict=avk_dict)

    # fig1 = plt.figure()

    # for i in range(avk_dict['num_lev']):
    #     plt.plot(avk[:, i], pre_lvl[::-1])

    # fig2 = plt.figure()
    # cols = list(range(0, avk_dict['num_lev'] + 1, 5))
    #
    # lev = list(range(0, avk_dict['num_lev']))
    # for i in cols:
    #     plt.plot(avk[:, i], alt_lev[::-1], label=f'{alt_lev[::-1][i]}')
    #     plt.legend()
    #     plt.title(f'avk_row{row}, lat:{lat}, lon:{lon}')

    # outpath = f'avk_pics/avk_row{row}.png'
    # plt.savefig(outpath)




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

    # plt.plot(n_lay, gosat_pre_lay[:-1])
    # plt.scatter(cum_sum_par, pre_lay)

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

    #print(old_prior)
    #print(new_prior)
    #print(avk)

    #print(new_prior - old_prior)
    #print(np.log(new_prior) - np.log(old_prior))
    #print(np.exp(np.log(new_prior) - np.log(old_prior)))
    #print(np.matmul((np.eye(avk.shape[0]) - avk), (new_prior[::-1] - old_prior[::-1]).T))
    #print(np.matmul((np.eye(avk.shape[0]) - avk), (new_prior[::-1] - old_prior[::-1])))


    tc_cor = np.exp(np.matmul((np.eye(avk.shape[0])-avk), (np.log(new_prior[::-1]) - np.log(old_prior[::-1])).T))
    #tc_cor_lev = np.matmul((1 - avk), (new_prior[::-1] - old_prior[::-1]))
    #print(tc_cor)
    #print((1 - avk))

    #tc_cor_lay = lev2lay(tc_cor_lev)
    #print(total_column(tc_cor_lay, dry_col))

    # plotting priors
    fig, axes = plt.subplots(2, 3, figsize=(10, 8), sharex=False, sharey=True)
    y = [i for i in range(len(new_prior))]

    axes[0, 0].plot(old_prior, y)
    axes[0, 0].set_title("IASI prior")
    axes[0, 0].set_xlabel("ppm")
    axes[0, 0].set_ylabel("Height level")
    axes[0, 0].set_xlim(0, 0.35)

    axes[0, 1].plot(new_prior, y)
    axes[0, 1].set_title("GOSAT prior")
    axes[0, 1].set_xlabel("ppm")
    axes[0, 1].set_xlim(0, 0.35)

    axes[0, 2].plot(new_prior - old_prior, y)
    axes[0, 2].set_title("Diff prior GOSAT - IASI")
    axes[0, 2].set_xlabel("ppm")

    axes[1, 0].plot(np.log(old_prior), y)
    axes[1, 0].set_title("ln IASI prior")
    axes[1, 0].set_xlabel("ln ppm")
    axes[1, 0].set_ylabel("Height level")
    axes[1, 0].set_xlim(-6.5, -1)

    axes[1, 1].plot(np.log(new_prior), y)
    axes[1, 1].set_title("ln GOSAT prior")
    axes[1, 1].set_xlabel("ln ppm")
    axes[1, 1].set_ylabel("Height level")
    axes[1, 1].set_xlim(-6.5, -1)

    axes[1, 2].plot(np.log(new_prior) - np.log(old_prior), y)
    axes[1, 2].set_title("Diff prior ln GOSAT - ln IASI")
    axes[1, 2].set_xlabel("ln ppm")
    plt.tight_layout()

    # plotting avk
    fig, axes = plt.subplots(1, 3, figsize=(10, 5), sharex=True, sharey=True)

    y = [i for i in range(avk.shape[0])]
    x0 = np.diagonal(avk)[::-1]

    max_row = np.unravel_index(np.argmax(avk), avk.shape)[0]
    max_col = np.unravel_index(np.argmax(avk), avk.shape)[1]

    axes[0].plot(x0, y)
    axes[0].set_title("AVK diag")
    axes[0].set_xlabel("AVK")
    axes[0].set_ylabel("Height level")
    axes[0].legend()

    for i in range(avk.shape[0]):
        x1 = avk[:, i][::-1]

        if i == max_col:
            axes[1].plot(x1, y, label=f'height level: {i}')
        else:
            axes[1].plot(x1, y)
        axes[1].set_title(f"AVK cols height levels {avk.shape[0]}")
        axes[1].set_xlabel("AVK")
        axes[1].legend()

    for i in range(avk.shape[0]):
        x2 = avk[i, :][::-1]

        if i == max_row:
            axes[2].plot(x2, y, label=f'height level: {i}')
        else:
            axes[2].plot(x2, y)
        axes[2].set_title(f"AVK rows height levels {avk.shape[0]}")
        axes[2].set_xlabel("AVK")
        axes[2].legend()

    plt.tight_layout()
    plt.show()

    # print(n2o_prof)
    # print(gosat_pre_lay)
    #
    # ypre = np.linspace(100, 100, 10)
    # x = np.interp(ypre, gosat_pre_lay[::-1], n2o_prof[::-1])

    # fig = plt.figure()
    # plt.scatter(n2o_prof, gosat_pre_lay)
    # plt.scatter(n2o_gosat_lay, pre_lay)
    # plt.show()

    breakpoint()


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