import matplotlib.pyplot as plt
import numpy as np


def change_prior(tc, dry_col, pre_lvl):

    pre_lay = np.zeros_like(pre_lvl[:-1])

    for i in range(len(pre_lay)):
        pre_lay[i] = (pre_lvl[i] + pre_lvl[i+1]) / 2

    cum_sum_par = np.cumsum(dry_col)

    plt.scatter(pre_lay, cum_sum_par)
    plt.show()

    breakpoint()