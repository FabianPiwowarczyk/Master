#!/usr/bin/env python3


def main():
    from CAMS.mean_cols import plot_mean_columns_multi
    plot_mean_columns_multi(months=[6, 7, 8], coord=[(0, 60), (-65, -50)], idx=1)

    # 1) Minimum over Canada: 60–65°N, 90–110°W, AMJ
    plot_mean_columns_multi(months=[4, 5, 6], coord=[(-110, -90), (60, 65)], idx=2)

    # 2) Atlantic minumum pops up and vanishes: 40–45°N, 30–40°W, MAM
    plot_mean_columns_multi(months=[3, 4, 5], coord=[(0, 30), (30, 40)], idx=3)

    # 3) W. Pacific jet exit: 30–40°N, 130–160°E, SON
    plot_mean_columns_multi(months=[9, 10, 11], coord=[(130, 160), (30, 40)], idx=4)

    # 4) Tibetan Plateau / W. China: 30–35°N, 80–95°E, JAS
    plot_mean_columns_multi(months=[9, 10, 11], coord=[(80, 95), (30, 35)], idx=5)

    # 7) SE Pacific stratocumulus: 15–25°S, 90–75°W, SON
    plot_mean_columns_multi(months=[7, 8, 9], coord=[(-90, -75), (-25, -10)], idx=6)

if __name__ == "__main__":
    main()