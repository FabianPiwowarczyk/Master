#!/usr/bin/env python3


def main():
    from CAMS.mean_cols import plot_mean_columns_multi

    # 1) Rocky mountains: 35–40°N, 110–115°W, JAS
    plot_mean_columns_multi(months=[7, 8, 9], coord=[(-115, -110), (35, 40)], idx=1)

    # 2) low over Siberia: 50–55°N, 75–80°E, OND
    plot_mean_columns_multi(months=[10, 11, 12], coord=[(75, 80), (50, 55)], idx=2)

    # 3) anomalies between Africa and south America: 20–25°S, 30–35°W, AMJ
    plot_mean_columns_multi(months=[4, 5, 6], coord=[(-35, -30), (-25, 20)], idx=3)


if __name__ == "__main__":
    main()