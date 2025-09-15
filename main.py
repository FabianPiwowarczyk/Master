import netCDF4 as nc
import numpy as np
from os.path import join
import os
import xarray as xr


file_path = "data/IASIA_MUSICA_030300_L2_AllTargetProducts_20200101000557_68496.nc"

# Open the dataset
dataset = nc.Dataset(file_path, mode='r')

# Print a list of all variable names
# print("Variable names:")
# for var_name in dataset.variables:
#    print(var_name)


def printvar(var_name):
    # Access the variable
    variable = dataset.variables[var_name]

    # Print the variable's name
    print("Variable:", var_name)

    # Print the variable's description (if available)
    if hasattr(variable, 'description'):
        print("Description:", variable.description)

    # Print the variable's data
    print("Data:")
    if var_name in ['musica_ghg']:
        print(len(variable[:, 0][0]))
    else:
        print(variable[:])
        #print(variable[0][0][:][:])

    if hasattr(variable, 'dimensions'):
        print("Dimensions:", variable.dimensions)

from IASI.mean_column import mean_apris

mean_apris()

#printvar('musica_fit_quality_flag')
#printvar('musica_ghg')

from CAMS.mean_cols import plot_mean_columns_multi

plot_mean_columns_multi(months=[6, 7, 8], coord=[(0, 60), (-65, -50)], idx=1)

# # 1) Minimum over Canada: 60–65°N, 90–110°W, AMJ
# plot_mean_columns_multi(months=[4, 5, 6], coord=[(-110, -90), (60, 65)], idx=2)
#
# # 2) Atlantic minumum pops up and vanishes: 40–45°N, 30–40°W, MAM
# plot_mean_columns_multi(months=[3, 4, 5], coord=[(0, 30), (30, 40)], idx=3)
#
# # 3) W. Pacific jet exit: 30–40°N, 130–160°E, SON
# plot_mean_columns_multi(months=[9, 10, 11], coord=[(130, 160), (30, 40)], idx=4)
#
# # 4) Tibetan Plateau / W. China: 30–35°N, 80–95°E, JAS
# plot_mean_columns_multi(months=[9, 10, 11], coord=[(80, 95), (30, 35)], idx=5)
#
# # 7) SE Pacific stratocumulus: 15–25°S, 90–75°W, SON
# plot_mean_columns_multi(months=[7, 8, 9], coord=[(-90, -75), (-25, -10)], idx=6)

#from convert_func import convert
#convert()

import gridding
#gridding.regrid_iasi()
#gridding.replot_all()


#from concurrent.futures import ProcessPoolExecutor
#from functools import partial

#inputs = [(1, 6), (7, 12)]  # Two different inputs

#import CAMS
#CAMS.cams_monthly_means(x_res=5, y_res=5, months=[1, 12])

# Use partial to fix a and b
# partial_func = partial(CAMS.cams_monthly_means, 5, 5)
#
# with ProcessPoolExecutor(max_workers=2) as executor:
#     results = list(executor.map(partial_func, inputs))

#from final_plots import find_one_specific_profile
#find_one_specific_profile()

#import netCDF4 as nc
#import glob

# Find all matching files
# files = glob.glob("monthly_means/cams_iasi_*_5x5_th0.nc")
#
# for path in files:
#     print(f"Processing {path}...")
#     with nc.Dataset(path, mode="r+") as ds:  # r+ means read/write without overwriting file
#         ds.variables['mean_tot'][:] *= 1000
#
# import xarray as xr

# gosat_path = 'monthly_means/{}_{}_{}x{}_th{}.nc'
# ds_cams_iasi = xr.open_dataset(gosat_path.format('cams_iasi', '01', 5, 5, 0))
#
# variable_cams_iasi = ds_cams_iasi['mean_tot'].values




