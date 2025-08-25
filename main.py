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
        print(variable[:, 0][0])
    else:
        print(variable[:])
        #print(variable[0][0][:][:])

    if hasattr(variable, 'dimensions'):
        print("Dimensions:", variable.dimensions)

#printvar('musica_fit_quality_flag')


#from convert_func import convert
#convert()

import gridding
#gridding.regrid_iasi()
gridding.replot_all()

#from gridding.plot_lv3 import seasonal_plot
#seasonal_plot()

from concurrent.futures import ProcessPoolExecutor
from functools import partial

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

import netCDF4 as nc
import glob

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
#
# print(variable_cams_iasi)



