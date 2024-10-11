import netCDF4 as nc
import IASI
import numpy as np
from os.path import join
import os


# file_path = "data/IASIA_MUSICA_030300_L2_AllTargetProducts_20200101000557_68496.nc"

# Open the dataset
# dataset = nc.Dataset(file_path, mode='r')

# Print a list of all variable names
# print("Variable names:")
# for var_name in dataset.variables:
#     print(var_name)


# def printvar(var_name):
#     # Access the variable
#     variable = dataset.variables[var_name]
#
#     # Print the variable's name
#     print("Variable:", var_name)
#
#     # Print the variable's description (if available)
#     if hasattr(variable, 'description'):
#         print("Description:", variable.description)
#
#     # Print the variable's data
#     print("Data:")
#     if var_name in ['musica_ghg']:
#         print(variable[:, 0][0])
#     else:
#         print(variable[0].filled(np.nan))  # Print all the values of the variable
#
#     if hasattr(variable, 'dimensions'):
#         print("Dimensions:", variable.dimensions)


#printvar('musica_ghg')
#printvar('musica_altitude_levels')
#printvar('musica_pressure_levels')
#printvar('musica_fit_quality_flag')

#variable = dataset.variables['musica_altitude_levels'][:]
#print(variable[0])

#variable = dataset.variables['musica_ghg'][:, 0]
#print(variable[0])

# count = 0
# for i in range(variable.shape[0]):
#     if float(0) not in variable[i, :]:
#         count += 1
#
# print(variable.shape[0], count)

# for i in range(5, 10):
#     print(variable[i, :])

#variable = dataset.variables['musica_fit_quality_flag']

# from convert_func import convert
# convert()

import gridding
# gridding.grid_func.monthly_mean(5, 5, 0,
#                                 'gosat_data/IUP-GHG-L2-N2O-GOSAT2-FOCAL-2020',
#                                 'xn2o', 'longitude', 'latitude',
#                                 'gosat')
#
# gridding.grid_func.monthly_mean(5, 5, 0,
#                                 'finished_data/2020_',
#                                 'total_column', 'lon', 'lat',
#                                 'iasi')

gridding.plot_lv3.plot_lv3_data(5, 5, 0,
                                'monthly_means/gosat_{}_{}x{}_th{}.nc', 'gosat', 280, 340)

gridding.plot_lv3.plot_lv3_data(5, 5, 0,
                                'monthly_means/iasi_{}_{}x{}_th{}.nc', 'iasi', 280, 340)

gridding.plot_dif(5, 5, 0,
                  'monthly_means/{}_{}_{}x{}_th{}.nc', -26, 13)

# gridding.grid_func.monthly_mean(1, 1, 0,
#                                 'finished_data/2020_',
#                                 'total_column', 'lon', 'lat',
#                                 'iasi')

gridding.plot_lv3.plot_lv3_data(1, 1, 0,
                                'monthly_means/iasi_{}_{}x{}_th{}.nc', 'iasi', 280, 340)

gridding.combined_plot(5, 5, 0, 'monthly_means/{}_{}_{}x{}_th{}.nc',
                       'mean_tot', 280, 340)

