import netCDF4 as nc
import IASI
import numpy as np
from os.path import join
import os


file_path = "data/IASIA_MUSICA_030300_L2_AllTargetProducts_20200101000557_68496.nc"

# Open the dataset
dataset = nc.Dataset(file_path, mode='r')

# Print a list of all variable names
#print("Variable names:")
#for var_name in dataset.variables:
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
        #print(variable[0])
        print(variable[0][0][:][:])
        #print(variable[0])


    if hasattr(variable, 'dimensions'):
        print("Dimensions:", variable.dimensions)

#printvar('musica_fit_quality_flag')
#printvar('musica_nol')
#printvar('musica_ghg_apriori')
#printvar('musica_ghg_avk_rank')
#printvar('musica_ghg_avk_val')
#printvar('musica_ghg_avk_lvec')
#printvar('musica_ghg_avk_rvec')

#data = IASI.read_iasi('data/IASIA_MUSICA_030300_L2_AllTargetProducts_20200101000557_68496.nc')


# rank = int(dataset.variables['musica_ghg_avk_rank'][0])
# eig_val = dataset.variables['musica_ghg_avk_val'][0][:rank]
#
# lvec = dataset.variables['musica_ghg_avk_lvec'][0][0][:rank][:].T
# rvec = dataset.variables['musica_ghg_avk_rvec'][0][0][:rank][:].T
#
# rxr_val = np.diag(eig_val)
#
# avk = np.dot(np.dot(lvec, rxr_val), rvec.T)

#print(avk)


from convert_func import convert
convert()

import gridding
#gridding.replot_all()


# xn2o = 330.e-3  # ppm
# n2o_prof = np.array([0.31999999,
#                      0.31999999,
#                      0.31999999,
#                      0.31999999,
#                      0.31999999,
#                      0.31999999,
#                      0.31999999,
#                      0.31999999,
#                      0.31999999,
#                      0.31999999,
#                      0.31999999,
#                      0.31999999,
#                      0.31999999,
#                      0.31967527,
#                      0.31827208,
#                      0.31372729,
#                      0.30633333,
#                      0.29544127,
#                      0.26624754,
#                      0.17457867])
#
# n2o_prof = xn2o * n2o_prof / np.mean(n2o_prof)
#
# print(n2o_prof)

