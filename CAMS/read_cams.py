import glob2
from os.path import join
import netCDF4 as nc
import numpy as np


def read_cams_data():
    dir = '/misc/ghgcci3/data/CAMS_N2O_v21r1/'
    paths = glob2.glob(join(f'{dir}', 'cams73_v21r1_n2o_conc_surface_inst_2020*.nc'))
    paths.sort()

    dataset = nc.Dataset(paths[0], mode='r')

    N2O = dataset.variables['N2O'][:, :, :, :]

    A = dataset.variables['ap'][:]
    B = dataset.variables['bp'][:]
    Psurf = dataset.variables['Psurf'][:, :, :]

    # Expand A and B to allow broadcasting
    A_expanded = A[None, :, None, None]  # Shape (1, 80, 1, 1)
    B_expanded = B[None, :, None, None]  # Shape (1, 80, 1, 1)

    # Expand Psurf to match
    Psurf_expanded = Psurf[:, None, :, :]  # Shape (248, 1, 143, 145)

    P = A_expanded + B_expanded * Psurf_expanded  # Shape (248, 80, 143, 145)

    del_P = P[:, :-1, :, :] - P[:, 1:, :, :]

    tc = np.sum(del_P * N2O, axis=1) / Psurf

    print(tc.shape)
    print(tc[0, 0, 0])

