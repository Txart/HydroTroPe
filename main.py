# %%
from classes.parameterizations import ExponentialBelowOneAboveStorageExpoTrans
from classes.peatland_hydrology import PeatlandHydroParameters, set_up_peatland_hydrology
from classes.peatland import Peatland
from classes.channel_hydrology import set_up_channel_hydrology, CWLHydroParameters
from classes.channel_network import ChannelNetwork
import preprocess_data
import argparse
import multiprocessing as mp
from scipy.spatial import distance
import sys
import pandas as pd
import platform
import pickle
import fipy as fp
import copy
import networkx as nx
from pathlib import Path
import numpy as np

import hydro_masters # contains the main hydro functions at the highest level of abstraction
import utilities

import os
# necessary to set the same solver also in csc
os.environ["FIPY_SOLVERS"] = "scipy"


# %% Parse cmd line arguments

parser = argparse.ArgumentParser(description='Run hydrology')

parser.add_argument('--yes-blocks', help='Use status quo blocks. This is the default behaviour.',
                    dest='blockOpt', action='store_true')
parser.add_argument('--no-blocks', help='Do not use any block',
                    dest='blockOpt', action='store_false')

parser.add_argument('--ncpu', required=True,
                    help='(int) Number of processors', type=int)

args = parser.parse_args()

parser.set_defaults(blockOpt=True)
blockOpt = args.blockOpt
N_CPU = args.ncpu

parent_directory = Path(r"C:\Users\03125327\github\hydrotrope")
data_parent_folder = Path(r"C:\Users\03125327\Dropbox\PhD\Computation\ForestCarbon\2022 Kalimantan customer work\0. Raw Data")
fn_pointers = parent_directory.joinpath(r'file_pointers.xlsx')

if N_CPU != 1:
    raise ValueError(
        'Multiprocessing not impletmeented in Windows. n_cpus in must be equal to 1')

# %% Prepare data
filenames_df = pd.read_excel(fn_pointers, header=2, dtype=str, engine='openpyxl')

# check that weather station locations and sourcesink data have the same names
if not hydro_masters._is_same_weather_station_names_in_sourcesink_and_coords(fn_pointers):
    raise ValueError('Weather station names in the weather station coords file and the sourcesink file must be equal. ABORTING.')

def read_graph(filenames_df):
    graph_fn = Path(filenames_df[filenames_df.Content ==
                'channel_network_graph_pickle'].Path.values[0])
    return pickle.load(open((graph_fn), "rb"))

def read_sourcesink(filenames_df):
    fn_sourcesink = Path(filenames_df[filenames_df.Content == 'sourcesink'].Path.values[0])
    return pd.read_excel(fn_sourcesink)
    
graph = read_graph(filenames_df)
sourcesink_df = read_sourcesink(filenames_df)

#%% Set up all the classes needed for the model

channel_network = ChannelNetwork(
    graph=graph, block_height_from_surface=0.0, block_coeff_k=2000.0,
    y_ini_below_DEM=-0.0, Q_ini_value=0.0, channel_bottom_below_DEM=8.0,
    y_BC_below_DEM=-0.0, Q_BC=0.0, channel_width=3.5, work_without_blocks=not blockOpt)

peatland = Peatland(cn=channel_network, fn_pointers=fn_pointers)

peat_hydro_params = PeatlandHydroParameters(
    dt=1/24,  # dt in days
    dx=50,  # dx in meters, only used if structured mesh
    max_sweeps=1000, # max number of iterations for convergence of FiPy's numerical solution
    fipy_desired_residual=1e-5, # Residual for FiPy
    s1=0.0, s2=0.0, t1=0, t2=0, # Peat hydraulic properties. Storage coefficient and transmissivity.
    use_several_weather_stations=True # If False, sets single P - ET everywhere. If true, uses different weather stations' position.
    )

# Set up cwl computation
cwl_params = CWLHydroParameters(dt=3600,  # s 
                                dx=100,  # m
                                ntimesteps=1,  # outer loop. If =1, it simulates dt seconds into the future. If =2, 2*dt, etc. So it is usually =1.
                                preissmann_a=0.6, # Parameter in the Preissmann method
                                porous_threshold_below_dem=3.0, # threshold position for n_manning in metres below dem. See shape of n_Manning()
                                n1=5,  # params for n_manning
                                n2=1, # params for n_manning
                                max_niter_newton=int(1e5), # Maximum number of iterations for the Newton method
                                max_niter_inexact=int(50), # Max number of iters for the calculation of the Jacobian in the inexact (aka accelerated) Newton method
                                rel_tol=1e-7, abs_tol=1e-7, # Tolerances for convergence
                                weight_A=5e-2, weight_Q=5e-2, # Scale the solution of the Newton method to update the solution variables
                                ncpus=1,
                                downstream_diri_BC=False # Downstream BC for the channel reaches. If False, no flux Neumann are imposed.
                                )

cwl_hydro = set_up_channel_hydrology(model_type='diff-wave-implicit-inexact',
                                     cwl_params=cwl_params,
                                     cn=channel_network)

# If you change this, change also other occurrences below!!
# parameterization = ExponentialBelowOneAboveStorageWithDepth(peat_hydro_params)
# parameterization = ExponentialBelowOneAboveStorage(peat_hydro_params)
parameterization = ExponentialBelowOneAboveStorageExpoTrans(peat_hydro_params)

hydro = set_up_peatland_hydrology(mesh_fn=Path(filenames_df[filenames_df.Content == 'mesh'].Path.values[0]),
                                  model_coupling='darcy',
                                  use_scaled_pde=False, zeta_diri_bc=-0.2,
                                  force_ponding_storage_equal_one=False,
                                  peatland=peatland, peat_hydro_params=peat_hydro_params,
                                  parameterization=parameterization,
                                  channel_network=channel_network, cwl_params=cwl_params)



# %% Params
hydro.ph_params.dt = 1/24  # dt in days
hydro.cn_params.dt = 3600  # dt in seconds
NDAYS = 364

# Read params
params_fn = Path.joinpath(parent_directory, '2d_calibration_parameters.xlsx')
PARAMS = pd.read_excel(params_fn, engine='openpyxl')
N_PARAMS = N_CPU

#%% Initial WTD
# Read from pickle
initial_zeta_from_pickle = True
if initial_zeta_from_pickle:
    initial_zeta_pickle_fn = Path(
        filenames_df[filenames_df.Content == 'initial_zeta_pickle'].Path.values[0])
    initial_zeta = pickle.load(open(initial_zeta_pickle_fn, 'rb'))
else:
    initial_zeta_fn = r"C:\Users\03125327\github\fc_hydro_kalimantan_2022\output\params_number_1\no_blocks\first96days/zeta_after_96_DAYS.tif"
    mesh_centroids_coords = np.column_stack(hydro.mesh.cellCenters.value)
    initial_zeta = utilities.sample_raster_from_coords(
            raster_filename=initial_zeta_fn,
            coords=mesh_centroids_coords)
    initial_zeta = np.nan_to_num(initial_zeta, -0.2) # if mesh larger than raster, fill with value

# %% Run multiprocessing csc
if platform.system() == 'Linux':

    if N_PARAMS > 1:
        # It only needs small adjustsments, see the 1 parameter case and extend
        raise NotImplementedError('Multiprocessing not implemented fully')
        hydro.verbose = True
        param_numbers = [1]
        multiprocessing_arguments = [(param_number, PARAMS, hydro, cwl_hydro, net_daily_source,
                                      parent_directory) for param_number in param_numbers]
        with mp.Pool(processes=N_CPU) as pool:
            pool.starmap(produce_family_of_rasters, multiprocessing_arguments)

    elif N_PARAMS == 1:
        hydro.verbose = True
        param_numbers = [1] 
        # arguments = [(param_number, PARAMS, hydro, cwl_hydro, net_daily_source,
        #               parent_directory) for param_number in param_numbers]
        # for args in arguments:
        #     produce_family_of_rasters(*args)

        for param_number in param_numbers:
            # Set initial zeta
            hydro.zeta = fp.CellVariable(
                name='zeta', mesh=hydro.mesh, value=initial_zeta, hasOld=True)

            hydro_masters.set_hydrological_params(hydro, cwl_hydro, PARAMS, param_number)

            hydro_masters.produce_family_of_rasters(param_number, hydro, cwl_hydro, NDAYS, sourcesink_df,
                    parent_directory)


# %% Run Windows
if platform.system() == 'Windows':
    hydro.verbose = True
    N_PARAMS = 1
    param_numbers = [1]

    for param_number in param_numbers:
        # Set initial zeta
        hydro.zeta = fp.CellVariable(
            name='zeta', mesh=hydro.mesh, value=initial_zeta, hasOld=True)

        hydro_masters.set_hydrological_params(hydro, cwl_hydro, PARAMS, param_number)

        hydro_masters.produce_family_of_rasters(param_number, hydro, cwl_hydro, NDAYS, sourcesink_df,
                  parent_directory)

# %%
