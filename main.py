#%%
from classes.parameterizations import ExponentialBelowOneAboveStorageExpoTrans
from classes.peatland_hydrology import set_up_peatland_hydrology
from classes.peat_hydro_params import PeatlandHydroParameters
from classes.peatland import Peatland
from classes.channel_hydrology import set_up_channel_hydrology, CWLHydroParameters
from classes.channel_network import ChannelNetwork

import argparse
import multiprocessing as mp
import pandas as pd
import platform
import pickle
import fipy as fp
from pathlib import Path
import numpy as np

import hydro_masters # contains the main hydro functions at the highest level of abstraction
import utilities

import os
# necessary to set the same solver also in csc
os.environ["FIPY_SOLVERS"] = "scipy"


#%% Parse cmd line arguments

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

#%% Folders and data
parent_directory = Path()
# data_parent_folder = Path("~/Dropbox/PhD/Computation/ForestCarbon/2022 Kalimantan customer work/0. Raw Data")
fn_pointers = parent_directory.joinpath(r'file_pointers.xlsx')

filenames_df = pd.read_excel(fn_pointers, header=2, dtype=str, engine='openpyxl')

# Check  output folder exists
output_folder = Path(filenames_df[filenames_df.Content ==
            'output_directory'].Path.values[0])
utilities.is_output_folder(output_folder)

# check that weather station locations and sourcesink data have the same names
if not hydro_masters._is_same_weather_station_names_in_sourcesink_and_coords(fn_pointers):
    raise ValueError('Weather station names in the weather station coords file and the sourcesink file must be equal. ABORTING.')

    
graph = utilities.read_graph(filenames_df)
sourcesink_df = utilities.read_sourcesink(filenames_df)

#%% Read params
params_fn = Path.joinpath(parent_directory, 'parameters.xlsx')
file_params_php = pd.read_excel(params_fn, engine='openpyxl', sheet_name='peat_hydro_prop')
file_params_general = pd.read_excel(params_fn, engine='openpyxl', sheet_name='general')
file_params_channel = pd.read_excel(params_fn, engine='openpyxl', sheet_name='channel')
file_params_peat = pd.read_excel(params_fn, engine='openpyxl', sheet_name='peat')

#%% Set up all the classes needed for the model

channel_network = ChannelNetwork(
    graph=graph,
    params_channel=file_params_channel,
    work_without_blocks=not blockOpt)

peatland = Peatland(cn=channel_network,
                    fn_pointers=fn_pointers)

peat_hydro_params = PeatlandHydroParameters(
                            params_peat=file_params_peat,
                            dt=1/24,  # dt in days
                            dx=50,  # dx in meters, only used if structured mesh
                            use_several_weather_stations=True # If False, sets single P - ET everywhere. If true, uses different weather stations' position.
                            )

# Set up cwl computation
cwl_params = CWLHydroParameters(
                                params_channel=file_params_channel,
                                dt=3600,  # s 
                                dx=100,  # m
                                downstream_diri_BC=False 
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
                                  use_scaled_pde=False,
                                  zeta_diri_bc=float(file_params_peat['zeta_diri_BC']),
                                  force_ponding_storage_equal_one=False,
                                  peatland=peatland, peat_hydro_params=peat_hydro_params,
                                  parameterization=parameterization,
                                  channel_network=channel_network, cwl_params=cwl_params)

#%% Params
hydro.ph_params.dt = 1/24  # dt in days
hydro.cn_params.dt = 3600  # dt in seconds
NDAYS = int(file_params_general['NDAYS'])

N_PARAMS = N_CPU

#%% Initial WTD
# Read from pickle
initial_zeta_origin = 'pickle' # pickle, tiff or constant

initial_zeta = utilities.get_initial_zeta(origin=initial_zeta_origin, filenames_df=filenames_df)

#%% Run hydro
if N_PARAMS > 1:
    # It only needs small adjustsments, see the 1 parameter case and extend
    raise Warning('Multiprocessing is not 100 % tested.')
    hydro.verbose = True
    param_numbers = [1]
    multiprocessing_arguments = [(param_number, file_params_php, hydro, cwl_hydro, net_daily_source,
                                  output_folder) for param_number in param_numbers]
    with mp.Pool(processes=N_CPU) as pool:
        pool.starmap(produce_family_of_rasters, multiprocessing_arguments)

elif N_PARAMS == 1:
    hydro.verbose = True
    param_numbers = [1] 

    for param_number in param_numbers:
        # Set initial zeta
        hydro.zeta = fp.CellVariable(
            name='zeta', mesh=hydro.mesh, value=initial_zeta, hasOld=True)

        hydro_masters.set_hydrological_params(hydro,
                                              params_hydro=file_params_php,
                                              param_number=param_number)

        hydro_masters.produce_family_of_rasters(param_number, hydro, cwl_hydro, NDAYS, sourcesink_df,
                output_folder)
