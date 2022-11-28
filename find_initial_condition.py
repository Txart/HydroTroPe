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
import geopandas as gpd
import os
from tqdm import tqdm

import hydro_masters # contains the main hydro functions at the highest level of abstraction
import read_preprocess_data

# necessary to set the same solver also in csc
os.environ["FIPY_SOLVERS"] = "scipy"


blockOpt = False
N_CPU = 1

parent_directory = Path(r"C:\Users\03125327\github\fc_hydro_kalimantan_2022")
data_parent_folder = Path(r"C:\Users\03125327\Dropbox\PhD\Computation\ForestCarbon\2022 Kalimantan customer work\0. Raw Data")
fn_pointers = parent_directory.joinpath(r'file_pointers.xlsx')

# %% Prepare data
filenames_df = pd.read_excel(fn_pointers, header=2, dtype=str, engine='openpyxl')

# check that weather station locations and sourcesink data have the same names
if not hydro_masters._is_same_weather_station_names_in_sourcesink_and_coords(fn_pointers):
    raise ValueError('Weather station names in the weather station coords file and the sourcesink file must be equal. ABORTING.')

graph_fn = Path(filenames_df[filenames_df.Content ==
                'channel_network_graph_pickle'].Path.values[0])
graph = pickle.load(open((graph_fn), "rb"))

fn_sourcesink = Path(filenames_df[filenames_df.Content == 'sourcesink'].Path.values[0])
sourcesink_df = pd.read_excel(fn_sourcesink)


#%% Set up all the classes needed for the model

channel_network = ChannelNetwork(
    graph=graph, block_height_from_surface=0.0, block_coeff_k=2000.0,
    y_ini_below_DEM=1.0, Q_ini_value=0.0, channel_bottom_below_DEM=8.0,
    y_BC_below_DEM=-0.5, Q_BC=0.0, channel_width=3.5, work_without_blocks=not blockOpt)

peatland = Peatland(cn=channel_network, fn_pointers=fn_pointers)

peat_hydro_params = PeatlandHydroParameters(
    dt=1/24,  # dt in days
    dx=50,  # dx in meters, only used if structured mesh
    max_sweeps=1000, fipy_desired_residual=1e-5,
    s1=0.0, s2=0.0, t1=0, t2=0,
    use_several_weather_stations=True)

# Set up cwl computation
cwl_params = CWLHydroParameters(g=9.8,  # g in meters per second
                                dt=3600,  # dt in seconds
                                dx=100,  # dx in meters
                                ntimesteps=1,  # outer loop. Number of steps simulated
                                preissmann_a=0.6,
                                # threshold position for n_manning in metres below dem
                                porous_threshold_below_dem=3.0,
                                n1=5,  # params for n_manning
                                n2=1,
                                max_niter_newton=int(1e5), max_niter_inexact=int(50),
                                rel_tol=1e-7, abs_tol=1e-7, weight_A=5e-2, weight_Q=5e-2, ncpus=1,
                                downstream_diri_BC=False)

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

# %% Function
def find_best_initial_condition(initial_zeta_value, param_number, PARAMS, hydro, cwl_hydro, sensor_coords, sensor_measurements, output_folder_path):
    hydro.ph_params.s1 = float(PARAMS[PARAMS.number == param_number].s1)
    hydro.ph_params.s2 = float(PARAMS[PARAMS.number == param_number].s2)
    hydro.ph_params.t1 = float(PARAMS[PARAMS.number == param_number].t1)
    hydro.ph_params.t2 = float(PARAMS[PARAMS.number == param_number].t2)
    cwl_hydro.cwl_params.porous_threshold_below_dem = float(PARAMS[PARAMS.number == param_number].porous_threshold)
    cwl_hydro.cwl_params.n1 = float(PARAMS[PARAMS.number == param_number].n1)
    cwl_hydro.cwl_params.n2 = float(PARAMS[PARAMS.number == param_number].n2)

    hydro.parameterization = ExponentialBelowOneAboveStorageExpoTrans(hydro.ph_params)
    
    hydro.zeta = hydro.create_uniform_fipy_var(
        uniform_value=initial_zeta_value, var_name='zeta')

    # Initialize returned variable
    best_initial_zeta = hydro.zeta.value
    best_fitness = np.inf

    
    MEAN_P_MINUS_ET = -0.003 # m/day. 
    hydro.ph_params.use_several_weather_stations = False
    hydro.set_sourcesink_variable(value=MEAN_P_MINUS_ET)

    N_DAYS = 100
    day = 0
    # If True, start day0 with a small timestep to smooth things
    needs_smaller_timestep = True
    NORMAL_TIMESTEP = 24  # Hourly
    SMALLER_TIMESTEP = 10000
    while day < N_DAYS:
        # Variables from current timestep for flexible solve
        hydro_old = copy.deepcopy(hydro)
        cwl_hydro_old = copy.deepcopy(cwl_hydro)
        if not needs_smaller_timestep:
            internal_timesteps = NORMAL_TIMESTEP
        elif needs_smaller_timestep:
            internal_timesteps = SMALLER_TIMESTEP

        hydro.ph_params.dt = 1/internal_timesteps  # dt in days
        hydro.cn_params.dt = 86400/internal_timesteps  # dt in seconds
        # Add pan ET
        # we need to insert this here, otherwise the panET is added every day.
        hydro.set_sourcesink_variable(value=MEAN_P_MINUS_ET)
        hydro.sourcesink = hydro.sourcesink - \
            hydro.compute_pan_ET_from_ponding_water(hydro.zeta)

        try:
            solution_function = hydro_masters.simulate_one_timestep_simple_two_step

            for i in tqdm(range(internal_timesteps)):  # internal timestep
                hydro, cwl_hydro = solution_function(hydro, cwl_hydro)

        except Exception as e:
            if internal_timesteps == NORMAL_TIMESTEP:
                print(f'Exception in INITIAL RASTER computation of param number {param_number} at day {day}: ', e,
                      ' I will retry with a smaller timestep')
                needs_smaller_timestep = True
                hydro = hydro_old
                cwl_hydro = cwl_hydro_old
                continue

            elif internal_timesteps == SMALLER_TIMESTEP:
                print(
                    f"Another exception caught in INITIAL RASTER computation with smaller timestep: {e}. ABORTING")
                break

        else:  # try was successful
            # Compute fitness
            zeta_values_at_sensor_locations = hydro.sample_zeta_values_at_locations(coords=sensor_coords)
            current_fitness = np.linalg.norm(
                sensor_measurements - zeta_values_at_sensor_locations)

            # If fitness is improved, update initial zeta. And pickle it for future use
            is_fitness_improved = current_fitness < best_fitness
            if is_fitness_improved:
                print(
                    f"zeta at day {day} is better than the previous best, with a fitness difference of {best_fitness - current_fitness} points")
                best_initial_zeta = hydro.zeta.value
                best_fitness = current_fitness

                fn_pickle = output_folder_path.joinpath('best_initial_zeta.p')

                pickle.dump(best_initial_zeta, open(fn_pickle, 'wb'))

            # go to next day
            day = day + 1
            needs_smaller_timestep = False

    return best_initial_zeta

# %% Params
hydro.ph_params.dt = 1/24  # dt in days
hydro.cn_params.dt = 3600  # dt in seconds

# Read params
params_fn = Path.joinpath(parent_directory, '2d_calibration_parameters.xlsx')
PARAMS = pd.read_excel(params_fn, engine='openpyxl')
N_PARAMS = N_CPU

#%% Initial condition sensor measurements and locations
ini_dipwell_fn = Path(filenames_df[filenames_df.Content ==
                'initial_dipwell_measurements'].Path.values[0])
ini_dipwell_data = pd.read_csv(ini_dipwell_fn)

ini_dipwell_coords = np.column_stack((ini_dipwell_data['x'].values,ini_dipwell_data['y'].values))
ini_dipwell_WTD_meters = ini_dipwell_data['WTD(m)'].values

if ini_dipwell_coords.shape[0] != len(ini_dipwell_WTD_meters):
    raise ValueError('Different number of dipwell locations and WTD measurements. ABORTING')

# %% Run Windows
if platform.system() == 'Windows':
    hydro.verbose = True
    param_number = 3
    initial_zeta_value = 1.0 # Meters. Negative means below the surface. Zero means complete saturation everywhere. 
    output_folder_path = parent_directory.joinpath(r'initial_condition')
    find_best_initial_condition(initial_zeta_value,
                                param_number, PARAMS, hydro, cwl_hydro,
                                sensor_coords=ini_dipwell_coords,
                                sensor_measurements=ini_dipwell_WTD_meters,
                                output_folder_path=parent_directory.joinpath('initial_condition'))



# %%
