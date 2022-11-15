#%% imports
import numpy as np
import pandas as pd
import fipy as fp
from pathlib import Path
import copy
from tqdm import tqdm

from classes.parameterizations import ExponentialBelowOneAboveStorageExpoTrans
#%%

def _is_same_weather_station_names_in_sourcesink_and_coords(fn_pointers:str)->bool:
    filenames_df = pd.read_excel(fn_pointers, header=2, dtype=str, engine='openpyxl')
    # read sourcesink
    fn_sourcesink = Path(filenames_df[filenames_df.Content == 'sourcesink'].Path.values[0])
    sourcesink_df = pd.read_excel(fn_sourcesink)
    # read ws coords
    fn_wscoords = Path(filenames_df[filenames_df.Content == 'weather_station_coords'].Path.values[0])
    wscoords = pd.read_excel(fn_wscoords) 

    return set(wscoords.Name.values) == set(sourcesink_df.columns)

# %% Simulation functions

def simulate_one_timestep_simple_two_step(hydro, cwl_hydro):
    # Alternative to simulate_one_timestep that only uses 1 fipy solution per timestep
    zeta_before = hydro.zeta
    # Simulate WTD
    hydro.theta = hydro.create_theta_from_zeta(hydro.zeta)

    # solution is saved into class attribute
    hydro.theta = hydro.run(mode='theta', fipy_var=hydro.theta)
    hydro.zeta = hydro.create_zeta_from_theta(hydro.theta)

    new_y_at_canals = hydro.convert_zeta_to_y_at_canals(hydro.zeta)
    y_prediction_array = hydro.cn.from_nodedict_to_nparray(new_y_at_canals)
    theta_prediction_array = hydro.parameterization.theta_from_zeta(
        zeta=y_prediction_array - hydro.cn.dem, dem=hydro.cn.dem)
    theta_previous_array = hydro.parameterization.theta_from_zeta(
        zeta=hydro.cn.y - hydro.cn.dem, dem=hydro.cn.dem)
    theta_diff = theta_prediction_array - theta_previous_array
    q = cwl_hydro.predict_q_for_next_timestep(
        theta_difference=theta_diff, seconds_per_timestep=hydro.cn_params.dt)

    # Add q to each component graph
    q_nodedict = hydro.cn.from_nparray_to_nodedict(q)
    for component_cn in cwl_hydro.component_channel_networks:
        component_cn.q = component_cn.from_nodedict_to_nparray(q_nodedict)

    # CWL solution step
    # cwl_hydro.run() takes y, q and others from each component's attribute variables
    cwl_hydro.run()  # solution stored in cwl_hydro.y_solution (or Q_solution also if Preissmann)
    # y is water height above common ref point
    y_solution_at_canals = cwl_hydro.y_solution

    # Distribute canal ponding water throughout the cell
    zeta_at_canals = cwl_hydro.convert_y_nodedict_to_zeta_at_canals(
        y_solution_at_canals)
    ponding_canal_nodes = cwl_hydro._find_zeta_in_ponding_canal_nodes(
        y_solution_at_canals)
    zeta_at_canals = hydro.distribute_canal_ponding_water_throughout_cell(
        zeta_at_canals, ponding_canal_nodes)
    y_solution_at_canals = cwl_hydro.convert_zeta_nodedict_to_y_at_canals(
        zeta_at_canals)

    # Append solution value at canals to each graph in the components
    for component_cn in cwl_hydro.component_channel_networks:
        component_cn.y = component_cn.from_nodedict_to_nparray(
            y_solution_at_canals)

    # Here I append the computed cwl to the global channel_network, both as channel_network.y,
    #  and also to each of the components in cwl_hydro.component_channel_networks
    # The components part is needed for the cwl computation, and the global graph part
    # is needed to taake into account nodes in the channel network outside the components
    zeta_at_canals = cwl_hydro.convert_y_nodedict_to_zeta_at_canals(
        y_solution_at_canals)
    mean_zeta_value = np.mean(list(zeta_at_canals.values()))
    for n in cwl_hydro.nodes_not_in_components:
        y_solution_at_canals[n] = mean_zeta_value + \
            cwl_hydro.cn.graph.nodes[n]['DEM']
    y_sol_all_canals = hydro.cn.from_nodedict_to_nparray(y_solution_at_canals)
    hydro.cn.y = y_sol_all_canals

    # Set solution of CWL simulation into fipy variable zeta for WTD simulation
    zeta_at_canals = cwl_hydro.convert_y_nodedict_to_zeta_at_canals(
        y_solution_at_canals)
    zeta_array = hydro.cn.from_nodedict_to_nparray(zeta_at_canals)
    hydro.set_canal_values_in_fipy_var(
        channel_network_var=zeta_array, fipy_var=hydro.zeta)

    return hydro, cwl_hydro

# %% main function

def run_daily_computations(hydro, cwl_hydro, net_daily_source, internal_timesteps, day):

    hydro.ph_params.dt = 1/internal_timesteps  # dt in days
    cwl_hydro.cwl_params.dt = 86400/internal_timesteps  # dt in seconds

    hydro.ph_params.use_several_weather_stations = True
    # # P-ET sink/source term at each weather station
    sources_dict = net_daily_source.loc[day].to_dict()
    hydro.set_sourcesink_variable(value=sources_dict)

    # Add pan ET  
    hydro.sourcesink = hydro.sourcesink - \
            hydro.compute_pan_ET_from_ponding_water(hydro.zeta)
    zeta_t0 = hydro.zeta.value

    solution_function = simulate_one_timestep_simple_two_step

    for hour in tqdm(range(internal_timesteps)):
        hydro, cwl_hydro = solution_function(hydro, cwl_hydro)

    zeta_t1 = hydro.zeta.value

    return hydro, cwl_hydro


def write_output_zeta_raster(zeta, hydro, full_folder_path, day):
    out_raster_fn = Path.joinpath(
        full_folder_path, f"zeta_after_{day}_DAYS.tif")
    hydro.save_fipy_var_in_raster_file(
        fipy_var=zeta, out_raster_fn=out_raster_fn, interpolation='linear')

    return None

def set_hydrological_params(hydro, cwl_hydro, PARAMS, param_number):
    hydro.ph_params.s1 = float(PARAMS[PARAMS.number == param_number].s1)
    hydro.ph_params.s2 = float(PARAMS[PARAMS.number == param_number].s2)
    hydro.ph_params.t1 = float(PARAMS[PARAMS.number == param_number].t1)
    hydro.ph_params.t2 = float(PARAMS[PARAMS.number == param_number].t2)
    cwl_hydro.cwl_params.porous_threshold_below_dem = float(PARAMS[PARAMS.number == param_number].porous_threshold)
    cwl_hydro.cwl_params.n1 = float(PARAMS[PARAMS.number == param_number].n1)
    cwl_hydro.cwl_params.n2 = float(PARAMS[PARAMS.number == param_number].n2)
    
    return None


def produce_family_of_rasters(param_number, hydro, cwl_hydro, N_DAYS,
                              net_daily_source, parent_directory):

    # hydro.parameterization = ExponentialBelowOneAboveStorage(
    #     hydro.ph_params)
    hydro.parameterization = ExponentialBelowOneAboveStorageExpoTrans(hydro.ph_params)
     
    # Outputs will go here
    output_directory = Path.joinpath(parent_directory, 'output')
    out_rasters_folder_name = f"params_number_{param_number}"
    full_folder_path = Path.joinpath(output_directory, out_rasters_folder_name)

    day = 0
    needs_smaller_timestep = False
    NORMAL_TIMESTEP = 24  # Hourly
    SMALLER_TIMESTEP = 1000

    while day < N_DAYS:
        print(f'\n computing day {day}')

        if not needs_smaller_timestep:
            internal_timesteps = NORMAL_TIMESTEP
        elif needs_smaller_timestep:
            internal_timesteps = SMALLER_TIMESTEP

        hydro_test = copy.deepcopy(hydro)
        cwl_hydro_test = copy.deepcopy(cwl_hydro)

        try:
            hydro_test, cwl_hydro_test = run_daily_computations(
                hydro_test, cwl_hydro_test, net_daily_source, internal_timesteps, day)

        except Exception as e:
            if internal_timesteps == NORMAL_TIMESTEP:
                print(f'Exception in computation of param number {param_number} at day {day}: ', e,
                      ' I will retry with a smaller timestep')
                needs_smaller_timestep = True
                continue

            elif internal_timesteps == SMALLER_TIMESTEP:
                print(
                    f"Another exception caught with smaller timestep: {e}. ABORTING")
                return 0

        else:  # try was successful
            # Update variables
            hydro = copy.deepcopy(hydro_test)
            cwl_hydro = copy.deepcopy(cwl_hydro_test)

            day += 1
            needs_smaller_timestep = False

            # write zeta to file
            if hydro.cn.work_without_blocks:
                foldername = 'no_blocks'
            elif not hydro.cn.work_without_blocks:
                foldername = 'yes_blocks'

            # if weather_stations, then the folder is yes_blocks by default

            full_write_foldername = full_folder_path.joinpath(foldername)
            print(f' writing output raster to {full_write_foldername}')

            write_output_zeta_raster(
                hydro.zeta, hydro, full_write_foldername, day)

            continue
    return 0
