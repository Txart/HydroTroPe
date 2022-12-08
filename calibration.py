"""
Compares modelled to measured WTD
"""
#%% Imports
from pathlib import Path
import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt

import utilities
#%% Funcs

def get_dipwell_names_and_locations(dipwell_info:gpd.geodataframe)->tuple[list,np.ndarray]:
    dipwell_names = list(dipwell_info['Dipwell ID'].values)
    dipwell_locations = np.zeros(shape=(len(dipwell_names), 2))
    for i, row in dipwell_info.iterrows():
        dipwell_locations[i] = np.array([row.geometry.x, row.geometry.y])
    
    return dipwell_names, dipwell_locations


def get_one_param_modelled_WTD_filepaths(n_param:int, n_days:int, is_blocked:bool)-> list:
    blocked_foldername = 'yes_blocks' if is_blocked else 'no_blocks'
    modelled_WTD_directory = Path(f'output/params_number_{n_param}/{blocked_foldername}')

    return [Path(modelled_WTD_directory.joinpath(f'zeta_after_{day}_DAYS.tif')) for day in range(1, n_days+1)]
    

def mse_two_vectors_with_nans(v: np.ndarray, w: np.ndarray)-> float:
    # NaN values in any of the two vectors get masked away and do not contribute to the mse
    mask = np.logical_or(np.isnan(v), np.isnan(w)) # Any of the two vectors can have NaNs
    v_masked = np.ma.masked_array(v, mask)
    w_masked = np.ma.masked_array(w, mask)
    
    return np.sum((v_masked - w_masked))**2 / np.sum(~mask)



def sort_dipwell_WTD_according_to_names(dipwell_names: list, dipwell_WTD:pd.DataFrame) -> np.ndarray:
    # Sort measured WTD according to names
    dipwell_WTD_sorted = np.zeros(shape=(NDAYS, len(dipwell_names)))
    for nday in range(NDAYS):
        dipwell_WTD_sorted[nday] = np.array([dipwell_WTD.drop('Date', axis='columns').loc[nday+1][dipwell_name] for dipwell_name in dipwell_names]) 

    return dipwell_WTD_sorted


def mse_one_parameter_modelled_and_measured_WTD(n_param:int, measured_WTD: np.ndarray, dipwell_locations:np.ndarray) -> float:
    modelled_WTD_filepaths = get_one_param_modelled_WTD_filepaths(n_param, n_days=NDAYS, is_blocked=False)
    mse = 0
    for day, fpath in enumerate(modelled_WTD_filepaths):
        zeta_modelled_at_sensor_locations = utilities.sample_raster_from_coords(raster_filename=fpath, coords=dipwell_locations)
        mse += mse_two_vectors_with_nans(v=measured_WTD[day], w=zeta_modelled_at_sensor_locations)
     
    return mse


def remove_dipwells_without_measurements(dipwell_names:list, dipwell_locations: np.ndarray, dipwell_names_with_some_measurement:list):
    dipwell_names_with_no_measurements =  set(dipwell_names) - set(dipwell_names_with_some_measurement)
    
    # Delete names
    dipwell_names_new = [name for name in dipwell_names if name not in dipwell_names_with_no_measurements]

    # Delete locations preserving the order. Index items first:
    indices_to_remove = []
    for name in dipwell_names_with_no_measurements:
        indices_to_remove.append(dipwell_names.index(name))
    # Delete
    dipwell_locations_new = np.delete(dipwell_locations, indices_to_remove, axis=0)

    return dipwell_names_new, dipwell_locations_new

#%% Main script
if '__name__' == '__main__':
    import utilities
#%%
    # Globals
    parent_directory = Path(r"C:\Users\03125327\github\fc_hydro_kalimantan_2022")
    fn_pointers = parent_directory.joinpath(r'file_pointers.xlsx')
    filenames_df = pd.read_excel(fn_pointers, header=2, dtype=str, engine='openpyxl')

    NDAYS = 96 # Number of computed days
    PARAMS = [1, 4, 5]

    # Read measured WTD and dipwell locations
    dipwell_fn = Path(filenames_df[filenames_df.Content == 'dipwell_measurements'].Path.values[0])
    dipwell_WTD = pd.read_csv(dipwell_fn)
    dipwell_location_fn = Path(filenames_df[filenames_df.Content == 'sensor_locations'].Path.values[0])
    dipwell_info = gpd.read_file(dipwell_location_fn)
    dipwell_names, dipwell_locations = get_dipwell_names_and_locations(dipwell_info)
    # Only take those dipwells for which we have at least one measurement
    dipwell_names_with_some_measurement = list(dipwell_WTD.drop('Date', axis='columns').columns)
    dipwell_names, dipwell_locations = remove_dipwells_without_measurements(dipwell_names, dipwell_locations, dipwell_names_with_some_measurement)
    # Sort dipwells to compare properly with measured
    dipwell_WTD_sorted = sort_dipwell_WTD_according_to_names(dipwell_names, dipwell_WTD) # shape: [day, WTD_meas]
    dipwell_WTD_sorted_meters = dipwell_WTD_sorted * 0.01 # cm -> m

    # Compute mse for all params
    mse_per_param = []
    for param in PARAMS:
        mse = mse_one_parameter_modelled_and_measured_WTD(n_param=param, measured_WTD=dipwell_WTD_sorted_meters, dipwell_locations=dipwell_locations)
        mse_per_param.append(mse)

    print(f"CALIBRATION RESULT: the best parameter is param {np.argmin(mse_per_param) + 1}, with a MSE of {min(mse_per_param)}")
        
# %% Plots for visual identification

def get_measured_WTD_at_dipwell_locations(n_param:int, dipwell_locations:np.ndarray) -> np.ndarray:
    modelled_WTD_at_dipwells = np.zeros(shape=(NDAYS, len(dipwell_locations))) # Output var
    modelled_WTD_filepaths = get_one_param_modelled_WTD_filepaths(n_param, n_days=NDAYS, is_blocked=False)
    for day, fpath in enumerate(modelled_WTD_filepaths):
        modelled_WTD_at_dipwells[day] = utilities.sample_raster_from_coords(raster_filename=fpath, coords=dipwell_locations)
     
    return modelled_WTD_at_dipwells


#%% 
for p in PARAMS:
    modelled_WTD_at_dipwells = get_measured_WTD_at_dipwell_locations(n_param=p, dipwell_locations=dipwell_locations)
    plt.figure()
    plt.title(f'param {p} - Modelled vs Measured')
    plt.plot(modelled_WTD_at_dipwells)
    plt.plot(dipwell_WTD_sorted_meters, marker='.', color='gray', linestyle='None')
    plt.xlabel('Time (days)')
    plt.ylabel('WT (m)')
    
    plt.savefig(f'output/plots/param_{p}_calibration.png')

# %%
