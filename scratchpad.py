#%% Imports
import numpy as np
import pandas as pd
import geopandas as gpd
from pathlib import Path
import matplotlib.pyplot as plt
from pathlib import Path

import read_preprocess_data
#%% defs
parent_folder = Path(r"C:\Users\03125327\Dropbox\PhD\Computation\ForestCarbon\2022 Kalimantan customer work\0. Raw Data")
parent_directory = Path(r"C:\Users\03125327\github\fc_hydro_kalimantan_2022")
fn_pointers = parent_directory.joinpath(r'file_pointers.xlsx')
filenames_df = pd.read_excel(
    fn_pointers, header=2, dtype=str, engine='openpyxl')
fn_sourcesink = Path(filenames_df[filenames_df.Content == 'sourcesink'].Path.values[0])

# %% Create sourcesink Excel from weather data
# Weather 
p, ET = read_preprocess_data.get_precip_and_ET(parent_folder=parent_folder)
sourcesink_df = (p.drop('Date', axis='columns') - ET) / 1000 # mm -> m
sourcesink_df.to_excel(fn_sourcesink, index=False)

# Dipwell
WTD_data = read_preprocess_data.read_dipwell_data(parent_folder)

#%% Remove dipwells close to canals and create file. Useful for calibration
fn_dipwell_data = Path(filenames_df[filenames_df.Content == 'dipwell_measurements'].Path.values[0])
dipwell_data = pd.read_csv(fn_dipwell_data)
dipwell_names_close_to_canals = [name for name in dipwell_data.columns if ('-1'  in name  or '-2'  in name)]
dip_data_no_canals = dipwell_data.drop(labels=dipwell_names_close_to_canals, axis='columns')
OUTPUT_FN = parent_directory.joinpath('data/revised_dipwell_data_from_first_rainfall_record_without_canal_sensors.csv')
dip_data_no_canals.to_csv(OUTPUT_FN, index=False)

#%% Plot dipwell data
plt.figure()
plt.title('Dipwell measured WTD')
plt.plot(WTD_data.drop('Date', axis='columns'), marker='.')
plt.hlines(y=0, xmin=0, xmax=len(WTD_data), colors='black', linestyles='dashed')
plt.ylabel('WT (cm)')
plt.xlabel('Time (days)')

# %% Set same starting date in rainfall and dipwell
first_rainfall_date = p.Date.min()
WTD_data_from_June = WTD_data[WTD_data['Date'] >= first_rainfall_date].reset_index(drop=True)
# Save WTD data to file
fn_june_dipwell_data = parent_directory.joinpath('data/dipwell_data_from_first_rainfall_record.csv')
WTD_data_from_June.to_csv(fn_june_dipwell_data, index=False)

#%% Plot data to choose initial condition
# 3rd day (23rd of June) has a very big amount of data (only 20 records missing)
# I choose this as initial condition
plt.figure()
plt.plot(WTD_data_from_June.isnull().sum(axis=1))

#%% Initial day dipwell values and locations
dipwell_fn = Path(filenames_df[filenames_df.Content ==
                'sensor_locations'].Path.values[0])
ini_day = 3
# Values
ini_day_dipwell_WTD = WTD_data_from_June.drop('Date', axis='columns').loc[ini_day]
ini_day_dipwell_names = list(WTD_data_from_June.loc[ini_day][WTD_data_from_June.loc[ini_day].notnull()].index)
ini_day_dipwell_names.remove('Date')
# Locations
dipwell_fn = Path(filenames_df[filenames_df.Content == 'sensor_locations'].Path.values[0])
dipwell_location = gpd.read_file(dipwell_fn)
ini_day_dipwell_locations = dipwell_location[dipwell_location['Dipwell ID'].isin(ini_day_dipwell_names)]

# Joint  dataframe of locations and WTD measurements
ini_day_dipwell_coords_and_WTD = pd.DataFrame(index=ini_day_dipwell_locations.index,
                                              columns=['ID', 'x', 'y', 'WTD(m)'])
for dipwell_name in ini_day_dipwell_names:
    dipwell_index = ini_day_dipwell_locations[ini_day_dipwell_locations['Dipwell ID'] == dipwell_name].index.values[0]
    dipwell_xcoord = ini_day_dipwell_locations[ini_day_dipwell_locations['Dipwell ID'] == dipwell_name].geometry.x.values[0]
    dipwell_ycoord = ini_day_dipwell_locations[ini_day_dipwell_locations['Dipwell ID'] == dipwell_name].geometry.y.values[0]
    dipwell_WTD_meters = ini_day_dipwell_WTD[dipwell_name]/100 # meters
    ini_day_dipwell_coords_and_WTD.loc[dipwell_index, 'ID'] = dipwell_name
    ini_day_dipwell_coords_and_WTD.loc[dipwell_index, 'x'] = dipwell_xcoord
    ini_day_dipwell_coords_and_WTD.loc[dipwell_index, 'y'] = dipwell_ycoord
    ini_day_dipwell_coords_and_WTD.loc[dipwell_index, 'WTD(m)'] = dipwell_WTD_meters

# Write info to files
fn_ini_day_dipwell_coords_and_WTD = parent_directory.joinpath('initial_condition/initial_day_dipwell_coords_and_measurements.csv')
ini_day_dipwell_coords_and_WTD.to_csv(fn_ini_day_dipwell_coords_and_WTD, index=False)
