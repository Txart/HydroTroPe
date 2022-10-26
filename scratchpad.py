#%% Imports
import pandas as pd
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

# %% Set same starting date in rainfall and dipwell
first_rainfall_date = sourcesink_df.Date.min()
WTD_data_from_June = WTD_data[WTD_data['Date'] >= first_rainfall_date].reset_index(drop=True)

# %%

