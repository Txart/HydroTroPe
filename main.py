#%% Imports
import pandas as pd
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt

import read_preprocess_data
#%% defs
parent_folder = Path(r"C:\Users\03125327\Dropbox\PhD\Computation\ForestCarbon\2022 Kalimantan customer work\0. Raw Data")

# %% read  data
# Weather 
p, ET = read_preprocess_data.get_precip_and_ET(parent_folder=parent_folder)
sourcesink_df = p.drop('Date', axis='columns') - ET
sourcesink_df['Date'] = p['Date']

# Dipwell
WTD_data = read_preprocess_data.read_dipwell_data(parent_folder)

# %% Set same starting date in rainfall and dipwell
first_rainfall_date = sourcesink_df.Date.min()
WTD_data_from_June = WTD_data[WTD_data['Date'] >= first_rainfall_date].reset_index(drop=True)

# %%

