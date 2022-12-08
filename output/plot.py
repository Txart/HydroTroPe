#%%
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

import read_preprocess_data

#%% defs
parent_folder = Path(r"C:\Users\03125327\Dropbox\PhD\Computation\ForestCarbon\2022 Kalimantan customer work\0. Raw Data")
parent_directory = Path(r"C:\Users\03125327\github\fc_hydro_kalimantan_2022")
fn_pointers = parent_directory.joinpath(r'file_pointers.xlsx')
filenames_df = pd.read_excel(
    fn_pointers, header=2, dtype=str, engine='openpyxl')
fn_sourcesink = Path(filenames_df[filenames_df.Content == 'sourcesink'].Path.values[0])
fn_revised_dipwells = parent_folder.joinpath(r"08. Dipwell and Weather Measurements\Revised_dipwell\Dipwell Measurements (May-November 2022).csv")

# %% Create sourcesink Excel from weather data
# Weather 
p, ET = read_preprocess_data.get_precip_and_ET(parent_folder=parent_folder)
sourcesink_df = (p.drop('Date', axis='columns') - ET) / 1000 # mm -> m
sourcesink_df.to_excel(fn_sourcesink, index=False)

# Dipwell
dipwell_data = read_preprocess_data.read_dipwell_data(parent_folder)
revised_dipwell_data = pd.read_csv(fn_revised_dipwells)

# Remove dipwells close to canals
dipwell_names_close_to_canals = [name for name in revised_dipwell_data.columns if ('-1'  in name  or '-2'  in name)]
rev_dip_data_no_canals = revised_dipwell_data.drop(labels=dipwell_names_close_to_canals, axis='columns')
# %% Plot weather
plt.figure()
p.drop('Date', axis=1).plot(title='P')
plt.savefig('output/precipitation.png')

plt.figure()
ET.plot(title='ET')
plt.savefig('output/ET.png')


#%% Plot dipwell data
plt.figure()
plt.title('Dipwell measured WTD')
plt.plot(dipwell_data.drop('Date', axis='columns'), marker='.')
plt.hlines(y=0, xmin=0, xmax=len(WTD_data), colors='black', linestyles='dashed')
plt.ylabel('WT (cm)')
plt.xlabel('Time (days)')


plt.figure()
plt.title('Revised dipwells (outliers trimmed)')
plt.plot(revised_dipwell_data.drop('Date', axis='columns'), marker='.')
plt.hlines(y=0, xmin=0, xmax=len(revised_dipwell_data), colors='black', linestyles='dashed')
plt.ylabel('WT (cm)')
plt.xlabel('Time (days)')


plt.figure()
plt.title('Revised dipwells far from canals (-1 and -2 dipwells removed)')
plt.plot(rev_dip_data_no_canals.drop('Date', axis='columns'), marker='.')
plt.hlines(y=0, xmin=0, xmax=len(rev_dip_data_no_canals), colors='black', linestyles='dashed')
plt.ylabel('WT (cm)')
plt.xlabel('Time (days)')

# %%

# %%

# %%
