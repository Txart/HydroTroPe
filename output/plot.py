#%%
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

import read_preprocess_data
#%% defs
parent_folder = Path(r"C:\Users\03125327\Dropbox\PhD\Computation\ForestCarbon\2022 Kalimantan customer work\0. Raw Data")

# %% read  data
# Weather 
p, ET = read_preprocess_data.get_precip_and_ET(parent_folder=parent_folder)
sourcesink_df = (p.drop('Date', axis='columns') - ET) / 1000 # mm -> m
sourcesink_df['Date'] = p['Date']

# %% Plot weather
plt.figure()
p.drop('Date', axis=1).plot(title='P')
plt.savefig('output/precipitation.png')

plt.figure()
ET.plot(title='ET')
plt.savefig('output/ET.png')
# %%
