#%% imports

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