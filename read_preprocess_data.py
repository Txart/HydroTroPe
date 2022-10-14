#%%
import pandas as pd
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from ET.estimate_daily_et import compute_ET
# %% read dipwell data
def _read_dipwell_csv(csv_path:Path)->pd.DataFrame:
    df = pd.read_csv(csv_path, sep=";")
    df.Date = pd.to_datetime(df.Date, dayfirst=True)

    return df

# def _order_dipwell_df_chronologically(df: pd.DataFrame)->pd.DataFrame:
#     df.sort_values(by='Date')
#     df.reset_index(drop=True, inplace=True)
#     return df.set_index(keys='Date')

def read_dipwell_data(parent_folder:Path)->pd.DataFrame:
    fn_WTD_daily = parent_folder.joinpath(r"08. Dipwell and Weather Measurements\Dipwell\WTL_Daily.csv")
    fn_WTD_notdaily = parent_folder.joinpath(r"08. Dipwell and Weather Measurements\Dipwell\WTL_Others.csv")

    WTD_daily_data = _read_dipwell_csv(csv_path=fn_WTD_daily)
    WTD_notdaily_data = _read_dipwell_csv(csv_path=fn_WTD_notdaily)

    # Order chronologically before merging

    # Merge all dipwell into one dataframe 
    WTD_data = pd.concat([WTD_notdaily_data.set_index(keys='Date'), WTD_daily_data.set_index(keys='Date')], axis=1)

    return WTD_data.reset_index()

# %% read weather data; compute ET. Return precip and ET

def fill_weather_nan_with_other_patrol_post_mean(df: pd.DataFrame) -> pd.DataFrame:
    patrol_post_names = ['KM 33', 'KM 14', 'RK']
    df_filled = df.copy(deep=True)

    patrol_post_mean = df.mean(axis=1)
    for pname in patrol_post_names:
        df[pname].fillna(value=patrol_post_mean, inplace=True)

    return df_filled

def fill_weather_nan_with_temporal_mean(df:pd.DataFrame) -> pd.DataFrame:
    patrol_post_names = ['KM 33', 'KM 14', 'RK']
    df_filled = df.copy(deep=True)

    for pname in patrol_post_names:
        pmean = df_filled.mean()[pname] # Fill with weather station mean
        df_filled[pname].fillna(value=pmean, inplace=True)

    return df_filled

def fill_weather_nan_with_mean(df: pd.DataFrame)->pd.DataFrame:
    # First, fill nan with mean of other patrol posts
    df = fill_weather_nan_with_other_patrol_post_mean(df)

    # Then, fill remaining nans with temporal avgs at each patrol post
    return fill_weather_nan_with_temporal_mean(df)

def fill_weather_nan_with_zeros(df: pd.DataFrame)->pd.DataFrame:
    # First, fill nan with mean of other patrol posts
    df = fill_weather_nan_with_other_patrol_post_mean(df)

    # Then, fill remaining nans with zeroes
    return df.fillna(value=0)

def compute_ET_df(julian_days, Tmax_data, Tmean_data, Tmin_data, rel_humidity_data, windspeed_data):
    patrol_post_names = ['KM 33', 'KM 14', 'RK']
    ET_df = pd.DataFrame(index=range(0,365), columns=[patrol_post_names])
    for i in range(0,365):
        for pname in patrol_post_names:
            ET_df.iloc[i][pname], *_ = compute_ET(julian_days[i],
                                                Tmax_data.iloc[i][pname],
                                                Tmean_data.iloc[i][pname],
                                                Tmin_data.iloc[i][pname],
                                                rel_humidity_data.iloc[i][pname],
                                                windspeed_data.iloc[i][pname])
    ET_df.columns = ET_df.columns.get_level_values(0) # Fix weird MultiIndex column thing
    return ET_df

def get_precip_and_ET(parent_folder):
    fn_weather_data = parent_folder.joinpath(r"08. Dipwell and Weather Measurements\Weather\csv\BMKG_combined_with_inside_data_for_1_year.xlsx")
    # Read weather data. Patrol posts inside the area + BMKG data from 2021 as proxy for the future
    precipitation_data = pd.read_excel(fn_weather_data, sheet_name='Rain')
    Tmean_data = pd.read_excel(fn_weather_data, sheet_name='Tmean')
    Tmin_data = pd.read_excel(fn_weather_data, sheet_name='Tmin')
    Tmax_data = pd.read_excel(fn_weather_data, sheet_name='Tmax')
    rel_humidity_data = pd.read_excel(fn_weather_data, sheet_name='Humidity')
    windspeed_data = pd.read_excel(fn_weather_data, sheet_name='WindSpeed')
    pressure_data = pd.read_excel(fn_weather_data, sheet_name='RelPressure')

    # Fill Nan with mean
    precipitation_data = fill_weather_nan_with_zeros(precipitation_data)
    Tmean_data = fill_weather_nan_with_zeros(Tmean_data)
    Tmin_data = fill_weather_nan_with_zeros(Tmin_data)
    Tmax_data = fill_weather_nan_with_zeros(Tmax_data)
    rel_humidity_data = fill_weather_nan_with_zeros(rel_humidity_data)
    windspeed_data = fill_weather_nan_with_zeros(windspeed_data)
    pressure_data = fill_weather_nan_with_zeros(pressure_data)

    julian_days = (171 + np.arange(0, 365))%365 # 171 is 20th of June, first day of computation

    ET_df = compute_ET_df(julian_days, Tmax_data, Tmean_data, Tmin_data, rel_humidity_data, windspeed_data)

    return precipitation_data, ET_df

#%% Utils
def visualize_dataframe_nodata(df):
    plt.figure(figsize=(16,10))
    sns.heatmap(df.isnull(), cbar=False, cmap="YlGnBu")
    plt.show()