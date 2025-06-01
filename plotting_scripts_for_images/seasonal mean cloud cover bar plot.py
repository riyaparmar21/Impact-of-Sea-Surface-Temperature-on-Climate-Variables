#!/usr/bin/env python
# coding: utf-8

# In[1]:


import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.dates import DateFormatter
import scipy.stats as stats


# In[2]:


import seaborn as sns
import cartopy.crs as ccrs
import matplotlib.dates as mdates
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator


# In[98]:


def plot_ecmwf_vars(file_path):
    ds = xr.open_dataset(file_path)
    # For Mean Sea Level Pressure (msl)
    msl_mean = ds.msl.mean(dim=['latitude', 'longitude'])
   
    # Convert time to a more readable format for x-axis
    times = pd.to_datetime(ds.valid_time.values)
    years = [t.strftime('%Y') for t in times]
   
    # Get unique years for bars
    unique_years = sorted(set(years))
    yearly_msl = {year: [] for year in unique_years}
    
    # For Total Cloud Cover (tcc)
    # Seasonal cloud cover
    tcc_mean = ds.tcc.mean(dim=['latitude', 'longitude'])
    seasons = {'December-January-February': [12, 1, 2], 'March-April-May': [3, 4, 5],
               'June-July-August': [6, 7, 8], 'September-October-November': [9, 10, 11]}
   
    seasonal_tcc = {season: [] for season in seasons}
    for i, t in enumerate(times):
        for season, months in seasons.items():
            if t.month in months:
                seasonal_tcc[season].append(tcc_mean.values[i])
   
    seasonal_tcc_mean = {season: np.mean(vals) for season, vals in seasonal_tcc.items()}
   
    plt.figure(figsize=(8, 5))
    plt.bar(seasonal_tcc_mean.keys(), seasonal_tcc_mean.values(), color='skyblue')
    plt.xlabel('Season')
    plt.ylabel('Mean Total Cloud Cover')
    plt.title('Seasonal Mean Cloud Cover')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()
    #plt.savefig('ecmwf_tcc_seasonal.png')
    plt.close()


# In[99]:


plot_ecmwf_vars(r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\ECMWF - skt , v10 ,u10 , msl ,tcc , 2m temp.nc")


# In[ ]:




