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


# In[94]:


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
    
    # For Temperature at 2m (t2m)
    # Convert Kelvin to Celsius
    t2m_mean = ds.t2m.mean(dim=['latitude', 'longitude']) - 273.15
   
    # Plot temperature as monthly bars
    months = [t.strftime('%Y-%m') for t in times]
   
    # Select just the first 24 months for better readability
    plt.figure(figsize=(14,9))
    plt.bar(months[:60], t2m_mean[:60], color='red')
    plt.xlabel('Month')
    plt.ylabel('Mean Temperature (°C)')
    plt.title('Monthly Mean Temperature')
    plt.xticks(rotation=90)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()
    #plt.savefig('ecmwf_t2m_monthly.png')
    plt.close()


# In[95]:


plot_ecmwf_vars(r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\ECMWF - skt , v10 ,u10 , msl ,tcc , 2m temp.nc")


# In[ ]:




