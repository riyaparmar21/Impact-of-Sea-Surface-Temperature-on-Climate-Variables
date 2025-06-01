#!/usr/bin/env python
# coding: utf-8

# In[4]:


import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.dates import DateFormatter
import scipy.stats as stats


# In[15]:


import seaborn as sns
import cartopy.crs as ccrs
import matplotlib.dates as mdates
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator


# In[18]:


# Set plotting style
plt.style.use('seaborn-whitegrid')

# 1. NOAA Precipitation Dataset - Average over time
def plot_noaa_precip_average(ds_precip):
    """
    Plot average NOAA precipitation heatmap over entire time range
    """
    # Calculate time average
    precip_avg = ds_precip.precip.mean(dim='time')
   
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))
   
    # Custom blue colormap for precipitation
    cmap_precip = LinearSegmentedColormap.from_list('precip_blues',
                                                    ['#FFFFFF', '#D1EEFC', '#75C6EF', '#1E90FF', '#0066CC', '#003366'])
   
    # Create heatmap
    precip_plot = ax.pcolormesh(ds_precip.lon, ds_precip.lat, precip_avg,
                               cmap=cmap_precip, shading='auto')
   
    # Add colorbar and labels
    cbar = plt.colorbar(precip_plot, ax=ax, pad=0.05)
    cbar.set_label('Average Precipitation (mm)', fontsize=12)
   
    ax.set_xlabel('Longitude', fontsize=12)
    ax.set_ylabel('Latitude', fontsize=12)
   
    # Add time range info to title
    start_date = np.datetime_as_string(ds_precip.time[0].values, unit="D")
    end_date = np.datetime_as_string(ds_precip.time[-1].values, unit="D")
    ax.set_title(f'NOAA Average Precipitation\n{start_date} to {end_date}', fontsize=14)
   
    plt.tight_layout()
    plt.show()
    return fig

# Load your datasets
ds_precip = xr.open_dataset( r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\precip_konkan.nc")
ds_ecmwf = xr.open_dataset( r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\ECMWF - skt , v10 ,u10 , msl ,tcc , 2m temp.nc")
ds_sst = xr.open_dataset( r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\noaa_sst_masked.nc")


# In[19]:


fig_precip_avg = plot_noaa_precip_average(ds_precip)


# In[ ]:




