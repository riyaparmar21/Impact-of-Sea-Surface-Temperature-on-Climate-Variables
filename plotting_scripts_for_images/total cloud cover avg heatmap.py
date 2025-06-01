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


# In[25]:


# Load your datasets
ds_precip = xr.open_dataset( r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\precip_konkan.nc")
ds_ecmwf = xr.open_dataset( r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\ECMWF - skt , v10 ,u10 , msl ,tcc , 2m temp.nc")
ds_sst = xr.open_dataset( r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\noaa_sst_masked.nc")


# In[26]:


# ECMWF Dataset - Each variable separately (averaged over time)
def plot_ecmwf_variable_average(ds_ecmwf, var_name, cmap=None):
    """
    Plot average ECMWF variables as heatmaps over entire time range
    var_name options: 'msl', 'tcc', 'u10', 'v10', 't2m'
    """
    # Set default colormap if not provided
    if cmap is None:
        cmaps = {
            'msl': 'viridis',
            'tcc': 'Blues',
            'u10': 'RdBu_r',
            'v10': 'RdBu_r',
            't2m': 'RdYlBu_r'
        }
        cmap = cmaps.get(var_name, 'viridis')
   
    # Calculate time average
    var_avg = ds_ecmwf[var_name].mean(dim='valid_time')
   
    # Get variable-specific title and colorbar label
    var_labels = {
        'msl': ('Mean Sea Level Pressure', 'Pressure (hPa)'),
        'tcc': ('Total Cloud Cover', 'Cloud Cover (%)'),
        'u10': ('U-Component of Wind at 10m', 'Wind Speed (m/s)'),
        'v10': ('V-Component of Wind at 10m', 'Wind Speed (m/s)'),
        't2m': ('Temperature at 2m', 'Temperature (K)')
    }
    title, cbar_label = var_labels.get(var_name, (var_name, var_name))
   
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))
   
    # Create heatmap
    var_plot = ax.pcolormesh(ds_ecmwf.longitude, ds_ecmwf.latitude, var_avg,
                            cmap=cmap, shading='auto')
   
    # Add colorbar and labels
    cbar = plt.colorbar(var_plot, ax=ax, pad=0.05)
    cbar.set_label(f'Average {cbar_label}', fontsize=12)
   
    ax.set_xlabel('Longitude', fontsize=12)
    ax.set_ylabel('Latitude', fontsize=12)
   
    # Add time range info to title
    start_date = np.datetime_as_string(ds_ecmwf.valid_time[0].values, unit="D")
    end_date = np.datetime_as_string(ds_ecmwf.valid_time[-1].values, unit="D")
    ax.set_title(f'Average {title}\n{start_date} to {end_date}', fontsize=14)
   
    plt.tight_layout()
    plt.show()
    return fig


# In[27]:


for var_name in ['tcc']:
    fig_var_avg = plot_ecmwf_variable_average(ds_ecmwf, var_name)


# In[ ]:




