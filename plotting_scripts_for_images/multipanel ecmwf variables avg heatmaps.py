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


# In[41]:


# Multi-Variable ECMWF Panel (all in one figure) - Averaged
def plot_ecmwf_panel_average(ds_ecmwf):
    """
    Plot all averaged ECMWF variables in a single panel
    """
    # Define variable names and their associated cmaps
    vars_and_cmaps = {
        'msl': 'viridis',
        'tcc': 'Blues',
        't2m': 'RdYlBu_r',
    }
   
    # Create figure with subplots
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 14))
    axes = axes.flatten()
   
    # Plot each variable
    for i, (var_name, cmap) in enumerate(vars_and_cmaps.items()):
        if i < len(axes) - 1:  # One spot reserved for wind vectors
            var_avg = ds_ecmwf[var_name].mean(dim='valid_time')
            var_plot = axes[i].pcolormesh(ds_ecmwf.longitude, ds_ecmwf.latitude, var_avg,
                                       cmap=cmap, shading='auto')
            cbar = plt.colorbar(var_plot, ax=axes[i], pad=0.05)
           
            # Get variable labels
            var_labels = {
                'msl': ('Mean Sea Level Pressure', 'Pressure (hPa)'),
                'tcc': ('Total Cloud Cover', 'Cloud Cover (%)'),
                't2m': ('Temperature at 2m', 'Temperature (K)')
            }
            title, cbar_label = var_labels.get(var_name, (var_name, var_name))
           
            cbar.set_label(f'Average {cbar_label}', fontsize=10)
            axes[i].set_title(f'Average {title}', fontsize=12)
            axes[i].set_xlabel('Longitude', fontsize=10)
            axes[i].set_ylabel('Latitude', fontsize=10)
   
    # Plot average wind vectors in the last panel
    u10_avg = ds_ecmwf.u10.mean(dim='valid_time')
    v10_avg = ds_ecmwf.v10.mean(dim='valid_time')
    wind_speed_avg = np.sqrt(u10_avg**2 + v10_avg**2)
   
    # Plot wind speed as background
    speed_plot = axes[-1].pcolormesh(ds_ecmwf.longitude, ds_ecmwf.latitude, wind_speed_avg,
                                  cmap='viridis', shading='auto')
   
    # Add wind vectors
    quiver_stride = 1  # Smaller stride for higher vector density
    lon_mesh, lat_mesh = np.meshgrid(ds_ecmwf.longitude[::quiver_stride], ds_ecmwf.latitude[::quiver_stride])
    axes[-1].quiver(lon_mesh, lat_mesh,
                 u10_avg[::quiver_stride, ::quiver_stride],
                 v10_avg[::quiver_stride, ::quiver_stride],
                 scale=30, color='white', alpha=0.8)
   
    cbar = plt.colorbar(speed_plot, ax=axes[-1], pad=0.05)
    cbar.set_label('Average Wind Speed (m/s)', fontsize=10)
    axes[-1].set_title('Average Wind Vectors (u10, v10)', fontsize=12)
    axes[-1].set_xlabel('Longitude', fontsize=10)
    axes[-1].set_ylabel('Latitude', fontsize=10)
   
    # Add time range to main title
    start_date = np.datetime_as_string(ds_ecmwf.valid_time[0].values, unit="D")
    end_date = np.datetime_as_string(ds_ecmwf.valid_time[-1].values, unit="D")
    plt.suptitle(f'ECMWF Average Variables\n{start_date} to {end_date}', fontsize=16)
   
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()
    return fig


# In[42]:


fig_panel_avg = plot_ecmwf_panel_average(ds_ecmwf)


# In[ ]:




