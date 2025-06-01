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


# In[51]:


# Load your datasets
ds_precip = xr.open_dataset( r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\precip_konkan.nc")
ds_ecmwf = xr.open_dataset( r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\ECMWF - skt , v10 ,u10 , msl ,tcc , 2m temp.nc")
ds_sst = xr.open_dataset( r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\noaa_sst_masked.nc")


# In[60]:


# Seasonal Averages 
def plot_seasonal_averages(ds, var_name, time_dim='valid_time', lat_dim='latitude', lon_dim='longitude'):
    """
    Plot seasonal averages for a given variable
    """
    # Convert dataset time to datetime if not already
    if not np.issubdtype(ds[time_dim].dtype, np.datetime64):
        ds[time_dim] = ds[time_dim].astype('datetime64[ns]')
   
    # Create a copy of the dataset with a 'season' coordinate
    ds_with_season = ds.copy()
   
    # Add season information
    month_to_season = {1: 'Dec-Jan-Feb', 2: 'Dec-Jan-Feb,', 3: 'Mar-Apr-May', 4: 'Mar-Apr-May', 5: 'Mar-Apr-May',
                       6: 'Jun-Jul-Aug', 7: 'Jun-Jul-Aug', 8: 'Jun-Jul-Aug', 9: 'Sep-Oct-Nov', 10: 'Sep-Oct-Nov',
                       11: 'Sep-Oct-Nov', 12: 'Dec-Jan-Feb,'}
   
    # Extract month and map to season
    months = ds[time_dim].dt.month.values
    seasons = [month_to_season[m] for m in months]
    ds_with_season = ds_with_season.assign_coords(season=('valid_time', seasons))
   
    # Calculate seasonal averages
    seasonal_avg = ds_with_season[var_name].groupby('season').mean(dim=time_dim)
   
    # Create figure with 4 subplots (one for each season)
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 14))
    axes = axes.flatten()
   
    # Define color maps based on variable
    if 'precip' in var_name:
        cmap = 'Blues'
        label = 'Precipitation (mm)'
    elif 'sst' in var_name:
        cmap = 'RdYlBu_r'
        label = 'Temperature (°C)'
    elif 't2m' in var_name:
        cmap = 'RdYlBu_r'
        label = 'Temperature (K)'
    else:
        cmap = 'viridis'
        label = var_name
   
    # Plot each season
    for i, season in enumerate(['Dec-Jan-Feb', 'Mar-Apr-May', 'Jun-Jul-Aug', 'Sep-Oct-Nov']):
        if season in seasonal_avg.season.values:
            season_data = seasonal_avg.sel(season=season)
           
            # Get coordinate names for proper plotting
            lon_coords = ds[lon_dim]
            lat_coords = ds[lat_dim]
           
            # Create heatmap
            im = axes[i].pcolormesh(lon_coords, lat_coords, season_data,
                                  cmap=cmap, shading='auto')
           
            # Add colorbar
            cbar = plt.colorbar(im, ax=axes[i], pad=0.05)
            cbar.set_label(f'Average {label}', fontsize=10)
           
            axes[i].set_title(f'{season} Average {var_name.upper()}', fontsize=12)
            axes[i].set_xlabel('Longitude', fontsize=10)
            axes[i].set_ylabel('Latitude', fontsize=10)
   
    plt.suptitle(f'Seasonal Average {var_name.upper()}', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()
    return fig


# In[61]:


fig_u10_seasonal = plot_seasonal_averages(ds_ecmwf, 'u10')


# In[ ]:




