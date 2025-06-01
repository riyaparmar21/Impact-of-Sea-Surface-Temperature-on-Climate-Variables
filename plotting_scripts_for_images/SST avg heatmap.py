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


# In[35]:


# NOAA SST Dataset - Average over time
def plot_noaa_sst_average(ds_sst):
    """
    Plot average NOAA Sea Surface Temperature heatmap over entire time range
    """
    # Calculate time average
    sst_avg = ds_sst.sst.mean(dim='time')
   
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))
   
    # Create heatmap with RdYlBu_r colormap for temperature
    sst_plot = ax.pcolormesh(ds_sst.lon, ds_sst.lat, sst_avg,
                           cmap='RdYlBu_r', shading='auto')
   
    # Add colorbar and labels
    cbar = plt.colorbar(sst_plot, ax=ax, pad=0.05)
    cbar.set_label('Average Sea Surface Temperature (°C)', fontsize=12)
   
    ax.set_xlabel('Longitude', fontsize=12)
    ax.set_ylabel('Latitude', fontsize=12)
   
    # Add time range info to title
    start_date = np.datetime_as_string(ds_sst.time[0].values, unit="D")
    end_date = np.datetime_as_string(ds_sst.time[-1].values, unit="D")
    ax.set_title(f'NOAA Average Sea Surface Temperature\n{start_date} to {end_date}', fontsize=14)
   
    plt.tight_layout()
    plt.show()
    return fig


# In[36]:


fig_sst_avg = plot_noaa_sst_average(ds_sst)


# In[ ]:




