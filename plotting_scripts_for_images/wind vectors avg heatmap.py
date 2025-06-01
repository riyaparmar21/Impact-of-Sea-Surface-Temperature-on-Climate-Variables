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


# In[33]:


# ECMWF Wind Vector Plot (Average U10 and V10)
def plot_ecmwf_wind_vectors_average(ds_ecmwf):
    """
    Plot average wind vectors using u10 and v10 components over entire time range
    """
    # Calculate time average
    u10_avg = ds_ecmwf.u10.mean(dim='valid_time')
    v10_avg = ds_ecmwf.v10.mean(dim='valid_time')
    wind_speed_avg = np.sqrt(u10_avg**2 + v10_avg**2)
   
    # Create plot with wind vectors and wind speed as background
    fig, ax = plt.subplots(figsize=(12, 10))
   
    # Plot wind speed as background
    speed_plot = ax.pcolormesh(ds_ecmwf.longitude, ds_ecmwf.latitude, wind_speed_avg,
                              cmap='viridis', shading='auto')
   
    # Add wind vectors (subsample for clarity)
    quiver_stride = 1  # adjust for desired density (smaller stride for higher density)
    lon_mesh, lat_mesh = np.meshgrid(ds_ecmwf.longitude[::quiver_stride], ds_ecmwf.latitude[::quiver_stride])
    ax.quiver(lon_mesh, lat_mesh,
             u10_avg[::quiver_stride, ::quiver_stride],
             v10_avg[::quiver_stride, ::quiver_stride],
             scale=30, color='white', alpha=0.8)
   
    # Add colorbar and labels
    cbar = plt.colorbar(speed_plot, ax=ax, pad=0.05)
    cbar.set_label('Average Wind Speed (m/s)', fontsize=12)
   
    ax.set_xlabel('Longitude', fontsize=12)
    ax.set_ylabel('Latitude', fontsize=12)
   
    # Add time range info to title
    start_date = np.datetime_as_string(ds_ecmwf.valid_time[0].values, unit="D")
    end_date = np.datetime_as_string(ds_ecmwf.valid_time[-1].values, unit="D")
    ax.set_title(f'Average Wind Vectors\n{start_date} to {end_date}', fontsize=14)
   
    plt.tight_layout()
    plt.show()
    return fig


# In[34]:


fig_wind_avg = plot_ecmwf_wind_vectors_average(ds_ecmwf)


# In[ ]:




