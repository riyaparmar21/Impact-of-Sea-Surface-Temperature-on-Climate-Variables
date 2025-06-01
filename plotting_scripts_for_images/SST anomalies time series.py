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


# In[11]:


#############################################
# NOAA SST Dataset Plotting Code
#############################################

def plot_noaa_sst(sst_file):
    """
    Create timeseries plots for NOAA Sea Surface Temperature data
    """
    # Load the dataset
    ds = xr.open_dataset(sst_file)
   
    # Create area-averaged SST timeseries
    plt.figure(figsize=(12, 6))
    sst_mean = ds.sst.mean(dim=['lat', 'lon'])
   
    # Convert to Celsius if it appears to be in Kelvin
    if sst_mean.mean() > 100:  # Likely in Kelvin
        sst_mean = sst_mean - 273.15
        units = '°C'
    else:
        units = '°C'  # Assuming it's already in Celsius
   
    # Create anomaly plot (difference from climatological mean)
    climatology = ds.sst.groupby('time.month').mean()
    anomalies = ds.sst.groupby('time.month') - climatology
    anom_mean = anomalies.mean(dim=['lat', 'lon'])
   
    plt.figure(figsize=(12, 6))
    anom_mean.plot(linewidth=1.5, color='purple')
   
    plt.title('NOAA SST Anomalies (Difference from Climatological Mean)', fontsize=14)
    plt.xlabel('Time', fontsize=12)
    plt.ylabel(f'Temperature Anomaly ({units})', fontsize=12)
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    #plt.savefig('noaa_sst_anomalies.png', dpi=300)
    plt.close()


# In[12]:


if __name__ == "__main__":

    noaa_sst_file =  r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\noaa_sst_masked.nc"
    
    plot_noaa_sst(noaa_sst_file)


# In[ ]:




