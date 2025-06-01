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


# In[9]:


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
   
    # Create seasonal cycle plot
    sst_seasonal = ds.sst.groupby('time.month').mean()
   
    plt.figure(figsize=(10, 6))
    sst_seasonal.mean(dim=['lat', 'lon']).plot(marker='o', color='darkred')
   
    plt.title('NOAA SST Seasonal Cycle', fontsize=14)
    plt.xlabel('Month', fontsize=12)
    plt.ylabel(f'Temperature ({units})', fontsize=12)
    plt.xticks(range(1, 13), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                            'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    #plt.savefig('noaa_sst_seasonal.png', dpi=300)
    plt.close()


# In[10]:


if __name__ == "__main__":

    noaa_sst_file =  r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\noaa_sst_masked.nc"
    
    plot_noaa_sst(noaa_sst_file)


# In[ ]:




