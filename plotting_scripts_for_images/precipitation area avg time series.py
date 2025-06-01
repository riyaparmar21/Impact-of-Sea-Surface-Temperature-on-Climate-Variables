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


# In[3]:


def plot_noaa_precip_area_average(precip_file):
    """
    Create timeseries plot for NOAA Precipitation data - Area average only
    """
    # Load the dataset
    ds = xr.open_dataset(precip_file)
   
    # Create a figure for area-averaged precipitation
    plt.figure(figsize=(12, 5))
   
    # Plot: Area-averaged precipitation timeseries
    precip_mean = ds.precip.mean(dim=['lat', 'lon'])
    precip_mean.plot(color='blue', linewidth=1.5)
   
    plt.title('NOAA Area-averaged Precipitation (2019-2025)', fontsize=14)
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Precipitation (mm/day)', fontsize=12)
    plt.gca().xaxis.set_major_locator(mdates.YearLocator())
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
   
    plt.tight_layout()
    #plt.savefig('noaa_precip_area_average.png', dpi=300)
    plt.show()
    plt.close()


# In[5]:


if __name__ == "__main__":
    # Define your file paths here
    noaa_precip_file = r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\precip_konkan.nc"

    # Generate individual dataset plots
    plot_noaa_precip_area_average(noaa_precip_file)


# In[ ]:




