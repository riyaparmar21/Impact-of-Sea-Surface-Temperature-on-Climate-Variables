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


def plot_noaa_precip_grid_points(precip_file):
    """
    Create timeseries plot for NOAA Precipitation data - Grid points only
    """
    # Load the dataset
    ds = xr.open_dataset(precip_file)
   
    # Create a figure for grid points precipitation
    plt.figure(figsize=(12, 5))
   
    # Plot: Precipitation at specific points
    for i, lat in enumerate(ds.lat.values):
        for j, lon in enumerate(ds.lon.values):
            if i == 0 and j == 0:  # Only plot a selection of grid points
                label = f"lat={lat:.2f}, lon={lon:.2f}"
                ds.precip.sel(lat=lat, lon=lon).plot(label=label)
   
    plt.title('NOAA Precipitation at Selected Grid Points', fontsize=14)
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Precipitation (mm/day)', fontsize=12)
    plt.legend(loc='upper right', fontsize=10)
    plt.gca().xaxis.set_major_locator(mdates.YearLocator())
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
   
    plt.tight_layout()
    #plt.savefig('noaa_precip_grid_points.png', dpi=300)
    plt.show()
    plt.close()
   


# In[12]:


if __name__ == "__main__":
    # Define your file paths here
    noaa_precip_file = r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\precip_konkan.nc"

    # Generate individual dataset plots
    plot_noaa_precip_grid_points(noaa_precip_file)


# In[ ]:




