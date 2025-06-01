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



def plot_noaa_precip_seasonal_cycle(precip_file):
    """
    Create seasonal cycle plot for NOAA Precipitation data
    """
    # Load the dataset
    ds = xr.open_dataset(precip_file)
   
    # Create seasonal cycle plot
    precip_seasonal = ds.precip.groupby('time.month').mean()
   
    plt.figure(figsize=(10, 6))
    precip_seasonal.mean(dim=['lat', 'lon']).plot(marker='o')
    plt.title('NOAA Precipitation Seasonal Cycle', fontsize=14)
    plt.xlabel('Month', fontsize=12)
    plt.ylabel('Precipitation (mm/day)', fontsize=12)
    plt.xticks(range(1, 13), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                            'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    plt.grid(True)
    plt.tight_layout()
    #plt.savefig('noaa_precip_seasonal.png', dpi=300)
    plt.show()
    plt.close()
   


# In[10]:


if __name__ == "__main__":
    # Define your file paths here
    noaa_precip_file = r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\precip_konkan.nc"

    # Generate individual dataset plots
    plot_noaa_precip_seasonal_cycle(noaa_precip_file)


# In[ ]:




