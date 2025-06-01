#!/usr/bin/env python
# coding: utf-8

# In[4]:


import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.dates import DateFormatter
import scipy.stats as stats


# In[13]:


sst_data = xr.open_dataset(r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\noaa_sst_masked.nc")
ecmwf_data = xr.open_dataset(r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\ECMWF - skt , v10 ,u10 , msl ,tcc , 2m temp.nc")

# Ensure same time period
sst = sst_data.sst.sel(time=slice('2019-01-01', '2024-12-01'))
tcc = ecmwf_data.tcc.sel(valid_time=slice('2019-01-01', '2024-12-01'))

# Interpolate to common grid
tcc_interp = tcc.interp(latitude=sst.lat, longitude=sst.lon)
tcc_interp = tcc_interp.rename({'valid_time': 'time'})

# Compute correlation at each grid point
correlation = xr.corr(sst, tcc_interp, dim='time')

# Plot
plt.figure(figsize=(12, 10))
im = plt.pcolormesh(sst.lon, sst.lat, correlation, cmap='RdBu_r', vmin=-1, vmax=1)
plt.colorbar(im, label='Correlation Coefficient')
plt.title('Spatial Correlation: SST vs Total Cloud Cover')

# Add correlation values as text on the plot
for i in range(len(sst.lat)):
    for j in range(len(sst.lon)):
        # Only show values at some grid points to avoid overcrowding
        if i % 2 == 0 and j % 2 == 0:  # Show every other point
            plt.text(sst.lon[j], sst.lat[i], 
                     f'{correlation.values[i, j]:.2f}',  # Format to 2 decimal places
                     ha='center', va='center', 
                     color='black', fontsize=8,
                     bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', pad=1))

plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.tight_layout()
plt.show()


# In[ ]:




