#!/usr/bin/env python
# coding: utf-8

# In[1]:


import xarray as xr
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# In[2]:


noaa_precip_ds = xr.open_dataset(r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\precip_konkan.nc")
ecmwf_ds = xr.open_dataset(r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\ECMWF - skt , v10 ,u10 , msl ,tcc , 2m temp.nc")
noaa_sst_ds = xr.open_dataset(r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\noaa_sst_masked.nc")

# Compute spatial averages:
# NOAA Precipitation (dimensions: time, lat, lon)
noaa_precip = noaa_precip_ds['precip'].mean(dim=['lat', 'lon'])
# NOAA SST (dimensions: time, lat, lon)
noaa_sst = noaa_sst_ds['sst'].mean(dim=['lat', 'lon'])
# ECMWF variables (dimensions: valid_time, latitude, longitude)
# Rename the coordinate 'valid_time' to 'time' for merging without altering variable names
ecmwf_msl = ecmwf_ds['msl'].mean(dim=['latitude', 'longitude']).rename({'valid_time': 'time'})
ecmwf_tcc = ecmwf_ds['tcc'].mean(dim=['latitude', 'longitude']).rename({'valid_time': 'time'})
ecmwf_u10 = ecmwf_ds['u10'].mean(dim=['latitude', 'longitude']).rename({'valid_time': 'time'})
ecmwf_v10 = ecmwf_ds['v10'].mean(dim=['latitude', 'longitude']).rename({'valid_time': 'time'})
ecmwf_t2m = ecmwf_ds['t2m'].mean(dim=['latitude', 'longitude']).rename({'valid_time': 'time'})


# Convert each DataArray to a pandas Series
series_noaa_precip = noaa_precip.to_series()
series_noaa_sst = noaa_sst.to_series()
series_ecmwf_msl = ecmwf_msl.to_series()
series_ecmwf_tcc = ecmwf_tcc.to_series()
series_ecmwf_u10 = ecmwf_u10.to_series()
series_ecmwf_v10 = ecmwf_v10.to_series()
series_ecmwf_t2m = ecmwf_t2m.to_series()


# Merge all series into a single DataFrame using an outer join on the time index
df = pd.concat([
    series_noaa_precip,
    series_noaa_sst,
    series_ecmwf_msl,
    series_ecmwf_tcc,
    series_ecmwf_u10,
    series_ecmwf_v10,
    series_ecmwf_t2m,
    
], axis=1)

# Set column names; ECMWF variables keep their original names.
df.columns = ['noaa_precip', 'noaa_sst', 'msl', 'tcc', 'u10', 'v10', 't2m']

# Compute the correlation matrix (NaNs will be handled pairwise)
corr_matrix = df.corr()

print("Correlation Matrix:")
print(corr_matrix)


# In[3]:


print("The mosdac rows and columns are empty because it has 100% NaN values")
# Plot the correlation matrix using seaborn heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', fmt=".2f")
plt.title("Correlation Matrix Heatmap")
plt.show()


# In[ ]:




