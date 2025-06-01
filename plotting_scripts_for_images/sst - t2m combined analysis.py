#!/usr/bin/env python
# coding: utf-8

# In[4]:


import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.dates import DateFormatter
import scipy.stats as stats


# In[5]:


noaa_precip = xr.open_dataset( r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\precip_konkan.nc")
noaa_sst = xr.open_dataset( r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\noaa_sst_masked.nc")
ecmwf = xr.open_dataset( r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\ECMWF - skt , v10 ,u10 , msl ,tcc , 2m temp.nc")


# In[7]:


# Calculate spatial means for each dataset
precip_mean = noaa_precip.precip.mean(dim=['lat', 'lon'])
sst_mean = noaa_sst.sst.mean(dim=['lat', 'lon'])
ecmwf_means = {
    'msl': ecmwf.msl.mean(dim=['latitude', 'longitude']),
    'tcc': ecmwf.tcc.mean(dim=['latitude', 'longitude']),
    't2m': ecmwf.t2m.mean(dim=['latitude', 'longitude']),
    'u10': ecmwf.u10.mean(dim=['latitude', 'longitude']),
    'v10': ecmwf.v10.mean(dim=['latitude', 'longitude'])
}

# Create time series DataFrames
noaa_precip_df = pd.DataFrame({
    'time': pd.to_datetime(noaa_precip.time.values),
    'precip': precip_mean.values
}).set_index('time')

noaa_sst_df = pd.DataFrame({
    'time': pd.to_datetime(noaa_sst.time.values),
    'sst': sst_mean.values
}).set_index('time')

ecmwf_df = pd.DataFrame({
    'time': pd.to_datetime(ecmwf.valid_time.values),
    'msl': ecmwf_means['msl'].values,
    'tcc': ecmwf_means['tcc'].values,
    't2m': ecmwf_means['t2m'].values,
    'u10': ecmwf_means['u10'].values,
    'v10': ecmwf_means['v10'].values
}).set_index('time')

# Function to plot time series with correlation
def plot_correlated_series(df1, var1, df2, var2, title=None, filename=None):
    # Match the time indices
    common_df = pd.merge(
        df1[[var1]], df2[[var2]],
        left_index=True, right_index=True,
        how='inner'
    )
   
    # Calculate correlation
    correlation, p_value = stats.pearsonr(common_df[var1], common_df[var2])
   
    # Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={'height_ratios': [3, 1]})
   
    # Time series plot
    ax1.plot(common_df.index, common_df[var1], 'b-', label=var1)
    ax1_twin = ax1.twinx()
    ax1_twin.plot(common_df.index, common_df[var2], 'r-', label=var2)
   
    # Add legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1_twin.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
   
    ax1.set_ylabel(var1, color='blue')
    ax1_twin.set_ylabel(var2, color='red')
   
    if title:
        ax1.set_title(f"{title}\nCorrelation: {correlation:.3f} (p-value: {p_value:.3e})")
    else:
        ax1.set_title(f"Correlation: {correlation:.3f} (p-value: {p_value:.3e})")
   
    # Scatter plot
    ax2.scatter(common_df[var1], common_df[var2], alpha=0.6)
    ax2.set_xlabel(var1)
    ax2.set_ylabel(var2)
   
    # Add trend line
    z = np.polyfit(common_df[var1], common_df[var2], 1)
    p = np.poly1d(z)
    ax2.plot(common_df[var1], p(common_df[var1]), "r--", alpha=0.8)
   
    plt.tight_layout()
   
    plt.show()


plot_correlated_series(
    noaa_sst_df, 'sst',
    ecmwf_df, 't2m',
    title='Combined Analysis: SST and 2m Temperature',
    filename='sst_t2m_analysis.png'
)


# In[ ]:




