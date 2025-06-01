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


# In[35]:


def plot_ecmwf_variable_timeseries(ecmwf_file, var):
    """
    Create timeseries plot for a specific ECMWF variable
   
    Parameters:
    -----------
    ecmwf_file : str
        Path to the ECMWF data file
    var : str
        Variable to plot ('msl', 'tcc', 'u10', 'v10', 't2m')
    """
    # Load the dataset
    ds = xr.open_dataset(ecmwf_file)
   
    # Dictionary of variable information
    titles = {
        'msl': 'Mean Sea Level Pressure',
        'tcc': 'Total Cloud Cover',
        'u10': '10m U Wind Component',
        'v10': '10m V Wind Component',
        't2m': '2m Temperature'
    }
    units = {
        'msl': 'hPa',
        'tcc': 'fraction',
        'u10': 'm/s',
        'v10': 'm/s',
        't2m': 'K'
    }
   
    if var not in ds:
        print(f"Variable '{var}' not found in dataset")
        return
       
    # Create figure for area average timeseries
    plt.figure(figsize=(12, 6))
   
    # Convert msl from Pa to hPa if needed
    if var == 'msl' and ds[var].max() > 10000:  # Likely in Pa
        area_mean = ds[var].mean(dim=['latitude', 'longitude']) / 100  # Convert to hPa
    else:
        area_mean = ds[var].mean(dim=['latitude', 'longitude'])
   
    # For temperature, convert to Celsius if it appears to be in Kelvin
    if var == 't2m' and area_mean.mean() > 100:  # Likely in Kelvin
        area_mean = area_mean - 273.15  # Convert to Celsius
        units[var] = '°C'  # Update unit
   
    area_mean.plot(linewidth=1.5, color=sns.color_palette("viridis", 5)[['msl', 'tcc', 'u10', 'v10', 't2m'].index(var)])
   
    plt.title(f'ECMWF {titles[var]} (Area Average)', fontsize=14)
    plt.xlabel('Time', fontsize=12)
   
    # Adjust y-axis label based on variable
    if var == 't2m' and units[var] == '°C':
        plt.ylabel(f'Temperature ({units[var]})', fontsize=12)
    else:
        plt.ylabel(f'{titles[var]} ({units[var]})', fontsize=12)
   
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    #plt.savefig(f'ecmwf_{var}_timeseries.png', dpi=300)
    plt.close()
   
    print(f"ECMWF {titles[var]} timeseries plot created successfully!")


def plot_ecmwf_variable_seasonal(ecmwf_file, var):
    """
    Create seasonal cycle plot for a specific ECMWF variable
   
    Parameters:
    -----------
    ecmwf_file : str
        Path to the ECMWF data file
    var : str
        Variable to plot ('msl', 'tcc', 'u10', 'v10', 't2m')
    """
    # Load the dataset
    ds = xr.open_dataset(ecmwf_file)
   
    # Dictionary of variable information
    titles = {
        'msl': 'Mean Sea Level Pressure',
        'tcc': 'Total Cloud Cover',
        'u10': '10m U Wind Component',
        'v10': '10m V Wind Component',
        't2m': '2m Temperature'
    }
    units = {
        'msl': 'hPa',
        'tcc': 'fraction',
        'u10': 'm/s',
        'v10': 'm/s',
        't2m': 'K'
    }
   
    if var not in ds:
        print(f"Variable '{var}' not found in dataset")
        return
   
    # Create seasonal cycle plot
    monthly_mean = ds[var].groupby('valid_time.month').mean()
   
    plt.figure(figsize=(10, 6))
   
    # Apply conversions if needed
    if var == 'msl' and monthly_mean.mean() > 10000:  # Likely in Pa
        plot_data = monthly_mean.mean(dim=['latitude', 'longitude']) / 100  # Convert to hPa
    else:
        plot_data = monthly_mean.mean(dim=['latitude', 'longitude'])
   
    # For temperature, convert to Celsius if it appears to be in Kelvin
    if var == 't2m' and plot_data.mean() > 100:  # Likely in Kelvin
        plot_data = plot_data - 273.15  # Convert to Celsius
        units[var] = '°C'  # Update unit
   
    plot_data.plot(marker='o', color=sns.color_palette("viridis", 5)[['msl', 'tcc', 'u10', 'v10', 't2m'].index(var)])
   
    plt.title(f'ECMWF {titles[var]} Seasonal Cycle', fontsize=14)
    plt.xlabel('Month', fontsize=12)
   
    if var == 't2m' and units[var] == '°C':
        plt.ylabel(f'Temperature ({units[var]})', fontsize=12)
    else:
        plt.ylabel(f'{titles[var]} ({units[var]})', fontsize=12)
   
    plt.xticks(range(1, 13), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                            'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    #plt.savefig(f'ecmwf_{var}_seasonal.png', dpi=300)
    plt.close()
   
    print(f"ECMWF {titles[var]} seasonal cycle plot created successfully!")


def plot_ecmwf_wind_speed(ecmwf_file):
    """
    Create wind speed plot from ECMWF u10 and v10 components
   
    Parameters:
    -----------
    ecmwf_file : str
        Path to the ECMWF data file
    """
    # Load the dataset
    ds = xr.open_dataset(ecmwf_file)
   
    if 'u10' not in ds or 'v10' not in ds:
        print("Wind components (u10, v10) not found in dataset")
        return
   
    # Calculate area-averaged wind speed and direction
    u_mean = ds['u10'].mean(dim=['latitude', 'longitude'])
    v_mean = ds['v10'].mean(dim=['latitude', 'longitude'])
    wind_speed = np.sqrt(u_mean**2 + v_mean**2)
   
    fig, ax1 = plt.subplots(figsize=(12, 6))
   
    # Plot wind speed
    ax1.plot(ds.valid_time, wind_speed, color='darkblue', linewidth=1.5)
    ax1.set_xlabel('Time', fontsize=12)
    ax1.set_ylabel('Wind Speed (m/s)', fontsize=12)
    ax1.set_title('ECMWF 10m Wind Speed (Area Average)', fontsize=14)
   
    plt.tight_layout()
    plt.show()
    #plt.savefig('ecmwf_wind_speed_timeseries.png', dpi=300)
    plt.close()
   
    print("ECMWF wind speed plot created successfully!")


def plot_ecmwf_data(ecmwf_file, var=None, plot_type=None):
    """
    Create plots for ECMWF variables
   
    Parameters:
    -----------
    ecmwf_file : str
        Path to the ECMWF data file
    var : str or None
        Variable to plot ('msl', 'tcc', 'u10', 'v10', 't2m', 'wind_speed', 'all')
        If None or 'all', plots all variables
    plot_type : str or None
        Type of plot ('timeseries', 'seasonal', 'all')
        If None or 'all', creates all plot types for the specified variable(s)
    """
    variables = ['msl', 'tcc', 'u10', 'v10', 't2m']
   
    # Handle the case when var is None or 'all'
    if var is None or var == 'all':
        # Plot all variables with all plot types
        for variable in variables:
            plot_ecmwf_variable_timeseries(ecmwf_file, variable)
            plot_ecmwf_variable_seasonal(ecmwf_file, variable)
       
        # Also plot wind speed
        plot_ecmwf_wind_speed(ecmwf_file)
        return
   
    # Handle wind_speed separately
    if var == 'wind_speed':
        plot_ecmwf_wind_speed(ecmwf_file)
        return
   
    # Handle individual variables
    if var in variables:
        if plot_type is None or plot_type == 'all':
            # Create both timeseries and seasonal plots
            plot_ecmwf_variable_timeseries(ecmwf_file, var)
            plot_ecmwf_variable_seasonal(ecmwf_file, var)
        elif plot_type == 'timeseries':
            plot_ecmwf_variable_timeseries(ecmwf_file, var)
        elif plot_type == 'seasonal':
            plot_ecmwf_variable_seasonal(ecmwf_file, var)
        else:
            print(f"Invalid plot_type: {plot_type}")
            print("Valid options: 'timeseries', 'seasonal', 'all'")
    else:
        print(f"Invalid variable: {var}")
        print(f"Valid options: {variables + ['wind_speed', 'all']}")


# In[43]:


if __name__ == "__main__":
    
    ecmwf_file =  r"C:\Users\Admin\RIYA PROJECT\DATASETS\main datasets to use\ECMWF - skt , v10 ,u10 , msl ,tcc , 2m temp.nc"

    plot_ecmwf_data(ecmwf_file, var='tcc', plot_type='seasonal')


# In[ ]:


# several ways to use it:

# 1. Plot a specific variable with a specific plot type:

#    plot_ecmwf_data(ecmwf_file, var='t2m', plot_type='timeseries')  # Only 2m temperature timeseries
#    plot_ecmwf_data(ecmwf_file, var='msl', plot_type='seasonal')    # Only mean sea level pressure seasonal cycle


# 2. Plot all types for a specific variable:

#    plot_ecmwf_data(ecmwf_file, var='tcc')  # Both timeseries and seasonal for total cloud cover


# 3. Plot just the wind speed:

#    plot_ecmwf_data(ecmwf_file, var='wind_speed')  # Only wind speed plot
#    ```

# 4. Plot everything (original behavior):

#    plot_ecmwf_data(ecmwf_file)  # Plots all variables with all plot types


# 5. Or call the specific functions directly:

#    plot_ecmwf_variable_timeseries(ecmwf_file, 'u10')  # Just u10 timeseries
#    plot_ecmwf_variable_seasonal(ecmwf_file, 'v10')    # Just v10 seasonal cycle
#    plot_ecmwf_wind_speed(ecmwf_file)                  # Just wind speed

