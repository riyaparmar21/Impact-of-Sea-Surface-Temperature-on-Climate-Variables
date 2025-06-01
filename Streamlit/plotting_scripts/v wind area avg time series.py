# #!/usr/bin/env python
# # coding: utf-8

# # In[1]:


# import xarray as xr
# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd
# from matplotlib.dates import DateFormatter
# import scipy.stats as stats


# # In[2]:


# import seaborn as sns
# import cartopy.crs as ccrs
# import matplotlib.dates as mdates
# from matplotlib.colors import LinearSegmentedColormap
# from matplotlib.ticker import MaxNLocator


# # In[35]:


# def plot_ecmwf_variable_timeseries(ecmwf_file, var):
#     """
#     Create timeseries plot for a specific ECMWF variable
   
#     Parameters:
#     -----------
#     ecmwf_file : str
#         Path to the ECMWF data file
#     var : str
#         Variable to plot ('msl', 'tcc', 'u10', 'v10', 't2m')
#     """
#     # Load the dataset
#     ds = xr.open_dataset(ecmwf_file)
   
#     # Dictionary of variable information
#     titles = {
#         'msl': 'Mean Sea Level Pressure',
#         'tcc': 'Total Cloud Cover',
#         'u10': '10m U Wind Component',
#         'v10': '10m V Wind Component',
#         't2m': '2m Temperature'
#     }
#     units = {
#         'msl': 'hPa',
#         'tcc': 'fraction',
#         'u10': 'm/s',
#         'v10': 'm/s',
#         't2m': 'K'
#     }
   
#     if var not in ds:
#         print(f"Variable '{var}' not found in dataset")
#         return
       
#     # Create figure for area average timeseries
#     plt.figure(figsize=(12, 6))
   
#     # Convert msl from Pa to hPa if needed
#     if var == 'msl' and ds[var].max() > 10000:  # Likely in Pa
#         area_mean = ds[var].mean(dim=['latitude', 'longitude']) / 100  # Convert to hPa
#     else:
#         area_mean = ds[var].mean(dim=['latitude', 'longitude'])
   
#     # For temperature, convert to Celsius if it appears to be in Kelvin
#     if var == 't2m' and area_mean.mean() > 100:  # Likely in Kelvin
#         area_mean = area_mean - 273.15  # Convert to Celsius
#         units[var] = '°C'  # Update unit
   
#     area_mean.plot(linewidth=1.5, color=sns.color_palette("viridis", 5)[['msl', 'tcc', 'u10', 'v10', 't2m'].index(var)])
   
#     plt.title(f'ECMWF {titles[var]} (Area Average)', fontsize=14)
#     plt.xlabel('Time', fontsize=12)
   
#     # Adjust y-axis label based on variable
#     if var == 't2m' and units[var] == '°C':
#         plt.ylabel(f'Temperature ({units[var]})', fontsize=12)
#     else:
#         plt.ylabel(f'{titles[var]} ({units[var]})', fontsize=12)
   
#     plt.grid(True)
#     plt.tight_layout()
#     plt.show()
#     #plt.savefig(f'ecmwf_{var}_timeseries.png', dpi=300)
#     plt.close()
   
#     print(f"ECMWF {titles[var]} timeseries plot created successfully!")


# def plot_ecmwf_variable_seasonal(ecmwf_file, var):
#     """
#     Create seasonal cycle plot for a specific ECMWF variable
   
#     Parameters:
#     -----------
#     ecmwf_file : str
#         Path to the ECMWF data file
#     var : str
#         Variable to plot ('msl', 'tcc', 'u10', 'v10', 't2m')
#     """
#     # Load the dataset
#     ds = xr.open_dataset(ecmwf_file)
   
#     # Dictionary of variable information
#     titles = {
#         'msl': 'Mean Sea Level Pressure',
#         'tcc': 'Total Cloud Cover',
#         'u10': '10m U Wind Component',
#         'v10': '10m V Wind Component',
#         't2m': '2m Temperature'
#     }
#     units = {
#         'msl': 'hPa',
#         'tcc': 'fraction',
#         'u10': 'm/s',
#         'v10': 'm/s',
#         't2m': 'K'
#     }
   
#     if var not in ds:
#         print(f"Variable '{var}' not found in dataset")
#         return
   
#     # Create seasonal cycle plot
#     monthly_mean = ds[var].groupby('valid_time.month').mean()
   
#     plt.figure(figsize=(10, 6))
   
#     # Apply conversions if needed
#     if var == 'msl' and monthly_mean.mean() > 10000:  # Likely in Pa
#         plot_data = monthly_mean.mean(dim=['latitude', 'longitude']) / 100  # Convert to hPa
#     else:
#         plot_data = monthly_mean.mean(dim=['latitude', 'longitude'])
   
#     # For temperature, convert to Celsius if it appears to be in Kelvin
#     if var == 't2m' and plot_data.mean() > 100:  # Likely in Kelvin
#         plot_data = plot_data - 273.15  # Convert to Celsius
#         units[var] = '°C'  # Update unit
   
#     plot_data.plot(marker='o', color=sns.color_palette("viridis", 5)[['msl', 'tcc', 'u10', 'v10', 't2m'].index(var)])
   
#     plt.title(f'ECMWF {titles[var]} Seasonal Cycle', fontsize=14)
#     plt.xlabel('Month', fontsize=12)
   
#     if var == 't2m' and units[var] == '°C':
#         plt.ylabel(f'Temperature ({units[var]})', fontsize=12)
#     else:
#         plt.ylabel(f'{titles[var]} ({units[var]})', fontsize=12)
   
#     plt.xticks(range(1, 13), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
#                             'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
#     plt.grid(True)
#     plt.tight_layout()
#     plt.show()
#     #plt.savefig(f'ecmwf_{var}_seasonal.png', dpi=300)
#     plt.close()
   
#     print(f"ECMWF {titles[var]} seasonal cycle plot created successfully!")


# def plot_ecmwf_wind_speed(ecmwf_file):
#     """
#     Create wind speed plot from ECMWF u10 and v10 components
   
#     Parameters:
#     -----------
#     ecmwf_file : str
#         Path to the ECMWF data file
#     """
#     # Load the dataset
#     ds = xr.open_dataset(ecmwf_file)
   
#     if 'u10' not in ds or 'v10' not in ds:
#         print("Wind components (u10, v10) not found in dataset")
#         return
   
#     # Calculate area-averaged wind speed and direction
#     u_mean = ds['u10'].mean(dim=['latitude', 'longitude'])
#     v_mean = ds['v10'].mean(dim=['latitude', 'longitude'])
#     wind_speed = np.sqrt(u_mean**2 + v_mean**2)
   
#     fig, ax1 = plt.subplots(figsize=(12, 6))
   
#     # Plot wind speed
#     ax1.plot(ds.valid_time, wind_speed, color='darkblue', linewidth=1.5)
#     ax1.set_xlabel('Time', fontsize=12)
#     ax1.set_ylabel('Wind Speed (m/s)', fontsize=12)
#     ax1.set_title('ECMWF 10m Wind Speed (Area Average)', fontsize=14)
   
#     plt.tight_layout()
#     plt.show()
#     #plt.savefig('ecmwf_wind_speed_timeseries.png', dpi=300)
#     plt.close()
   
#     print("ECMWF wind speed plot created successfully!")


# def plot_ecmwf_data(ecmwf_file, var=None, plot_type=None):
#     """
#     Create plots for ECMWF variables
   
#     Parameters:
#     -----------
#     ecmwf_file : str
#         Path to the ECMWF data file
#     var : str or None
#         Variable to plot ('msl', 'tcc', 'u10', 'v10', 't2m', 'wind_speed', 'all')
#         If None or 'all', plots all variables
#     plot_type : str or None
#         Type of plot ('timeseries', 'seasonal', 'all')
#         If None or 'all', creates all plot types for the specified variable(s)
#     """
#     variables = ['msl', 'tcc', 'u10', 'v10', 't2m']
   
#     # Handle the case when var is None or 'all'
#     if var is None or var == 'all':
#         # Plot all variables with all plot types
#         for variable in variables:
#             plot_ecmwf_variable_timeseries(ecmwf_file, variable)
#             plot_ecmwf_variable_seasonal(ecmwf_file, variable)
       
#         # Also plot wind speed
#         plot_ecmwf_wind_speed(ecmwf_file)
#         return
   
#     # Handle wind_speed separately
#     if var == 'wind_speed':
#         plot_ecmwf_wind_speed(ecmwf_file)
#         return
   
#     # Handle individual variables
#     if var in variables:
#         if plot_type is None or plot_type == 'all':
#             # Create both timeseries and seasonal plots
#             plot_ecmwf_variable_timeseries(ecmwf_file, var)
#             plot_ecmwf_variable_seasonal(ecmwf_file, var)
#         elif plot_type == 'timeseries':
#             plot_ecmwf_variable_timeseries(ecmwf_file, var)
#         elif plot_type == 'seasonal':
#             plot_ecmwf_variable_seasonal(ecmwf_file, var)
#         else:
#             print(f"Invalid plot_type: {plot_type}")
#             print("Valid options: 'timeseries', 'seasonal', 'all'")
#     else:
#         print(f"Invalid variable: {var}")
#         print(f"Valid options: {variables + ['wind_speed', 'all']}")


# # In[46]:


# if __name__ == "__main__":
    
#     ecmwf_file =  r"ecmwf datas.nc"

#     plot_ecmwf_data(ecmwf_file, var='v10', plot_type='timeseries')


# # In[ ]:


# # several ways to use it:

# # 1. Plot a specific variable with a specific plot type:

# #    plot_ecmwf_data(ecmwf_file, var='t2m', plot_type='timeseries')  # Only 2m temperature timeseries
# #    plot_ecmwf_data(ecmwf_file, var='msl', plot_type='seasonal')    # Only mean sea level pressure seasonal cycle


# # 2. Plot all types for a specific variable:

# #    plot_ecmwf_data(ecmwf_file, var='tcc')  # Both timeseries and seasonal for total cloud cover


# # 3. Plot just the wind speed:

# #    plot_ecmwf_data(ecmwf_file, var='wind_speed')  # Only wind speed plot
# #    ```

# # 4. Plot everything (original behavior):

# #    plot_ecmwf_data(ecmwf_file)  # Plots all variables with all plot types


# # 5. Or call the specific functions directly:

# #    plot_ecmwf_variable_timeseries(ecmwf_file, 'u10')  # Just u10 timeseries
# #    plot_ecmwf_variable_seasonal(ecmwf_file, 'v10')    # Just v10 seasonal cycle
# #    plot_ecmwf_wind_speed(ecmwf_file)                  # Just wind speed

#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd # For pd.to_datetime if time conversion is tricky
# from matplotlib.dates import DateFormatter # Not used
# import scipy.stats as stats # Not used
import seaborn as sns # For color palettes
# import cartopy.crs as ccrs # Not used
# import matplotlib.dates as mdates # Not used for this specific plot
# from matplotlib.colors import LinearSegmentedColormap # Not used
# from matplotlib.ticker import MaxNLocator # Not used

import argparse
import sys
from pathlib import Path

def plot_v10_area_avg_timeseries(ecmwf_file_path_str, nc_variable_name, output_path_str):
    """
    Create an area average timeseries plot for 10m V-Wind Component (V10).
   
    Parameters:
    -----------
    ecmwf_file_path_str : str
        Path to the ECMWF data file.
    nc_variable_name : str
        The NetCDF variable name for V10 (e.g., 'v10').
    output_path_str : str
        Path to save the output plot image.
    """
    try:
        ecmwf_file_path = Path(ecmwf_file_path_str)
        output_path = Path(output_path_str)

        if not ecmwf_file_path.exists():
            print(f"Error: ECMWF data file not found at {ecmwf_file_path}", file=sys.stderr)
            sys.exit(1)

        # Load the dataset
        actual_time_coord_name = None
        try:
            ds = xr.open_dataset(ecmwf_file_path, decode_times=True)
            if 'time' in ds.coords: actual_time_coord_name = 'time'
            elif 'valid_time' in ds.coords: actual_time_coord_name = 'valid_time'
            if not actual_time_coord_name:
                 print(f"Error: Could not identify 'time' or 'valid_time' coordinate after initial load in {ecmwf_file_path}", file=sys.stderr)
                 sys.exit(1)
            if actual_time_coord_name == 'valid_time':
                ds = ds.rename({'valid_time': 'time'})
                actual_time_coord_name = 'time'
        except Exception as e: # Fallback
            print(f"Warning: Initial time decoding failed: {e}. Trying with decode_times=False.", file=sys.stderr)
            try:
                ds = xr.open_dataset(ecmwf_file_path, decode_times=False)
                if 'time' in ds.coords: actual_time_coord_name = 'time'
                elif 'valid_time' in ds.coords: actual_time_coord_name = 'valid_time'
                else: print(f"Error: No time coord in fallback for {ecmwf_file_path}.", file=sys.stderr); sys.exit(1)
                
                time_ds_to_decode = xr.Dataset({actual_time_coord_name: ds[actual_time_coord_name]})
                decoded_time_ds = xr.decode_cf(time_ds_to_decode)
                ds[actual_time_coord_name] = decoded_time_ds[actual_time_coord_name]
                if actual_time_coord_name == 'valid_time':
                    ds = ds.rename({'valid_time': 'time'})
                    actual_time_coord_name = 'time'
            except Exception as e_fallback: print(f"Error: Failed to load {ecmwf_file_path} with fallback: {e_fallback}", file=sys.stderr); sys.exit(1)
       
        if nc_variable_name not in ds:
            print(f"Error: Variable '{nc_variable_name}' (expected for V10) not found in dataset {ecmwf_file_path}", file=sys.stderr)
            print(f"Available variables: {list(ds.data_vars)}", file=sys.stderr)
            sys.exit(1)
            
        data_var = ds[nc_variable_name]

        # Ensure the identified time coordinate is a dimension of the data variable
        dim_to_use_for_time = actual_time_coord_name
        if actual_time_coord_name not in data_var.dims:
            potential_time_dims = [d for d in data_var.dims if actual_time_coord_name in d or 'time' in d.lower()]
            if not potential_time_dims:
                print(f"Error: Time coordinate '{actual_time_coord_name}' is not a dimension of '{nc_variable_name}'. Dims: {data_var.dims}", file=sys.stderr)
                sys.exit(1)
            dim_to_use_for_time = potential_time_dims[0]

        # Determine spatial coordinate names for averaging
        lat_coord_name = 'latitude' if 'latitude' in data_var.dims else 'lat'
        lon_coord_name = 'longitude' if 'longitude' in data_var.dims else 'lon'

        if lat_coord_name not in data_var.dims or lon_coord_name not in data_var.dims:
            print(f"Error: Could not determine latitude/longitude coordinate names for '{nc_variable_name}'. Dims: {data_var.dims}", file=sys.stderr)
            sys.exit(1)
            
        # Calculate area average
        area_mean = data_var.mean(dim=[lat_coord_name, lon_coord_name], skipna=True)
        
        display_units = data_var.attrs.get('units', 'm/s') # Default to m/s for wind components

        plt.figure(figsize=(12, 6))
        
        # Select color for V10 (e.g., from a palette)
        # Original list: ['msl', 'tcc', 'u10', 'v10', 't2m'] -> v10 is index 3
        plot_color = sns.color_palette("viridis", 5)[3] 
        
        area_mean.plot(ax=plt.gca(), linewidth=1.5, color=plot_color)
       
        plt.title(f'10m V-Wind Component (Area Average)', fontsize=14)
        plt.xlabel('Time', fontsize=12)
        plt.ylabel(f'V-Wind ({display_units})', fontsize=12)
       
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        
        try:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved successfully to {output_path}", file=sys.stdout)
        except Exception as e:
            print(f"Error: Failed to save plot to {output_path}: {e}", file=sys.stderr)
            sys.exit(1)
        finally:
            plt.close()
            
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        # import traceback # For debugging
        # traceback.print_exc(file=sys.stderr) # For debugging
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate an area average timeseries plot for 10m V-Wind (V10).")
    parser.add_argument('--files', required=True, nargs=1, help='Path to the input ECMWF NetCDF file.')
    parser.add_argument('--varnames', required=True, nargs=1, help="NetCDF variable name for V10 (e.g., 'v10').")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type (unused by this script).')
    parser.add_argument('--vars', nargs='+', help='Internal variable keys (unused by this script).')

    args = parser.parse_args()
        
    ecmwf_file = args.files[0]
    nc_var_name = args.varnames[0]

    plot_v10_area_avg_timeseries(ecmwf_file, nc_var_name, args.output)