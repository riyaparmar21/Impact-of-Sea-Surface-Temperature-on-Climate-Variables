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


# # In[11]:


# #############################################
# # NOAA SST Dataset Plotting Code
# #############################################

# def plot_noaa_sst(sst_file):
#     """
#     Create timeseries plots for NOAA Sea Surface Temperature data
#     """
#     # Load the dataset
#     ds = xr.open_dataset(sst_file)
   
#     # Create area-averaged SST timeseries
#     plt.figure(figsize=(12, 6))
#     sst_mean = ds.sst.mean(dim=['lat', 'lon'])
   
#     # Convert to Celsius if it appears to be in Kelvin
#     if sst_mean.mean() > 100:  # Likely in Kelvin
#         sst_mean = sst_mean - 273.15
#         units = '°C'
#     else:
#         units = '°C'  # Assuming it's already in Celsius
   
#     # Create anomaly plot (difference from climatological mean)
#     climatology = ds.sst.groupby('time.month').mean()
#     anomalies = ds.sst.groupby('time.month') - climatology
#     anom_mean = anomalies.mean(dim=['lat', 'lon'])
   
#     plt.figure(figsize=(12, 6))
#     anom_mean.plot(linewidth=1.5, color='purple')
   
#     plt.title('NOAA SST Anomalies (Difference from Climatological Mean)', fontsize=14)
#     plt.xlabel('Time', fontsize=12)
#     plt.ylabel(f'Temperature Anomaly ({units})', fontsize=12)
#     plt.grid(True)
#     plt.tight_layout()
#     plt.show()
#     #plt.savefig('noaa_sst_anomalies.png', dpi=300)
#     plt.close()


# # In[12]:


# if __name__ == "__main__":

#     noaa_sst_file =  r"noaa_sst_masked.nc"
    
#     plot_noaa_sst(noaa_sst_file)


# # In[ ]:



#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd # For pd.to_datetime if time conversion is tricky
# from matplotlib.dates import DateFormatter # Not used
# import scipy.stats as stats # Not used
# import seaborn as sns # Not directly used for this plot
# import cartopy.crs as ccrs # Not used
# import matplotlib.dates as mdates # Not used
# from matplotlib.colors import LinearSegmentedColormap # Not used
# from matplotlib.ticker import MaxNLocator # Not used

import argparse
import sys
from pathlib import Path

def plot_sst_anomalies_timeseries(noaa_sst_file_path_str, nc_variable_name, output_path_str):
    """
    Create a timeseries plot for area-averaged NOAA SST anomalies.
    Anomalies are calculated as the difference from the monthly climatological mean.
   
    Parameters:
    -----------
    noaa_sst_file_path_str : str
        Path to the NOAA SST data file.
    nc_variable_name : str
        The NetCDF variable name for SST (e.g., 'sst').
    output_path_str : str
        Path to save the output plot image.
    """
    try:
        noaa_sst_file_path = Path(noaa_sst_file_path_str)
        output_path = Path(output_path_str)

        if not noaa_sst_file_path.exists():
            print(f"Error: NOAA SST data file not found at {noaa_sst_file_path}", file=sys.stderr)
            sys.exit(1)

        # Load the dataset
        actual_time_coord_name = None
        try:
            ds_sst = xr.open_dataset(noaa_sst_file_path, decode_times=True)
            if 'time' in ds_sst.coords: actual_time_coord_name = 'time'
            if not actual_time_coord_name:
                 print(f"Error: Could not identify 'time' coordinate after initial load in {noaa_sst_file_path}", file=sys.stderr)
                 sys.exit(1)
        except Exception as e: # Fallback
            print(f"Warning: Initial time decoding failed for {noaa_sst_file_path}: {e}. Trying with decode_times=False.", file=sys.stderr)
            try:
                ds_sst = xr.open_dataset(noaa_sst_file_path, decode_times=False)
                if 'time' in ds_sst.coords: actual_time_coord_name = 'time'
                else: print(f"Error: No 'time' coord in fallback for {noaa_sst_file_path}.", file=sys.stderr); sys.exit(1)
                
                time_ds_to_decode = xr.Dataset({actual_time_coord_name: ds_sst[actual_time_coord_name]})
                decoded_time_ds = xr.decode_cf(time_ds_to_decode)
                ds_sst[actual_time_coord_name] = decoded_time_ds[actual_time_coord_name]
            except Exception as e_fallback: print(f"Error: Failed to load {noaa_sst_file_path} with fallback: {e_fallback}", file=sys.stderr); sys.exit(1)

        if nc_variable_name not in ds_sst:
            print(f"Error: Variable '{nc_variable_name}' (expected for SST) not found in dataset {noaa_sst_file_path}", file=sys.stderr)
            print(f"Available variables: {list(ds_sst.data_vars)}", file=sys.stderr)
            sys.exit(1)
        
        sst_data_var = ds_sst[nc_variable_name]
        
        # Ensure the identified time coordinate is a dimension of the data variable
        dim_to_use_for_time = actual_time_coord_name
        if actual_time_coord_name not in sst_data_var.dims:
            potential_time_dims = [d for d in sst_data_var.dims if actual_time_coord_name in d or 'time' in d.lower()]
            if not potential_time_dims:
                print(f"Error: Time coordinate '{actual_time_coord_name}' is not a dimension of '{nc_variable_name}'. Dims: {sst_data_var.dims}", file=sys.stderr)
                sys.exit(1)
            dim_to_use_for_time = potential_time_dims[0]
            # print(f"Using dimension '{dim_to_use_for_time}' for time operations on variable '{nc_variable_name}'.", file=sys.stdout)

        # Determine spatial coordinate names (typically 'lat', 'lon')
        lat_coord_name = 'lat' if 'lat' in sst_data_var.dims else 'latitude'
        lon_coord_name = 'lon' if 'lon' in sst_data_var.dims else 'longitude'

        if lat_coord_name not in sst_data_var.dims or lon_coord_name not in sst_data_var.dims:
            print(f"Error: Could not determine lat/lon dimensions for '{nc_variable_name}'. Dims: {sst_data_var.dims}", file=sys.stderr)
            sys.exit(1)

        # Handle potential unit conversion for SST data before anomaly calculation
        # Anomalies should be in the same units as the base temperature (e.g., °C or K differences)
        units_for_anomaly = "°C" # Default target for anomaly display
        sst_processed_for_anomaly = sst_data_var.copy() # Work on a copy

        if sst_data_var.mean(skipna=True).item() > 100:  # Heuristic: if mean SST > 100, likely Kelvin
            sst_processed_for_anomaly = sst_data_var - 273.15
            print(f"SST data for anomaly calculation appears to be in Kelvin, converted to Celsius.", file=sys.stdout)
        else: # Assumed to be Celsius or already converted
            units_original = sst_data_var.attrs.get('units', '°C')
            if 'k' in units_original.lower(): # If units attribute explicitly says K but values are low
                units_for_anomaly = "K"
            print(f"SST data for anomaly calculation assumed to be in Celsius or original units ({units_for_anomaly}).", file=sys.stdout)

        # Calculate monthly climatology
        climatology = sst_processed_for_anomaly.groupby(f"{actual_time_coord_name}.month").mean(dim=dim_to_use_for_time, skipna=True)
        
        # Calculate anomalies
        anomalies = sst_processed_for_anomaly.groupby(f"{actual_time_coord_name}.month") - climatology
        
        # Calculate area-average of anomalies
        anom_mean = anomalies.mean(dim=[lat_coord_name, lon_coord_name], skipna=True)
       
        plt.figure(figsize=(12, 6))
        anom_mean.plot(linewidth=1.5, color='purple')
       
        plt.title('SST Anomalies (Difference from Monthly Climatology)', fontsize=14)
        plt.xlabel('Time', fontsize=12)
        plt.ylabel(f'Temperature Anomaly ({units_for_anomaly})', fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.axhline(0, color='black', linestyle='--', linewidth=0.8) # Add a zero line
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
    parser = argparse.ArgumentParser(description="Generate a timeseries plot for NOAA SST anomalies.")
    parser.add_argument('--files', required=True, nargs=1, help='Path to the input NOAA SST NetCDF file.')
    parser.add_argument('--varnames', required=True, nargs=1, help="NetCDF variable name for SST (e.g., 'sst').")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type (unused by this script).')
    parser.add_argument('--vars', nargs='+', help='Internal variable keys (unused by this script).')
    
    args = parser.parse_args()
        
    noaa_file = args.files[0]
    nc_var_name = args.varnames[0]

    plot_sst_anomalies_timeseries(noaa_file, nc_var_name, args.output)
