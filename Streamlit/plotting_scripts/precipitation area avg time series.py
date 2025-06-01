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


# # In[3]:


# def plot_noaa_precip_area_average(precip_file):
#     """
#     Create timeseries plot for NOAA Precipitation data - Area average only
#     """
#     # Load the dataset
#     ds = xr.open_dataset(precip_file)
   
#     # Create a figure for area-averaged precipitation
#     plt.figure(figsize=(12, 5))
   
#     # Plot: Area-averaged precipitation timeseries
#     precip_mean = ds.precip.mean(dim=['lat', 'lon'])
#     precip_mean.plot(color='blue', linewidth=1.5)
   
#     plt.title('NOAA Area-averaged Precipitation (2019-2025)', fontsize=14)
#     plt.xlabel('Time', fontsize=12)
#     plt.ylabel('Precipitation (mm/day)', fontsize=12)
#     plt.gca().xaxis.set_major_locator(mdates.YearLocator())
#     plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
   
#     plt.tight_layout()
#     #plt.savefig('noaa_precip_area_average.png', dpi=300)
#     plt.show()
#     plt.close()


# # In[5]:


# if __name__ == "__main__":
#     # Define your file paths here
#     noaa_precip_file = r"precip_konkan.nc"

#     # Generate individual dataset plots
#     plot_noaa_precip_area_average(noaa_precip_file)


# # In[ ]:




#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd # For pd.to_datetime if time conversion is tricky
# from matplotlib.dates import DateFormatter # mdates.DateFormatter is used
# import scipy.stats as stats # Not used
# import seaborn as sns # Not directly used for this plot
# import cartopy.crs as ccrs # Not used
import matplotlib.dates as mdates # For x-axis formatting
# from matplotlib.colors import LinearSegmentedColormap # Not used
# from matplotlib.ticker import MaxNLocator # Not used

import argparse
import sys
from pathlib import Path

def plot_precip_area_avg_timeseries(noaa_precip_file_path_str, nc_variable_name, output_path_str):
    """
    Create an area average timeseries plot for NOAA Precipitation.
   
    Parameters:
    -----------
    noaa_precip_file_path_str : str
        Path to the NOAA precipitation data file.
    nc_variable_name : str
        The NetCDF variable name for precipitation (e.g., 'precip').
    output_path_str : str
        Path to save the output plot image.
    """
    try:
        noaa_precip_file_path = Path(noaa_precip_file_path_str)
        output_path = Path(output_path_str)

        if not noaa_precip_file_path.exists():
            print(f"Error: NOAA precipitation data file not found at {noaa_precip_file_path}", file=sys.stderr)
            sys.exit(1)

        # Load the dataset
        actual_time_coord_name = None
        try:
            ds_precip = xr.open_dataset(noaa_precip_file_path, decode_times=True)
            if 'time' in ds_precip.coords: actual_time_coord_name = 'time'
            if not actual_time_coord_name:
                 print(f"Error: Could not identify 'time' coordinate after initial load in {noaa_precip_file_path}", file=sys.stderr)
                 sys.exit(1)
        except Exception as e: # Fallback
            print(f"Warning: Initial time decoding failed for {noaa_precip_file_path}: {e}. Trying with decode_times=False.", file=sys.stderr)
            try:
                ds_precip = xr.open_dataset(noaa_precip_file_path, decode_times=False)
                if 'time' in ds_precip.coords: actual_time_coord_name = 'time'
                else: print(f"Error: No 'time' coord in fallback for {noaa_precip_file_path}.", file=sys.stderr); sys.exit(1)
                
                time_ds_to_decode = xr.Dataset({actual_time_coord_name: ds_precip[actual_time_coord_name]})
                decoded_time_ds = xr.decode_cf(time_ds_to_decode)
                ds_precip[actual_time_coord_name] = decoded_time_ds[actual_time_coord_name]
            except Exception as e_fallback: print(f"Error: Failed to load {noaa_precip_file_path} with fallback: {e_fallback}", file=sys.stderr); sys.exit(1)

        if nc_variable_name not in ds_precip:
            print(f"Error: Variable '{nc_variable_name}' (expected for precipitation) not found in dataset {noaa_precip_file_path}", file=sys.stderr)
            print(f"Available variables: {list(ds_precip.data_vars)}", file=sys.stderr)
            sys.exit(1)
        
        data_var = ds_precip[nc_variable_name]
        
        # Ensure the identified time coordinate is a dimension of the data variable
        dim_to_use_for_time = actual_time_coord_name
        if actual_time_coord_name not in data_var.dims:
            potential_time_dims = [d for d in data_var.dims if actual_time_coord_name in d or 'time' in d.lower()]
            if not potential_time_dims:
                print(f"Error: Time coordinate '{actual_time_coord_name}' is not a dimension of '{nc_variable_name}'. Dims: {data_var.dims}", file=sys.stderr)
                sys.exit(1)
            dim_to_use_for_time = potential_time_dims[0]
            # print(f"Using dimension '{dim_to_use_for_time}' for time operations on variable '{nc_variable_name}'.", file=sys.stdout) # Can be verbose

        # Determine spatial coordinate names (typically 'lat', 'lon' for NOAA CPC data)
        lat_coord_name = 'lat' if 'lat' in data_var.dims else 'latitude'
        lon_coord_name = 'lon' if 'lon' in data_var.dims else 'longitude'

        if lat_coord_name not in data_var.dims or lon_coord_name not in data_var.dims:
            print(f"Error: Could not determine lat/lon dimensions for '{nc_variable_name}'. Dims: {data_var.dims}", file=sys.stderr)
            sys.exit(1)

        # Calculate area average
        precip_mean = data_var.mean(dim=[lat_coord_name, lon_coord_name], skipna=True)
        
        display_units = data_var.attrs.get('units', 'mm/day') # Default to mm/day
        plot_title = 'Area-Averaged Precipitation' # Generic title

        plt.figure(figsize=(12, 5))
        precip_mean.plot(ax=plt.gca(), color='blue', linewidth=1.5)
       
        plt.title(plot_title, fontsize=14)
        plt.xlabel('Time', fontsize=12)
        plt.ylabel(f'Precipitation ({display_units})', fontsize=12)
        
        # Format x-axis to show years
        plt.gca().xaxis.set_major_locator(mdates.YearLocator())
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        plt.gca().tick_params(axis='x', rotation=45) # Rotate for better readability if many years
       
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
    parser = argparse.ArgumentParser(description="Generate an area average timeseries plot for NOAA precipitation.")
    parser.add_argument('--files', required=True, nargs=1, help='Path to the input NOAA precipitation NetCDF file.')
    parser.add_argument('--varnames', required=True, nargs=1, help="NetCDF variable name for precipitation (e.g., 'precip').")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type (unused by this script).')
    parser.add_argument('--vars', nargs='+', help='Internal variable keys (unused by this script).')
    
    args = parser.parse_args()
        
    noaa_file = args.files[0]
    nc_var_name = args.varnames[0]

    plot_precip_area_avg_timeseries(noaa_file, nc_var_name, args.output)