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


# # In[94]:


# def plot_ecmwf_vars(file_path):
#     ds = xr.open_dataset(file_path)
#     # For Mean Sea Level Pressure (msl)
#     msl_mean = ds.msl.mean(dim=['latitude', 'longitude'])
   
#     # Convert time to a more readable format for x-axis
#     times = pd.to_datetime(ds.valid_time.values)
#     years = [t.strftime('%Y') for t in times]
   
#     # Get unique years for bars
#     unique_years = sorted(set(years))
#     yearly_msl = {year: [] for year in unique_years}
    
#     # For Temperature at 2m (t2m)
#     # Convert Kelvin to Celsius
#     t2m_mean = ds.t2m.mean(dim=['latitude', 'longitude']) - 273.15
   
#     # Plot temperature as monthly bars
#     months = [t.strftime('%Y-%m') for t in times]
   
#     # Select just the first 24 months for better readability
#     plt.figure(figsize=(14,9))
#     plt.bar(months[:60], t2m_mean[:60], color='red')
#     plt.xlabel('Month')
#     plt.ylabel('Mean Temperature (°C)')
#     plt.title('Monthly Mean Temperature')
#     plt.xticks(rotation=90)
#     plt.grid(axis='y', linestyle='--', alpha=0.7)
#     plt.tight_layout()
#     plt.show()
#     #plt.savefig('ecmwf_t2m_monthly.png')
#     plt.close()


# # In[95]:


# plot_ecmwf_vars(r"ecmwf datas.nc")


# # In[ ]:




#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import sys
from pathlib import Path
import traceback

def plot_2m_temp_monthly_mean_bar(ecmwf_file_path_str, nc_variable_name, output_path_str):
    """
    Plots the monthly mean 2m temperature as a bar chart.
   
    Parameters:
    -----------
    ecmwf_file_path_str : str
        Path to the ECMWF data file.
    nc_variable_name : str
        The NetCDF variable name for 2m temperature (e.g., 't2m').
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
        try:
            ds = xr.open_dataset(ecmwf_file_path, decode_times=True)
            print("Dataset loaded with decode_times=True", file=sys.stdout)
        except Exception as e:
            print(f"Warning: Initial time decoding failed: {e}. Trying with decode_times=False.", file=sys.stderr)
            try:
                ds = xr.open_dataset(ecmwf_file_path, decode_times=False)
                # Find time coordinate
                time_coord_name = None
                for coord in ds.coords:
                    if 'time' in coord.lower() or 'date' in coord.lower():
                        time_coord_name = coord
                        break
                if not time_coord_name:
                    print(f"Error: No recognizable time coordinate found. Available coordinates: {list(ds.coords)}", file=sys.stderr)
                    sys.exit(1)

                print(f"Found time coordinate: {time_coord_name}", file=sys.stdout)
                # Try to decode the time coordinate
                try:
                    time_ds = xr.Dataset({time_coord_name: ds[time_coord_name]})
                    decoded_time_ds = xr.decode_cf(time_ds)
                    ds[time_coord_name] = decoded_time_ds[time_coord_name]
                except Exception as e_decode:
                    print(f"Error: Failed to decode time coordinate {time_coord_name}: {e_decode}", file=sys.stderr)
                    sys.exit(1)
            except Exception as e_fallback:
                print(f"Error: Failed to load dataset {ecmwf_file_path} even with fallback: {e_fallback}", file=sys.stderr)
                sys.exit(1)

        if nc_variable_name not in ds:
            print(f"Error: Variable '{nc_variable_name}' (expected for 2m Temp) not found in dataset {ecmwf_file_path}", file=sys.stderr)
            print(f"Available variables: {list(ds.data_vars)}", file=sys.stderr)
            sys.exit(1)
            
        data_var = ds[nc_variable_name]

        # Determine the time dimension for t2m
        time_dim_name = None
        for dim in data_var.dims:
            if 'time' in dim.lower() or 'date' in dim.lower():
                time_dim_name = dim
                break
        if not time_dim_name:
            print(f"Error: No time dimension found in variable '{nc_variable_name}'. Dims: {data_var.dims}", file=sys.stderr)
            sys.exit(1)

        # Rename the time dimension to 'time' if necessary
        if time_dim_name != 'time':
            print(f"Renaming time dimension '{time_dim_name}' to 'time' for variable '{nc_variable_name}'", file=sys.stdout)
            ds = ds.rename({time_dim_name: 'time'})
            time_dim_name = 'time'

        # Verify the renaming
        if time_dim_name not in data_var.dims:
            data_var = ds[nc_variable_name]  # Refresh data_var after renaming
            if time_dim_name not in data_var.dims:
                print(f"Error: Time dimension '{time_dim_name}' still not found in variable '{nc_variable_name}' after renaming. Dims: {data_var.dims}", file=sys.stderr)
                sys.exit(1)

        # Determine spatial coordinate names
        lat_coord_name = 'latitude' if 'latitude' in data_var.dims else 'lat' if 'lat' in data_var.dims else None
        lon_coord_name = 'longitude' if 'longitude' in data_var.dims else 'lon' if 'lon' in data_var.dims else None
        if not lat_coord_name or not lon_coord_name:
            print(f"Error: Could not determine lat/lon coordinate names. Dims: {data_var.dims}", file=sys.stderr)
            sys.exit(1)

        # Calculate area average
        t2m_area_mean = data_var.mean(dim=[lat_coord_name, lon_coord_name], skipna=True)
        
        # Unit conversion: Convert Kelvin to Celsius if applicable
        display_units = 'K'
        mean_value = t2m_area_mean.mean(skipna=True).item()
        if mean_value > 100:  # Heuristic for Kelvin
            t2m_area_mean = t2m_area_mean - 273.15
            display_units = '°C'
            print(f"Data for '{nc_variable_name}' appears to be in Kelvin, converted to Celsius.", file=sys.stdout)
        else:
            original_units = data_var.attrs.get('units', 'K')
            display_units = original_units
            print(f"Data for '{nc_variable_name}' presented in original units: {display_units} (or already Celsius-like values).", file=sys.stdout)

        # Ensure time is a DatetimeIndex
        if not isinstance(t2m_area_mean[time_dim_name].to_index(), pd.DatetimeIndex):
            try:
                t2m_area_mean[time_dim_name] = pd.to_datetime(t2m_area_mean[time_dim_name].values)
            except Exception as e_time_convert:
                print(f"Error: Could not convert time coordinate to pandas DatetimeIndex: {e_time_convert}", file=sys.stderr)
                sys.exit(1)

        # Resample to monthly mean
        try:
            monthly_mean_temp = t2m_area_mean.resample({time_dim_name: 'ME'}).mean(skipna=True)  # 'ME' for month-end
        except Exception as e_resample:
            print(f"Error: Failed to resample data to monthly mean: {e_resample}", file=sys.stderr)
            sys.exit(1)
        
        # Prepare month labels for x-axis
        month_labels = monthly_mean_temp[time_dim_name].dt.strftime('%Y-%m').values
        
        # Plot temperature as monthly bars
        plt.figure(figsize=(14, 9))
        num_months_to_plot = len(month_labels)
        max_bars = 72
        if num_months_to_plot > max_bars:
            print(f"Warning: Data has {num_months_to_plot} months. Plotting the first {max_bars} for readability.", file=sys.stdout)
            month_labels_plot = month_labels[:max_bars]
            plot_values = monthly_mean_temp.data[:max_bars]
        else:
            month_labels_plot = month_labels
            plot_values = monthly_mean_temp.data

        plt.bar(month_labels_plot, plot_values, color='skyblue')
        plt.xlabel('Month (YYYY-MM)', fontsize=12)
        plt.ylabel(f'Mean Temperature ({display_units})', fontsize=12)
        plt.title('Monthly Mean 2m Temperature (Area Average)', fontsize=14)
        plt.xticks(rotation=90, fontsize=10)
        plt.yticks(fontsize=10)
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()
        
        try:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved successfully to {output_path}", file=sys.stdout)
        except Exception as e_save:
            print(f"Error: Failed to save plot to {output_path}: {e_save}", file=sys.stderr)
            sys.exit(1)
        finally:
            plt.close()

    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        print(traceback.format_exc(), file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a monthly mean bar plot for 2m temperature.")
    parser.add_argument('--files', required=True, nargs='+', help='Path to the input ECMWF NetCDF file(s). Expects one file.')
    parser.add_argument('--varnames', required=True, nargs='+', help="NetCDF variable name(s) for 2m temperature. Expects one (e.g., 't2m').")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type (unused by this script).')
    parser.add_argument('--vars', nargs='+', help='Internal variable keys (unused by this script).')

    args = parser.parse_args()

    if len(args.files) != 1:
        print("Error: This script expects exactly one input file path for --files.", file=sys.stderr)
        sys.exit(1)
    if len(args.varnames) != 1:
        print("Error: This script expects exactly one variable name for --varnames.", file=sys.stderr)
        sys.exit(1)
        
    ecmwf_file = args.files[0]
    nc_var_name = args.varnames[0]

    plot_2m_temp_monthly_mean_bar(ecmwf_file, nc_var_name, args.output)