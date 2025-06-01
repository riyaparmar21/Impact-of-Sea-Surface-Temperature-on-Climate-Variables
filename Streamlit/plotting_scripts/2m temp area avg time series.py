#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# from matplotlib.dates import DateFormatter # Not explicitly used in the relevant function
# import scipy.stats as stats # Not explicitly used in the relevant function
import seaborn as sns
# import cartopy.crs as ccrs # Not used for this plot type
# import matplotlib.dates as mdates # Not explicitly used in the relevant function
# from matplotlib.colors import LinearSegmentedColormap # Not used
# from matplotlib.ticker import MaxNLocator # Not used

import argparse
import sys
from pathlib import Path

def plot_2m_temp_area_avg_timeseries(ecmwf_file_path_str, nc_variable_name, output_path_str):
    """
    Create an area average timeseries plot for 2m temperature.
   
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
        # Use decode_times=False for robustness, then decode manually or let xarray handle it on access
        try:
            ds = xr.open_dataset(ecmwf_file_path, decode_times=True) 
        except Exception as e:
            # Fallback if default time decoding fails
            print(f"Warning: Initial time decoding failed: {e}. Trying with decode_times=False.", file=sys.stderr)
            try:
                ds = xr.open_dataset(ecmwf_file_path, decode_times=False)
                # Attempt to decode/use common time coordinates
                if 'time' in ds.coords:
                    ds['time'] = xr.decode_cf(ds[['time']])['time']
                elif 'valid_time' in ds.coords: # Common in some ECMWF GRIB conversions
                    ds = ds.rename({'valid_time': 'time'})
                    ds['time'] = xr.decode_cf(ds[['time']])['time']
                else: # Add other potential time coord names if necessary
                    print(f"Error: Could not find or decode a recognizable time coordinate in {ecmwf_file_path}", file=sys.stderr)
                    sys.exit(1)
            except Exception as e_fallback:
                print(f"Error: Failed to load dataset {ecmwf_file_path} even with fallback: {e_fallback}", file=sys.stderr)
                sys.exit(1)
       
        if nc_variable_name not in ds:
            print(f"Error: Variable '{nc_variable_name}' (expected for 2m Temp) not found in dataset {ecmwf_file_path}", file=sys.stderr)
            print(f"Available variables: {list(ds.data_vars)}", file=sys.stderr)
            sys.exit(1)
           
        data_var = ds[nc_variable_name]

        # Determine spatial coordinate names (latitude, longitude might vary)
        lat_coord_name = None
        lon_coord_name = None
        if 'latitude' in data_var.dims:
            lat_coord_name = 'latitude'
        elif 'lat' in data_var.dims:
            lat_coord_name = 'lat'
        
        if 'longitude' in data_var.dims:
            lon_coord_name = 'longitude'
        elif 'lon' in data_var.dims:
            lon_coord_name = 'lon'

        if not lat_coord_name or not lon_coord_name:
            print(f"Error: Could not determine latitude/longitude coordinate names in variable '{nc_variable_name}'. Found dims: {data_var.dims}", file=sys.stderr)
            sys.exit(1)

        # Create figure for area average timeseries
        plt.figure(figsize=(12, 6))
       
        area_mean = data_var.mean(dim=[lat_coord_name, lon_coord_name])
        current_units = data_var.attrs.get('units', 'Unknown')
        display_units = 'K' # Default to Kelvin

        # For temperature, convert to Celsius if it appears to be in Kelvin
        # Using a fixed threshold, adjust if necessary for your data's typical range
        if area_mean.mean(skipna=True).item() > 100:  # Check if mean is high, suggesting Kelvin
            area_mean = area_mean - 273.15  # Convert to Celsius
            display_units = '°C'
            print(f"Data for '{nc_variable_name}' appears to be in Kelvin, converted to Celsius.", file=sys.stdout)
        else:
            display_units = current_units if current_units != 'Unknown' else 'K'
            print(f"Data for '{nc_variable_name}' assumed to be in {display_units} (or already Celsius).", file=sys.stdout)
       
        # Use a specific color for t2m or a default
        plot_color = sns.color_palette("viridis", 5)[4] # Index 4 for 't2m' in original list
        area_mean.plot(linewidth=1.5, color=plot_color)
       
        plt.title(f'2m Temperature (Area Average)', fontsize=14) # Title is specific to this script
        plt.xlabel('Time', fontsize=12)
        plt.ylabel(f'Temperature ({display_units})', fontsize=12)
       
        plt.grid(True)
        plt.tight_layout() # Adjust layout to prevent labels from overlapping
        
        try:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved successfully to {output_path}", file=sys.stdout)
        except Exception as e:
            print(f"Error: Failed to save plot to {output_path}: {e}", file=sys.stderr)
            sys.exit(1)
        finally:
            plt.close() # Close the plot to free memory
   
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate an area average timeseries plot for 2m temperature.")
    parser.add_argument('--files', required=True, nargs='+', help='Path to the input ECMWF NetCDF file(s). Expects one file.')
    parser.add_argument('--varnames', required=True, nargs='+', help="NetCDF variable name(s) for 2m temperature. Expects one variable name (e.g., 't2m').")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type being generated (unused by this script but part of interface).')
    parser.add_argument('--vars', nargs='+', help='Internal variable keys (unused by this script but part of interface).')

    args = parser.parse_args()

    if len(args.files) != 1:
        print("Error: This script expects exactly one input file path for --files.", file=sys.stderr)
        sys.exit(1)
    if len(args.varnames) != 1:
        print("Error: This script expects exactly one variable name for --varnames.", file=sys.stderr)
        sys.exit(1)
        
    ecmwf_file = args.files[0]
    nc_var_name = args.varnames[0] # This should be 't2m' or equivalent as per Streamlit app config

    # Optional: Add a check if nc_var_name is indeed what's expected for 2m temp, e.g. 't2m'
    # if nc_var_name.lower() not in ['t2m', '2t']: # Common names for 2m temperature
    #     print(f"Warning: Received variable name '{nc_var_name}'. Expected something like 't2m'. Proceeding anyway.", file=sys.stderr)

    plot_2m_temp_area_avg_timeseries(ecmwf_file, nc_var_name, args.output)