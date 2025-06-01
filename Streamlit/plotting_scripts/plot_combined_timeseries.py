#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import argparse
import sys
from pathlib import Path
import traceback
import warnings
from packaging import version

def load_and_average_variable(file_path_str, nc_var_name, expected_time_coord, var_role="Variable"):
    """
    Helper function to load, process time, and spatially average a variable.
    Adjusts for known differences in coordinate names based on var_role.
    """
    file_path = Path(file_path_str)
    if not file_path.exists():
        print(f"Error: Data file not found at {file_path} for {var_role}", file=sys.stderr)
        return None, None # Return None for series and units

    # Determine typical spatial coordinate names based on var_role or file naming conventions
    if "noaa" in file_path.name.lower() or "precip" in nc_var_name.lower() or "sst" in nc_var_name.lower():
        lat_dim_options = ['lat', 'latitude']
        lon_dim_options = ['lon', 'longitude']
    elif "ecmwf" in file_path.name.lower():
        lat_dim_options = ['latitude', 'lat']
        lon_dim_options = ['longitude', 'lon']
    else: # Default
        lat_dim_options = ['latitude', 'lat']
        lon_dim_options = ['longitude', 'lon']

    actual_time_coord_name = None
    try:
        ds = xr.open_dataset(file_path, decode_times=True)
        # Try expected first, then common alternatives
        if expected_time_coord in ds.coords: actual_time_coord_name = expected_time_coord
        elif 'time' in ds.coords: actual_time_coord_name = 'time'
        elif 'valid_time' in ds.coords: actual_time_coord_name = 'valid_time'
        
        if not actual_time_coord_name:
            print(f"Error: Could not identify a known time coordinate for {var_role} in {file_path}", file=sys.stderr)
            return None, None
        if actual_time_coord_name != 'time': # Standardize to 'time'
            ds = ds.rename({actual_time_coord_name: 'time'})
        actual_time_coord_name = 'time'

    except Exception as e:
        print(f"Warning: Initial time decoding failed for {file_path} ({var_role}): {e}. Trying decode_times=False.", file=sys.stderr)
        try:
            ds = xr.open_dataset(file_path, decode_times=False)
            time_coord_fallback = None
            if expected_time_coord in ds.coords: time_coord_fallback = expected_time_coord
            elif 'time' in ds.coords: time_coord_fallback = 'time'
            elif 'valid_time' in ds.coords: time_coord_fallback = 'valid_time'
            else: print(f"Error: No time coord in fallback for {file_path} ({var_role}).", file=sys.stderr); return None, None
            
            time_ds_decode = xr.Dataset({time_coord_fallback: ds[time_coord_fallback]})
            decoded_time = xr.decode_cf(time_ds_decode)
            ds[time_coord_fallback] = decoded_time[time_coord_fallback]
            if time_coord_fallback != 'time':
                ds = ds.rename({time_coord_fallback: 'time'})
            actual_time_coord_name = 'time'
        except Exception as e_fb:
            print(f"Error: Failed to load {file_path} ({var_role}) with fallback: {e_fb}", file=sys.stderr)
            return None, None

    if nc_var_name not in ds:
        print(f"Error: Variable '{nc_var_name}' not found in {file_path} for {var_role}", file=sys.stderr)
        return None, None
    
    data_var = ds[nc_var_name]

    lat_dim = next((d for d in lat_dim_options if d in data_var.dims), None)
    lon_dim = next((d for d in lon_dim_options if d in data_var.dims), None)

    if not lat_dim or not lon_dim:
        print(f"Error: Could not find lat/lon dims for {nc_var_name} in {file_path}. Dims: {data_var.dims}", file=sys.stderr)
        return None, None
    
    spatial_avg = data_var.mean(dim=[lat_dim, lon_dim], skipna=True)
    display_units = data_var.attrs.get('units', 'unknown')
    
    if not isinstance(ds[actual_time_coord_name].to_index(), pd.DatetimeIndex):
        try:
            ds[actual_time_coord_name] = pd.to_datetime(ds[actual_time_coord_name].values)
        except Exception as e_time_convert:
            print(f"Warning: Could not convert time index to DatetimeIndex for {var_role} from {file_path}: {e_time_convert}", file=sys.stderr)

    return spatial_avg.to_series().rename(var_role), display_units


def plot_combined_two_variable_timeseries(
    file1_path_str, var1_nc_name, var1_display_name, var1_expected_time_coord,
    file2_path_str, var2_nc_name, var2_display_name, var2_expected_time_coord,
    output_path_str):
    """
    Plots combined time series and scatter plot for two variables, showing correlation.
    """
    try:
        output_path = Path(output_path_str)

        series1, units1 = load_and_average_variable(file1_path_str, var1_nc_name, var1_expected_time_coord, var1_display_name)
        if series1 is None: sys.exit(1)

        series2, units2 = load_and_average_variable(file2_path_str, var2_nc_name, var2_expected_time_coord, var2_display_name)
        if series2 is None: sys.exit(1)

        # Apply specific unit conversions if needed (example for T2M, MSL, TCC)
        if var1_display_name.lower() == 't2m' and series1.mean() > 100: # Kelvin to Celsius
            series1 = series1 - 273.15
            units1 = '°C'
        elif var1_display_name.lower() == 'msl' and series1.max() > 10000: # Pa to hPa
            series1 = series1 / 100
            units1 = 'hPa'
        elif var1_display_name.lower() == 'tcc' and series1.max() <= 1.0 and series1.min() >= 0.0: # Fraction to %
            series1 = series1 * 100
            units1 = '%'
        
        if var2_display_name.lower() == 't2m' and series2.mean() > 100:
            series2 = series2 - 273.15
            units2 = '°C'
        elif var2_display_name.lower() == 'msl' and series2.max() > 10000:
            series2 = series2 / 100
            units2 = 'hPa'
        elif var2_display_name.lower() == 'tcc' and series2.max() <= 1.0 and series2.min() >= 0.0:
            series2 = series2 * 100
            units2 = '%'
        
        # Ensure indices are DatetimeIndex
        for s_name, series_obj in [(var1_display_name, series1), (var2_display_name, series2)]:
            if not isinstance(series_obj.index, pd.DatetimeIndex):
                try: series_obj.index = pd.to_datetime(series_obj.index)
                except: print(f"Error: {s_name} series index could not be converted to DatetimeIndex.", file=sys.stderr); sys.exit(1)
        
        # Determine resampling frequency based on Pandas version (for PC compatibility)
        pandas_version = pd.__version__
        resample_freq = 'M' 
        if version.parse(pandas_version) < version.parse('2.0.0'):
            try: # Check if pd.Version exists (Pandas >= 1.0.0)
                if pd.Version(pandas_version) < pd.Version('2.0.0'):
                    resample_freq = 'ME'
            except AttributeError: # pd.Version does not exist (Pandas < 1.0.0)
                resample_freq = 'ME' 
        
        try:
            series1_aligned = series1.resample(resample_freq).mean()
            series2_aligned = series2.resample(resample_freq).mean()
            print(f"Data resampled to '{resample_freq}' frequency for alignment.", file=sys.stdout)
        except Exception as e_resample:
            print(f"Warning: Could not resample to '{resample_freq}': {e_resample}. Using original time resolution.", file=sys.stderr)
            series1_aligned = series1
            series2_aligned = series2

        common_df = pd.merge(
            series1_aligned.rename(var1_display_name),
            series2_aligned.rename(var2_display_name),
            left_index=True, right_index=True,
            how='inner' 
        ).dropna()

        if common_df.empty or len(common_df) < 2:
            print(f"Error: No common data points for {var1_display_name} and {var2_display_name}.", file=sys.stderr)
            fig, ax = plt.subplots(figsize=(12,6))
            ax.text(0.5,0.5, f"Insufficient overlapping data for {var1_display_name}-{var2_display_name} analysis", ha='center', va='center')
            ax.set_title(f"{var1_display_name} and {var2_display_name} Combined Analysis")
        else:
            correlation, p_value = stats.pearsonr(common_df[var1_display_name], common_df[var2_display_name])
            
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=False, gridspec_kw={'height_ratios': [3, 1]})
            
            ax1.plot(common_df.index, common_df[var1_display_name], 'b-', label=f'{var1_display_name} ({units1})')
            ax1.set_ylabel(f'{var1_display_name} ({units1})', color='blue', fontsize=12)
            ax1.tick_params(axis='y', labelcolor='blue')

            ax1_twin = ax1.twinx()
            ax1_twin.plot(common_df.index, common_df[var2_display_name], 'r-', label=f'{var2_display_name} ({units2})')
            ax1_twin.set_ylabel(f'{var2_display_name} ({units2})', color='red', fontsize=12)
            ax1_twin.tick_params(axis='y', labelcolor='red')
            
            ax1.set_xlabel("Time (Monthly Mean)", fontsize=12)
            ax1.set_title(f"{var1_display_name} and {var2_display_name} Time Series\nCorrelation: {correlation:.3f} (p-value: {p_value:.3e})", fontsize=14)
            
            lines, labels = ax1.get_legend_handles_labels()
            lines2, labels2 = ax1_twin.get_legend_handles_labels()
            ax1_twin.legend(lines + lines2, labels + labels2, loc='upper left', fontsize=10)
            ax1.grid(True, linestyle=':', alpha=0.7)

            ax2.scatter(common_df[var1_display_name], common_df[var2_display_name], alpha=0.6)
            ax2.set_xlabel(f'{var1_display_name} ({units1})', fontsize=12)
            ax2.set_ylabel(f'{var2_display_name} ({units2})', fontsize=12)
            ax2.set_title(f'{var1_display_name} vs. {var2_display_name} Scatter Plot with Trend Line', fontsize=12)
            
            z = np.polyfit(common_df[var1_display_name], common_df[var2_display_name], 1)
            p = np.poly1d(z)
            ax2.plot(common_df[var1_display_name], p(common_df[var1_display_name]), "k--", alpha=0.8, linewidth=1.5)
            ax2.grid(True, linestyle=':', alpha=0.7)

        plt.tight_layout(pad=2.0)
        
        try:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved successfully to {output_path}", file=sys.stdout)
        except Exception as e_save:
            print(f"Error: Failed to save plot to {output_path}: {e_save}", file=sys.stderr)
            sys.exit(1)
        finally:
            plt.close()

    except Exception as e_main:
        print(f"An unexpected error occurred in plot_combined_two_variable_timeseries: {e_main}", file=sys.stderr)
        sys.exit(1)
        
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate a combined analysis plot for two variables.")
    parser.add_argument('--files', required=True, nargs='+', 
                        help='Paths to one or two NetCDF files (if one file has both vars, pass it once).')
    parser.add_argument('--varnames', required=True, nargs=2, 
                        help="NetCDF variable names for var1 and var2 (var1_nc_name var2_nc_name).")
    parser.add_argument('--displaynames', required=True, nargs=2,
                        help="Display names for var1 and var2, used in plot titles/labels (var1_disp_name var2_disp_name).")
    parser.add_argument('--timecoords', required=True, nargs=2,
                        help="Expected time coordinate names for var1 and var2 (e.g., 'time' 'valid_time').")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type (e.g., "Combined Timeseries: Var1 and Var2").')
    parser.add_argument('--vars', nargs='*', help='Internal variable keys (optional, for context).')

    try:
        args = parser.parse_args()

        if len(args.files) == 1:
            file1 = args.files[0]
            file2 = args.files[0]
        elif len(args.files) == 2:
            file1, file2 = args.files
        else:
            raise ValueError("You must provide one or two files only.")

        plot_combined_two_variable_timeseries(
            file1, args.varnames[0], args.displaynames[0], args.timecoords[0],
            file2, args.varnames[1], args.displaynames[1], args.timecoords[1],
            args.output
        )

    except Exception as e:
        print("❌ Exception occurred:", e)
        traceback.print_exc()
