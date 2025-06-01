# #!/usr/bin/env python
# # coding: utf-8

# # In[4]:


# import xarray as xr
# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd
# from matplotlib.dates import DateFormatter
# import scipy.stats as stats


# # In[5]:


# noaa_precip = xr.open_dataset( r"precip_konkan.nc")
# noaa_sst = xr.open_dataset( r"noaa_sst_masked.nc")
# ecmwf = xr.open_dataset( r"ecmwf datas.nc")


# # In[6]:


# # Calculate spatial means for each dataset
# precip_mean = noaa_precip.precip.mean(dim=['lat', 'lon'])
# sst_mean = noaa_sst.sst.mean(dim=['lat', 'lon'])
# ecmwf_means = {
#     'msl': ecmwf.msl.mean(dim=['latitude', 'longitude']),
#     'tcc': ecmwf.tcc.mean(dim=['latitude', 'longitude']),
#     't2m': ecmwf.t2m.mean(dim=['latitude', 'longitude']),
#     'u10': ecmwf.u10.mean(dim=['latitude', 'longitude']),
#     'v10': ecmwf.v10.mean(dim=['latitude', 'longitude'])
# }

# # Create time series DataFrames
# noaa_precip_df = pd.DataFrame({
#     'time': pd.to_datetime(noaa_precip.time.values),
#     'precip': precip_mean.values
# }).set_index('time')

# noaa_sst_df = pd.DataFrame({
#     'time': pd.to_datetime(noaa_sst.time.values),
#     'sst': sst_mean.values
# }).set_index('time')

# ecmwf_df = pd.DataFrame({
#     'time': pd.to_datetime(ecmwf.valid_time.values),
#     'msl': ecmwf_means['msl'].values,
#     'tcc': ecmwf_means['tcc'].values,
#     't2m': ecmwf_means['t2m'].values,
#     'u10': ecmwf_means['u10'].values,
#     'v10': ecmwf_means['v10'].values
# }).set_index('time')

# # Function to plot time series with correlation
# def plot_correlated_series(df1, var1, df2, var2, title=None, filename=None):
#     # Match the time indices
#     common_df = pd.merge(
#         df1[[var1]], df2[[var2]],
#         left_index=True, right_index=True,
#         how='inner'
#     )
   
#     # Calculate correlation
#     correlation, p_value = stats.pearsonr(common_df[var1], common_df[var2])
   
#     # Plot
#     fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={'height_ratios': [3, 1]})
   
#     # Time series plot
#     ax1.plot(common_df.index, common_df[var1], 'b-', label=var1)
#     ax1_twin = ax1.twinx()
#     ax1_twin.plot(common_df.index, common_df[var2], 'r-', label=var2)
   
#     # Add legend
#     lines1, labels1 = ax1.get_legend_handles_labels()
#     lines2, labels2 = ax1_twin.get_legend_handles_labels()
#     ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
   
#     ax1.set_ylabel(var1, color='blue')
#     ax1_twin.set_ylabel(var2, color='red')
   
#     if title:
#         ax1.set_title(f"{title}\nCorrelation: {correlation:.3f} (p-value: {p_value:.3e})")
#     else:
#         ax1.set_title(f"Correlation: {correlation:.3f} (p-value: {p_value:.3e})")
   
#     # Scatter plot
#     ax2.scatter(common_df[var1], common_df[var2], alpha=0.6)
#     ax2.set_xlabel(var1)
#     ax2.set_ylabel(var2)
   
#     # Add trend line
#     z = np.polyfit(common_df[var1], common_df[var2], 1)
#     p = np.poly1d(z)
#     ax2.plot(common_df[var1], p(common_df[var1]), "r--", alpha=0.8)
   
#     plt.tight_layout()
   
#     plt.show()

# # Plot all requested combinations
# # 1. Precipitation and SST
# plot_correlated_series(
#     noaa_precip_df, 'precip',
#     noaa_sst_df, 'sst',
#     title='Combined Analysis: NOAA Precipitation and SST',
#     filename='precip_sst_analysis.png'
# )


# # In[ ]:




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

def load_and_average_variable(file_path_str, nc_var_name, default_time_coord='time', var_role="Variable"):
    """Helper function to load, process time, and spatially average a variable."""
    file_path = Path(file_path_str)
    if not file_path.exists():
        print(f"Error: Data file not found at {file_path} for {var_role}", file=sys.stderr)
        return None

    actual_time_coord_name = None
    try:
        ds = xr.open_dataset(file_path, decode_times=True)
        if default_time_coord in ds.coords: actual_time_coord_name = default_time_coord
        elif 'time' in ds.coords: actual_time_coord_name = 'time'
        # Add other common time coord names if necessary (e.g., 'valid_time')
        
        if not actual_time_coord_name:
            print(f"Error: Could not identify a known time coordinate for {var_role} in {file_path}", file=sys.stderr)
            return None
        if actual_time_coord_name != 'time': # Standardize to 'time'
            ds = ds.rename({actual_time_coord_name: 'time'})
        actual_time_coord_name = 'time'

    except Exception as e:
        print(f"Warning: Initial time decoding failed for {file_path} ({var_role}): {e}. Trying decode_times=False.", file=sys.stderr)
        try:
            ds = xr.open_dataset(file_path, decode_times=False)
            time_coord_fallback = None
            if default_time_coord in ds.coords: time_coord_fallback = default_time_coord
            elif 'time' in ds.coords: time_coord_fallback = 'time'
            elif 'valid_time' in ds.coords: time_coord_fallback = 'valid_time' # Though less common for NOAA precip/sst
            else: print(f"Error: No time coord in fallback for {file_path} ({var_role}).", file=sys.stderr); return None
            
            time_ds_decode = xr.Dataset({time_coord_fallback: ds[time_coord_fallback]})
            decoded_time = xr.decode_cf(time_ds_decode)
            ds[time_coord_fallback] = decoded_time[time_coord_fallback]
            if time_coord_fallback != 'time':
                ds = ds.rename({time_coord_fallback: 'time'})
            actual_time_coord_name = 'time'
        except Exception as e_fb:
            print(f"Error: Failed to load {file_path} ({var_role}) with fallback: {e_fb}", file=sys.stderr)
            return None

    if nc_var_name not in ds:
        print(f"Error: Variable '{nc_var_name}' not found in {file_path} for {var_role}", file=sys.stderr)
        return None
    
    data_var = ds[nc_var_name]

    # For NOAA data, lat/lon are typically 'lat', 'lon'
    lat_dim = next((d for d in ['lat', 'latitude'] if d in data_var.dims), None)
    lon_dim = next((d for d in ['lon', 'longitude'] if d in data_var.dims), None)


    if not lat_dim or not lon_dim:
        print(f"Error: Could not find lat/lon dims for {nc_var_name} in {file_path}. Dims: {data_var.dims}", file=sys.stderr)
        return None
    
    spatial_avg = data_var.mean(dim=[lat_dim, lon_dim], skipna=True)
    
    if not isinstance(ds[actual_time_coord_name].to_index(), pd.DatetimeIndex):
        try:
            ds[actual_time_coord_name] = pd.to_datetime(ds[actual_time_coord_name].values)
        except Exception as e_time_convert:
            print(f"Warning: Could not convert time index to DatetimeIndex for {var_role} from {file_path}: {e_time_convert}", file=sys.stderr)

    return spatial_avg.to_series().rename(nc_var_name)


def plot_sst_precip_combined_analysis(sst_file_path_str, sst_nc_var_name, 
                                      precip_file_path_str, precip_nc_var_name, 
                                      output_path_str):
    """
    Plots combined time series and scatter plot for SST and Precipitation, showing correlation.
    """
    try:
        output_path = Path(output_path_str)

        series_sst = load_and_average_variable(sst_file_path_str, sst_nc_var_name, 
                                               default_time_coord='time', var_role="SST")
        if series_sst is None: sys.exit(1)
        # SST typically in Celsius.

        series_precip = load_and_average_variable(precip_file_path_str, precip_nc_var_name, 
                                                  default_time_coord='time', var_role="Precipitation")
        if series_precip is None: sys.exit(1)
        
        # Get units for labels
        sst_display_units = xr.open_dataset(sst_file_path_str)[sst_nc_var_name].attrs.get('units', '°C')
        precip_display_units = xr.open_dataset(precip_file_path_str)[precip_nc_var_name].attrs.get('units', 'mm/day')


        # Ensure indices are DatetimeIndex before merging/resampling
        if not isinstance(series_sst.index, pd.DatetimeIndex):
            try: series_sst.index = pd.to_datetime(series_sst.index)
            except: print(f"Error: SST series index could not be converted to DatetimeIndex.", file=sys.stderr); sys.exit(1)
        if not isinstance(series_precip.index, pd.DatetimeIndex):
            try: series_precip.index = pd.to_datetime(series_precip.index)
            except: print(f"Error: Precipitation series index could not be converted to DatetimeIndex.", file=sys.stderr); sys.exit(1)
        
        try:
            series_sst_aligned = series_sst.resample('ME').mean()
            series_precip_aligned = series_precip.resample('ME').mean()
            print("SST and Precipitation data resampled to monthly frequency for alignment.", file=sys.stdout)
        except Exception as e_resample:
            print(f"Warning: Could not resample SST/Precipitation to monthly: {e_resample}. Using original time resolution.", file=sys.stderr)
            series_sst_aligned = series_sst
            series_precip_aligned = series_precip
            
        common_df = pd.merge(
            series_sst_aligned.rename('SST'),
            series_precip_aligned.rename('Precip'),
            left_index=True, right_index=True,
            how='inner' # Use 'inner' to keep only overlapping time points
        ).dropna() # Drop rows with NaNs that might result from resampling or initial data

        if common_df.empty or len(common_df) < 2: # Need at least 2 points for correlation
            print("Error: No common overlapping data points found between SST and Precipitation. Cannot compute correlation.", file=sys.stderr)
            fig, ax = plt.subplots(figsize=(12,6))
            ax.text(0.5,0.5, "Insufficient overlapping data for SST-Precipitation analysis", horizontalalignment='center', verticalalignment='center')
            ax.set_title("SST and Precipitation Combined Analysis")
        else:
            correlation, p_value = stats.pearsonr(common_df['SST'], common_df['Precip'])
            
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=False, gridspec_kw={'height_ratios': [3, 1]})
            
            ax1.plot(common_df.index, common_df['SST'], 'b-', label=f'SST ({sst_display_units})')
            ax1.set_ylabel(f'SST ({sst_display_units})', color='blue', fontsize=12)
            ax1.tick_params(axis='y', labelcolor='blue')

            ax1_twin = ax1.twinx()
            ax1_twin.plot(common_df.index, common_df['Precip'], 'r-', label=f'Precip ({precip_display_units})')
            ax1_twin.set_ylabel(f'Precip ({precip_display_units})', color='red', fontsize=12)
            ax1_twin.tick_params(axis='y', labelcolor='red')
            
            ax1.set_xlabel("Time", fontsize=12)
            ax1.set_title(f"SST and Precipitation Time Series\nCorrelation: {correlation:.3f} (p-value: {p_value:.3e})", fontsize=14)
            
            lines, labels = ax1.get_legend_handles_labels()
            lines2, labels2 = ax1_twin.get_legend_handles_labels()
            ax1_twin.legend(lines + lines2, labels + labels2, loc='upper left', fontsize=10)
            ax1.grid(True, linestyle=':', alpha=0.7)

            ax2.scatter(common_df['SST'], common_df['Precip'], alpha=0.6, color='forestgreen') # Changed scatter color
            ax2.set_xlabel(f'SST ({sst_display_units})', fontsize=12)
            ax2.set_ylabel(f'Precip ({precip_display_units})', fontsize=12)
            ax2.set_title('SST vs. Precipitation Scatter Plot with Trend Line', fontsize=12)
            
            z = np.polyfit(common_df['SST'], common_df['Precip'], 1)
            p = np.poly1d(z)
            ax2.plot(common_df['SST'], p(common_df['SST']), "k--", alpha=0.8, linewidth=1.5)
            ax2.grid(True, linestyle=':', alpha=0.7)

        plt.tight_layout(pad=2.0)
        
        try:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved successfully to {output_path}", file=sys.stdout)
        except Exception as e:
            print(f"Error: Failed to save plot to {output_path}: {e}", file=sys.stderr)
            sys.exit(1)
        finally:
            plt.close()

    except Exception as e:
        print(f"An unexpected error occurred in plot_sst_precip_combined_analysis: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a combined analysis plot for SST and Precipitation.")
    parser.add_argument('--files', required=True, nargs='+', 
                        help='Paths to the input NetCDF files (SST and Precipitation sources).')
    parser.add_argument('--varnames', required=True, nargs='+',
                        help="NetCDF variable names. Order should match --vars.")
    parser.add_argument('--vars', required=True, nargs=2, # Expects exactly 2 internal keys
                        help="Two internal variable keys (e.g., 'noaa_sst' 'noaa_precip').")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type (unused by this script).')

    args = parser.parse_args()

    if len(args.varnames) != len(args.vars) or len(args.vars) != 2:
        print("Error: This script expects exactly two --vars (internal keys) and a corresponding number of --varnames.", file=sys.stderr)
        sys.exit(1)
    
    SST_KEY_IDENTIFIER_SUBSTRING = 'sst' 
    PRECIP_KEY_IDENTIFIER_SUBSTRING = 'precip'

    SST_FILE_BASENAME_EXPECTED = "noaa_sst_masked.nc" # From your Streamlit config
    PRECIP_FILE_BASENAME_EXPECTED = "precip_konkan.nc" # From your Streamlit config

    sst_file_arg, precip_file_arg = None, None
    sst_varname_arg, precip_varname_arg = None, None
    
    key1, key2 = args.vars[0], args.vars[1]
    varname1, varname2 = args.varnames[0], args.varnames[1]

    if SST_KEY_IDENTIFIER_SUBSTRING in key1.lower() and PRECIP_KEY_IDENTIFIER_SUBSTRING in key2.lower():
        sst_varname_arg, precip_varname_arg = varname1, varname2
    elif PRECIP_KEY_IDENTIFIER_SUBSTRING in key1.lower() and SST_KEY_IDENTIFIER_SUBSTRING in key2.lower():
        precip_varname_arg, sst_varname_arg = varname1, varname2
    else:
        print(f"Warning: Could not definitively identify SST and Precipitation roles from --vars: {args.vars}. Ensure one key contains '{SST_KEY_IDENTIFIER_SUBSTRING}' and the other '{PRECIP_KEY_IDENTIFIER_SUBSTRING}'.", file=sys.stderr)
        sys.exit(1)

    for f_path_str in args.files: # args.files contains unique paths
        if Path(f_path_str).name == SST_FILE_BASENAME_EXPECTED:
            sst_file_arg = f_path_str
        if Path(f_path_str).name == PRECIP_FILE_BASENAME_EXPECTED:
            precip_file_arg = f_path_str
            
    # Handle cases where one or both might be missing or if they are the same file (less likely for SST/Precip)
    if not sst_file_arg and len(args.files) == 1: # If only one unique file was passed
        if Path(args.files[0]).name == SST_FILE_BASENAME_EXPECTED or SST_FILE_BASENAME_EXPECTED == PRECIP_FILE_BASENAME_EXPECTED:
            sst_file_arg = args.files[0]
    
    if not precip_file_arg and len(args.files) == 1:
        if Path(args.files[0]).name == PRECIP_FILE_BASENAME_EXPECTED or SST_FILE_BASENAME_EXPECTED == PRECIP_FILE_BASENAME_EXPECTED:
            precip_file_arg = args.files[0]

    if not sst_file_arg or not precip_file_arg or not sst_varname_arg or not precip_varname_arg:
        print("Error: Could not determine necessary file paths or variable names for SST and Precipitation.", file=sys.stderr)
        print(f"  SST File: {sst_file_arg}, SST Var: {sst_varname_arg}")
        print(f"  Precip File: {precip_file_arg}, Precip Var: {precip_varname_arg}")
        print(f"  Args.files: {args.files}, Args.vars: {args.vars}, Args.varnames: {args.varnames}")
        sys.exit(1)
    
    print(f"Plotting SST (file: '{sst_file_arg}', var: '{sst_varname_arg}') vs Precip (file: '{precip_file_arg}', var: '{precip_varname_arg}')", file=sys.stdout)

    plot_sst_precip_combined_analysis(sst_file_arg, sst_varname_arg, 
                                      precip_file_arg, precip_nc_var_name=precip_varname_arg, # Ensure correct varname passed
                                      output_path_str=args.output)