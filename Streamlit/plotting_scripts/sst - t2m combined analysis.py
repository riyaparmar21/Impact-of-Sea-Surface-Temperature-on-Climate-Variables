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


# # In[7]:


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


# plot_correlated_series(
#     noaa_sst_df, 'sst',
#     ecmwf_df, 't2m',
#     title='Combined Analysis: SST and 2m Temperature',
#     filename='sst_t2m_analysis.png'
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
        elif 'valid_time' in ds.coords: actual_time_coord_name = 'valid_time'
        
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
            elif 'valid_time' in ds.coords: time_coord_fallback = 'valid_time'
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

    lat_dim = next((d for d in ['latitude', 'lat'] if d in data_var.dims), None)
    lon_dim = next((d for d in ['longitude', 'lon'] if d in data_var.dims), None)

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

def plot_sst_t2m_combined_analysis(sst_file_path_str, sst_nc_var_name, 
                                   t2m_file_path_str, t2m_nc_var_name, 
                                   output_path_str):
    """
    Plots combined time series and scatter plot for SST and 2m Temperature, showing correlation.
    """
    try:
        output_path = Path(output_path_str)

        series_sst = load_and_average_variable(sst_file_path_str, sst_nc_var_name, 
                                               default_time_coord='time', var_role="SST")
        if series_sst is None: sys.exit(1)
        # SST typically in Celsius, no conversion assumed here.

        series_t2m = load_and_average_variable(t2m_file_path_str, t2m_nc_var_name, 
                                               default_time_coord='valid_time', var_role="T2M")
        if series_t2m is None: sys.exit(1)

        # Unit conversion for T2M (Kelvin to Celsius if needed)
        t2m_display_units = "°C" # Target for T2M
        if series_t2m.mean() > 100 : # Heuristic for Kelvin
            series_t2m = series_t2m - 273.15
            print(f"T2M data appears to be in Kelvin, converted to Celsius.", file=sys.stdout)
        
        # Ensure indices are DatetimeIndex before merging/resampling
        if not isinstance(series_sst.index, pd.DatetimeIndex):
            try: series_sst.index = pd.to_datetime(series_sst.index)
            except: print(f"Error: SST series index could not be converted to DatetimeIndex.", file=sys.stderr); sys.exit(1)
        if not isinstance(series_t2m.index, pd.DatetimeIndex):
            try: series_t2m.index = pd.to_datetime(series_t2m.index)
            except: print(f"Error: T2M series index could not be converted to DatetimeIndex.", file=sys.stderr); sys.exit(1)

        try:
            series_sst_aligned = series_sst.resample('ME').mean()
            series_t2m_aligned = series_t2m.resample('ME').mean()
            print("SST and T2M data resampled to monthly frequency for alignment.", file=sys.stdout)
        except Exception as e_resample:
            print(f"Warning: Could not resample SST/T2M to monthly: {e_resample}. Using original time resolution.", file=sys.stderr)
            series_sst_aligned = series_sst
            series_t2m_aligned = series_t2m

        common_df = pd.merge(
            series_sst_aligned.rename('SST'),
            series_t2m_aligned.rename('T2M'),
            left_index=True, right_index=True,
            how='inner'
        ).dropna()

        if common_df.empty or len(common_df) < 2:
            print("Error: No common overlapping data points found between SST and T2M. Cannot compute correlation.", file=sys.stderr)
            fig, ax = plt.subplots(figsize=(12,6))
            ax.text(0.5,0.5, "Insufficient overlapping data for SST-T2M analysis", horizontalalignment='center', verticalalignment='center')
            ax.set_title("SST and 2m Temperature Combined Analysis")
        else:
            correlation, p_value = stats.pearsonr(common_df['SST'], common_df['T2M'])
            
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=False, gridspec_kw={'height_ratios': [3, 1]})
            
            ax1.plot(common_df.index, common_df['SST'], 'b-', label='SST (°C)')
            ax1.set_ylabel('SST (°C)', color='blue', fontsize=12)
            ax1.tick_params(axis='y', labelcolor='blue')

            ax1_twin = ax1.twinx()
            ax1_twin.plot(common_df.index, common_df['T2M'], 'r-', label=f'T2M ({t2m_display_units})')
            ax1_twin.set_ylabel(f'T2M ({t2m_display_units})', color='red', fontsize=12)
            ax1_twin.tick_params(axis='y', labelcolor='red')
            
            ax1.set_xlabel("Time", fontsize=12)
            ax1.set_title(f"SST and 2m Temperature Time Series\nCorrelation: {correlation:.3f} (p-value: {p_value:.3e})", fontsize=14)
            
            lines, labels = ax1.get_legend_handles_labels()
            lines2, labels2 = ax1_twin.get_legend_handles_labels()
            ax1_twin.legend(lines + lines2, labels + labels2, loc='upper left', fontsize=10)
            ax1.grid(True, linestyle=':', alpha=0.7)

            ax2.scatter(common_df['SST'], common_df['T2M'], alpha=0.6, color='green')
            ax2.set_xlabel('SST (°C)', fontsize=12)
            ax2.set_ylabel(f'T2M ({t2m_display_units})', fontsize=12)
            ax2.set_title('SST vs. 2m Temperature Scatter Plot with Trend Line', fontsize=12)
            
            z = np.polyfit(common_df['SST'], common_df['T2M'], 1)
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
        print(f"An unexpected error occurred in plot_sst_t2m_combined_analysis: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a combined analysis plot for SST and 2m Temperature.")
    parser.add_argument('--files', required=True, nargs='+',
                        help='Paths to the input NetCDF files (SST and T2M sources).')
    parser.add_argument('--varnames', required=True, nargs='+',
                        help="NetCDF variable names. Order should match --vars.")
    parser.add_argument('--vars', required=True, nargs=2, # Expects exactly 2 internal keys
                        help="Two internal variable keys (e.g., 'noaa_sst' 'ecmwf_t2m').")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type (unused by this script).')

    args = parser.parse_args()

    if len(args.varnames) != len(args.vars) or len(args.vars) != 2:
        print("Error: This script expects exactly two --vars (internal keys) and a corresponding number of --varnames.", file=sys.stderr)
        sys.exit(1)
    
    SST_KEY_IDENTIFIER_SUBSTRING = 'sst' 
    T2M_KEY_IDENTIFIER_SUBSTRING = 't2m' # Could also be '2t' or similar if keys vary

    SST_FILE_BASENAME_EXPECTED = "noaa_sst_masked.nc"
    T2M_FILE_BASENAME_EXPECTED = "ecmwf datas.nc" # T2M is from ECMWF

    sst_file_arg, t2m_file_arg = None, None
    sst_varname_arg, t2m_varname_arg = None, None
    
    key1, key2 = args.vars[0], args.vars[1]
    varname1, varname2 = args.varnames[0], args.varnames[1]

    if SST_KEY_IDENTIFIER_SUBSTRING in key1.lower() and T2M_KEY_IDENTIFIER_SUBSTRING in key2.lower():
        sst_varname_arg, t2m_varname_arg = varname1, varname2
    elif T2M_KEY_IDENTIFIER_SUBSTRING in key1.lower() and SST_KEY_IDENTIFIER_SUBSTRING in key2.lower():
        t2m_varname_arg, sst_varname_arg = varname1, varname2
    else:
        print(f"Warning: Could not definitively identify SST and T2M roles from --vars: {args.vars}. Ensure one key contains '{SST_KEY_IDENTIFIER_SUBSTRING}' and the other '{T2M_KEY_IDENTIFIER_SUBSTRING}'.", file=sys.stderr)
        sys.exit(1)

    for f_path_str in args.files:
        if Path(f_path_str).name == SST_FILE_BASENAME_EXPECTED:
            sst_file_arg = f_path_str
        if Path(f_path_str).name == T2M_FILE_BASENAME_EXPECTED:
            t2m_file_arg = f_path_str
            
    if not sst_file_arg and len(args.files) == 1:
        if Path(args.files[0]).name == SST_FILE_BASENAME_EXPECTED or SST_FILE_BASENAME_EXPECTED == T2M_FILE_BASENAME_EXPECTED:
            sst_file_arg = args.files[0]
    
    if not t2m_file_arg and len(args.files) == 1:
        if Path(args.files[0]).name == T2M_FILE_BASENAME_EXPECTED or SST_FILE_BASENAME_EXPECTED == T2M_FILE_BASENAME_EXPECTED:
            t2m_file_arg = args.files[0]

    if not sst_file_arg or not t2m_file_arg or not sst_varname_arg or not t2m_varname_arg:
        print("Error: Could not determine necessary file paths or variable names for SST and T2M.", file=sys.stderr)
        print(f"  SST File: {sst_file_arg}, SST Var: {sst_varname_arg}")
        print(f"  T2M File: {t2m_file_arg}, T2M Var: {t2m_varname_arg}")
        print(f"  Args.files: {args.files}, Args.vars: {args.vars}, Args.varnames: {args.varnames}")
        sys.exit(1)
    
    print(f"Plotting SST (file: '{sst_file_arg}', var: '{sst_varname_arg}') vs T2M (file: '{t2m_file_arg}', var: '{t2m_varname_arg}')", file=sys.stdout)

    plot_sst_t2m_combined_analysis(sst_file_arg, sst_varname_arg, 
                                   t2m_file_arg, t2m_varname_arg, 
                                   args.output)