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


# # In[10]:


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
#     ecmwf_df, 'msl',
#     title='Combined Analysis: SST and Mean Sea Level Pressure',
#     filename='sst_msl_analysis.png'
# )


# # In[ ]:




#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# from matplotlib.dates import DateFormatter # Not used
import scipy.stats as stats
# import seaborn as sns # Not directly used for this plot
# import cartopy.crs as ccrs # Not used
# import matplotlib.dates as mdates # Not used
# from matplotlib.colors import LinearSegmentedColormap # Not used
# from matplotlib.ticker import MaxNLocator # Not used

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
        elif 'time' in ds.coords: actual_time_coord_name = 'time' # Common fallback
        elif 'valid_time' in ds.coords: actual_time_coord_name = 'valid_time' # Another common one
        
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
    
    # Ensure time index is DatetimeIndex for to_series()
    if not isinstance(ds[actual_time_coord_name].to_index(), pd.DatetimeIndex):
        try:
            ds[actual_time_coord_name] = pd.to_datetime(ds[actual_time_coord_name].values)
        except Exception as e_time_convert:
            print(f"Warning: Could not convert time index to DatetimeIndex for {var_role} from {file_path}: {e_time_convert}", file=sys.stderr)
            # Continue, but to_series() might have a non-DatetimeIndex

    return spatial_avg.to_series().rename(nc_var_name) # Rename series to nc_var_name

def plot_sst_msl_combined_analysis(sst_file_path_str, sst_nc_var_name, 
                                   msl_file_path_str, msl_nc_var_name, 
                                   output_path_str):
    """
    Plots combined time series and scatter plot for SST and MSL, showing correlation.
    """
    try:
        output_path = Path(output_path_str)

        # Load and process SST
        series_sst = load_and_average_variable(sst_file_path_str, sst_nc_var_name, 
                                               default_time_coord='time', var_role="SST")
        if series_sst is None: sys.exit(1)
        # SST is typically in Celsius, no major unit conversion usually needed here from common sources.

        # Load and process MSL
        series_msl = load_and_average_variable(msl_file_path_str, msl_nc_var_name, 
                                               default_time_coord='valid_time', var_role="MSL")
        if series_msl is None: sys.exit(1)

        # Unit conversion for MSL (Pa to hPa if needed)
        msl_display_units = "hPa" # Target
        if series_msl.max() > 10000 : # Heuristic for Pa
            series_msl = series_msl / 100
            print(f"MSL data appears to be in Pa, converted to hPa.", file=sys.stdout)
        
        # Merge dataframes on their time index (must be DatetimeIndex)
        # Ensure indices are DatetimeIndex before merging
        if not isinstance(series_sst.index, pd.DatetimeIndex):
            try: series_sst.index = pd.to_datetime(series_sst.index)
            except: print(f"Error: SST series index could not be converted to DatetimeIndex.", file=sys.stderr); sys.exit(1)
        if not isinstance(series_msl.index, pd.DatetimeIndex):
            try: series_msl.index = pd.to_datetime(series_msl.index)
            except: print(f"Error: MSL series index could not be converted to DatetimeIndex.", file=sys.stderr); sys.exit(1)

        # Resample to a common frequency (e.g., monthly) to ensure alignment
        # This helps if one dataset is daily and the other is monthly, for example.
        try:
            series_sst_aligned = series_sst.resample('ME').mean()
            series_msl_aligned = series_msl.resample('ME').mean()
            print("SST and MSL data resampled to monthly frequency for alignment.", file=sys.stdout)
        except Exception as e_resample:
            print(f"Warning: Could not resample SST/MSL to monthly: {e_resample}. Using original time resolution.", file=sys.stderr)
            series_sst_aligned = series_sst
            series_msl_aligned = series_msl

        common_df = pd.merge(
            series_sst_aligned.rename('SST'), # Rename for clarity in df
            series_msl_aligned.rename('MSL'), # Rename for clarity in df
            left_index=True, right_index=True,
            how='inner' # Use 'inner' to keep only overlapping time points
        ).dropna() # Drop rows with NaNs that might result from resampling or initial data

        if common_df.empty or len(common_df) < 2: # Need at least 2 points for correlation
            print("Error: No common overlapping data points found between SST and MSL after alignment and NaN removal. Cannot compute correlation.", file=sys.stderr)
            # Create an empty plot with a message
            fig, ax = plt.subplots(figsize=(12,6))
            ax.text(0.5,0.5, "Insufficient overlapping data for SST-MSL analysis", horizontalalignment='center', verticalalignment='center')
            ax.set_title("SST and MSL Combined Analysis")
        else:
            # Calculate correlation
            correlation, p_value = stats.pearsonr(common_df['SST'], common_df['MSL'])
            
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=False, gridspec_kw={'height_ratios': [3, 1]})
            
            # Time series plot
            ax1.plot(common_df.index, common_df['SST'], 'b-', label='SST (°C)') # Assuming SST is in Celsius
            ax1.set_ylabel('SST (°C)', color='blue', fontsize=12)
            ax1.tick_params(axis='y', labelcolor='blue')

            ax1_twin = ax1.twinx()
            ax1_twin.plot(common_df.index, common_df['MSL'], 'r-', label=f'MSL ({msl_display_units})')
            ax1_twin.set_ylabel(f'MSL ({msl_display_units})', color='red', fontsize=12)
            ax1_twin.tick_params(axis='y', labelcolor='red')
            
            ax1.set_xlabel("Time", fontsize=12)
            ax1.set_title(f"SST and MSL Time Series\nCorrelation: {correlation:.3f} (p-value: {p_value:.3e})", fontsize=14)
            
            # Add legend for twin axes
            lines, labels = ax1.get_legend_handles_labels()
            lines2, labels2 = ax1_twin.get_legend_handles_labels()
            ax1_twin.legend(lines + lines2, labels + labels2, loc='upper left', fontsize=10)
            ax1.grid(True, linestyle=':', alpha=0.7)

            # Scatter plot
            ax2.scatter(common_df['SST'], common_df['MSL'], alpha=0.6, color='purple')
            ax2.set_xlabel('SST (°C)', fontsize=12)
            ax2.set_ylabel(f'MSL ({msl_display_units})', fontsize=12)
            ax2.set_title('SST vs. MSL Scatter Plot with Trend Line', fontsize=12)
            
            # Add trend line
            z = np.polyfit(common_df['SST'], common_df['MSL'], 1)
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
        print(f"An unexpected error occurred in plot_sst_msl_combined_analysis: {e}", file=sys.stderr)
        # import traceback # For debugging
        # traceback.print_exc(file=sys.stderr) # For debugging
        sys.exit(1)


# Replace the if __name__ == "__main__": block in your
# plotting_scripts/sst - msl combined analysis.py script with this:

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a combined analysis plot for SST and MSL.")
    parser.add_argument('--files', required=True, nargs='+',
                        help='Paths to the unique input NetCDF files involved (SST and MSL sources).')
    parser.add_argument('--varnames', required=True, nargs='+',
                        help="NetCDF variable names. Order should match --vars.")
    parser.add_argument('--vars', required=True, nargs=2, # Expects exactly 2 internal keys
                        help="Two internal variable keys (e.g., 'noaa_sst' 'ecmwf_msl').")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type (unused by this script).')

    args = parser.parse_args()

    if len(args.varnames) != len(args.vars) or len(args.vars) != 2:
        print("Error: This script expects exactly two --vars (internal keys) and a corresponding number of --varnames.", file=sys.stderr)
        print(f"Received --vars: {args.vars}, --varnames: {args.varnames}", file=sys.stderr)
        sys.exit(1)
    
    # Expected internal key characteristics (from Streamlit config's VARIABLE_DETAILS_MAP['file_key'])
    # These are used to identify which var is SST and which is MSL from args.vars
    SST_KEY_IDENTIFIER_SUBSTRING = 'sst' 
    MSL_KEY_IDENTIFIER_SUBSTRING = 'msl'

    # Expected file basenames (from Streamlit config's VARIABLE_DETAILS_MAP['data_file'].name)
    # These help map the unique file paths in args.files back to SST or MSL roles.
    SST_FILE_BASENAME_EXPECTED = "noaa_sst_masked.nc"  # As defined in your Streamlit app for NOAA_SST_FILE
    MSL_FILE_BASENAME_EXPECTED = "ecmwf datas.nc"    # As defined in your Streamlit app for ECMWF_FILE

    sst_file_arg, msl_file_arg = None, None
    sst_varname_arg, msl_varname_arg = None, None
    
    # Assign roles (SST vs MSL) based on args.vars content
    # args.vars contains the two internal keys passed by Streamlit
    # args.varnames contains the corresponding nc_var_names in the same order
    
    key1_is_sst = SST_KEY_IDENTIFIER_SUBSTRING in args.vars[0].lower()
    key1_is_msl = MSL_KEY_IDENTIFIER_SUBSTRING in args.vars[0].lower()
    key2_is_sst = SST_KEY_IDENTIFIER_SUBSTRING in args.vars[1].lower()
    key2_is_msl = MSL_KEY_IDENTIFIER_SUBSTRING in args.vars[1].lower()

    if key1_is_sst and key2_is_msl:
        sst_varname_arg = args.varnames[0]
        msl_varname_arg = args.varnames[1]
    elif key1_is_msl and key2_is_sst:
        msl_varname_arg = args.varnames[0]
        sst_varname_arg = args.varnames[1]
    else:
        # Fallback or error if roles are ambiguous from keys
        print(f"Warning: Could not definitively identify SST and MSL roles from --vars: {args.vars}. Check key substrings.", file=sys.stderr)
        # As a simple fallback, assume order if specific identifiers are missing, though this is less robust.
        # For this specific SST-MSL script, it's better to be strict.
        print("Ensure one internal key in --vars contains 'sst' and the other 'msl'.", file=sys.stderr)
        sys.exit(1)

    # Now, find the actual file paths for SST and MSL from the unique list in args.files
    for f_path_str in args.files:
        if Path(f_path_str).name == SST_FILE_BASENAME_EXPECTED:
            sst_file_arg = f_path_str
        if Path(f_path_str).name == MSL_FILE_BASENAME_EXPECTED:
            msl_file_arg = f_path_str
            
    # If SST and MSL might come from the same physical file (e.g., if ECMWF_FILE also contained SST)
    # and only one unique file path was passed in args.files.
    if not sst_file_arg and len(args.files) == 1:
        # Check if the single passed file is the expected SST file OR if it's the MSL file (and MSL file is same as SST file)
        if Path(args.files[0]).name == SST_FILE_BASENAME_EXPECTED or SST_FILE_BASENAME_EXPECTED == MSL_FILE_BASENAME_EXPECTED:
            sst_file_arg = args.files[0]
    
    if not msl_file_arg and len(args.files) == 1:
        if Path(args.files[0]).name == MSL_FILE_BASENAME_EXPECTED or SST_FILE_BASENAME_EXPECTED == MSL_FILE_BASENAME_EXPECTED:
            msl_file_arg = args.files[0]

    if not sst_file_arg or not msl_file_arg:
        print("Error: Could not determine the file path for SST and/or MSL.", file=sys.stderr)
        print(f"  Expected SST file basename: {SST_FILE_BASENAME_EXPECTED}, Found: {sst_file_arg}")
        print(f"  Expected MSL file basename: {MSL_FILE_BASENAME_EXPECTED}, Found: {msl_file_arg}")
        print(f"  Unique files passed (--files): {args.files}")
        sys.exit(1)
        
    if not sst_varname_arg or not msl_varname_arg:
        print("Error: Could not determine NetCDF variable names for SST and/or MSL.", file=sys.stderr)
        sys.exit(1)

    print(f"Plotting SST (file: '{sst_file_arg}', var: '{sst_varname_arg}') vs MSL (file: '{msl_file_arg}', var: '{msl_varname_arg}')", file=sys.stdout)

    plot_sst_msl_combined_analysis(sst_file_arg, sst_varname_arg, 
                                   msl_file_arg, msl_varname_arg, 
                                   args.output)