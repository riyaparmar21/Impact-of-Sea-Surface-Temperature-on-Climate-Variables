# #!/usr/bin/env python
# # coding: utf-8

# # In[1]:


# import xarray as xr
# import pandas as pd
# import numpy as np
# import seaborn as sns
# import matplotlib.pyplot as plt


# # In[2]:


# noaa_precip_ds = xr.open_dataset(r"precip_konkan.nc")
# ecmwf_ds = xr.open_dataset(r"ecmwf datas.nc")
# noaa_sst_ds = xr.open_dataset(r"noaa_sst_masked.nc")

# # Compute spatial averages:
# # NOAA Precipitation (dimensions: time, lat, lon)
# noaa_precip = noaa_precip_ds['precip'].mean(dim=['lat', 'lon'])
# # NOAA SST (dimensions: time, lat, lon)
# noaa_sst = noaa_sst_ds['sst'].mean(dim=['lat', 'lon'])
# # ECMWF variables (dimensions: valid_time, latitude, longitude)
# # Rename the coordinate 'valid_time' to 'time' for merging without altering variable names
# ecmwf_msl = ecmwf_ds['msl'].mean(dim=['latitude', 'longitude']).rename({'valid_time': 'time'})
# ecmwf_tcc = ecmwf_ds['tcc'].mean(dim=['latitude', 'longitude']).rename({'valid_time': 'time'})
# ecmwf_u10 = ecmwf_ds['u10'].mean(dim=['latitude', 'longitude']).rename({'valid_time': 'time'})
# ecmwf_v10 = ecmwf_ds['v10'].mean(dim=['latitude', 'longitude']).rename({'valid_time': 'time'})
# ecmwf_t2m = ecmwf_ds['t2m'].mean(dim=['latitude', 'longitude']).rename({'valid_time': 'time'})


# # Convert each DataArray to a pandas Series
# series_noaa_precip = noaa_precip.to_series()
# series_noaa_sst = noaa_sst.to_series()
# series_ecmwf_msl = ecmwf_msl.to_series()
# series_ecmwf_tcc = ecmwf_tcc.to_series()
# series_ecmwf_u10 = ecmwf_u10.to_series()
# series_ecmwf_v10 = ecmwf_v10.to_series()
# series_ecmwf_t2m = ecmwf_t2m.to_series()


# # Merge all series into a single DataFrame using an outer join on the time index
# df = pd.concat([
#     series_noaa_precip,
#     series_noaa_sst,
#     series_ecmwf_msl,
#     series_ecmwf_tcc,
#     series_ecmwf_u10,
#     series_ecmwf_v10,
#     series_ecmwf_t2m,
    
# ], axis=1)

# # Set column names; ECMWF variables keep their original names.
# df.columns = ['noaa_precip', 'noaa_sst', 'msl', 'tcc', 'u10', 'v10', 't2m']

# # Compute the correlation matrix (NaNs will be handled pairwise)
# corr_matrix = df.corr()

# print("Correlation Matrix:")
# print(corr_matrix)


# # In[3]:


# print("The mosdac rows and columns are empty because it has 100% NaN values")
# # Plot the correlation matrix using seaborn heatmap
# plt.figure(figsize=(10, 8))
# sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', fmt=".2f")
# plt.title("Correlation Matrix Heatmap")
# plt.show()


# # In[ ]:


#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import sys
from pathlib import Path

def plot_correlation_matrix(input_files_dict, output_image_path_str, output_csv_path_str):
    """
    Computes and plots a correlation matrix for the provided variables.
    Saves the plot as an image and the correlation data as a CSV.

    Parameters:
    -----------
    input_files_dict : dict
        Dictionary where keys are internal variable keys (e.g., 'noaa_sst', 'ecmwf_msl')
        and values are dicts {'file_path': str, 'var_name_in_nc': str, 'user_friendly_name': str}.
        'user_friendly_name' will be used as column names in the DataFrame.
    output_image_path_str : str
        Path to save the output plot image.
    output_csv_path_str : str
        Path to save the output correlation data CSV.
    """
    all_series = []
    column_names = []

    try:
        for internal_key, var_info in input_files_dict.items():
            file_path = Path(var_info['file_path'])
            nc_var_name = var_info['var_name_in_nc']
            user_friendly_col_name = var_info['user_friendly_name'] # e.g., 'sst', 'precip'

            if not file_path.exists():
                print(f"Error: Data file not found at {file_path} for variable {internal_key}", file=sys.stderr)
                continue # Skip this variable but try others

            try:
                ds = xr.open_dataset(file_path, decode_times=True)
                actual_time_coord_name = None
                if 'time' in ds.coords: actual_time_coord_name = 'time'
                elif 'valid_time' in ds.coords: actual_time_coord_name = 'valid_time'
                
                if not actual_time_coord_name and 'time' not in ds.dims and 'valid_time' not in ds.dims: # Check dimensions too
                     print(f"Warning: No 'time' or 'valid_time' coordinate found for {internal_key} in {file_path}. Trying to find a time-like dimension.", file=sys.stderr)
                     time_like_dims = [d for d in ds[nc_var_name].dims if 'time' in d.lower()]
                     if time_like_dims:
                         actual_time_coord_name = time_like_dims[0]
                         print(f"Using dimension '{actual_time_coord_name}' as time for {internal_key}.", file=sys.stdout)
                     else:
                        print(f"Error: Could not identify time coordinate/dimension for {internal_key} in {file_path}. Skipping.", file=sys.stderr)
                        continue


                if actual_time_coord_name and actual_time_coord_name != 'time' and actual_time_coord_name in ds.coords:
                    ds = ds.rename({actual_time_coord_name: 'time'})
                # If actual_time_coord_name was a dimension but not a coordinate, renaming might not be needed or might fail.
                # The key is that the dimension used for time operations is consistent or correctly identified.
                
                # Ensure 'time' is now the time coordinate we work with for this dataset
                if 'time' not in ds.coords:
                     print(f"Error: 'time' coordinate not established for {internal_key} in {file_path} after attempting rename. Skipping.", file=sys.stderr)
                     continue


            except Exception as e:
                print(f"Warning: Could not load or decode time for {file_path} (var: {internal_key}): {e}. Trying decode_times=False.", file=sys.stderr)
                try:
                    ds = xr.open_dataset(file_path, decode_times=False)
                    # Manual decode logic (simplified, adapt if complex time units)
                    time_coord_name_fallback = None
                    if 'time' in ds.coords: time_coord_name_fallback = 'time'
                    elif 'valid_time' in ds.coords: time_coord_name_fallback = 'valid_time'
                    
                    if time_coord_name_fallback:
                        time_ds_to_decode = xr.Dataset({time_coord_name_fallback: ds[time_coord_name_fallback]})
                        decoded_time_ds = xr.decode_cf(time_ds_to_decode)
                        ds[time_coord_name_fallback] = decoded_time_ds[time_coord_name_fallback]
                        if time_coord_name_fallback == 'valid_time':
                            ds = ds.rename({'valid_time': 'time'})
                    else:
                        print(f"Error: Failed to find time coordinate even in fallback for {file_path}. Skipping {internal_key}.", file=sys.stderr)
                        continue
                    if 'time' not in ds.coords:
                         print(f"Error: 'time' coordinate still not established after fallback for {internal_key}. Skipping.", file=sys.stderr)
                         continue
                except Exception as e_fb:
                    print(f"Error: Failed to load dataset {file_path} (var: {internal_key}) even with fallback: {e_fb}. Skipping.", file=sys.stderr)
                    continue

            if nc_var_name not in ds:
                print(f"Error: Variable '{nc_var_name}' not found in {file_path} for {internal_key}. Skipping.", file=sys.stderr)
                continue

            data_var = ds[nc_var_name]

            # Determine spatial dimensions for averaging
            lat_dim = 'latitude' if 'latitude' in data_var.dims else 'lat'
            lon_dim = 'longitude' if 'longitude' in data_var.dims else 'lon'

            if lat_dim not in data_var.dims or lon_dim not in data_var.dims:
                print(f"Warning: Could not find standard lat/lon dimensions for {internal_key} in {file_path}. Dims: {data_var.dims}. Attempting mean over all non-time dims.", file=sys.stderr)
                # Fallback: mean over all dimensions that are not 'time'
                non_time_dims = [d for d in data_var.dims if d != 'time'] # Assumes time dim is 'time'
                if not non_time_dims:
                    print(f"Error: No non-time dimensions found to average over for {internal_key}. Skipping.", file=sys.stderr)
                    continue
                spatial_avg_var = data_var.mean(dim=non_time_dims, skipna=True)
            else:
                spatial_avg_var = data_var.mean(dim=[lat_dim, lon_dim], skipna=True)

            # Convert to pandas Series, ensure time index is DatetimeIndex
            series = spatial_avg_var.to_series()
            if not isinstance(series.index, pd.DatetimeIndex):
                try:
                    series.index = pd.to_datetime(series.index)
                except Exception as e_idx:
                    print(f"Warning: Could not convert index to DatetimeIndex for {internal_key}. Data might not align: {e_idx}", file=sys.stderr)
            
            all_series.append(series)
            column_names.append(user_friendly_col_name) # Use user-friendly name for column

        if not all_series:
            print("Error: No data series were successfully processed. Cannot create correlation matrix.", file=sys.stderr)
            sys.exit(1)

        # Merge all series into a single DataFrame using an outer join on the time index
        # This aligns data by time, filling NaNs where time points don't match.
        df_merged = pd.concat(all_series, axis=1)
        df_merged.columns = column_names

        # Resample to a common frequency (e.g., monthly 'M') if necessary, to handle slightly misaligned daily/hourly data
        # This step is crucial if time indices are not perfectly aligned for all variables.
        # Using monthly mean as a common strategy.
        try:
            # Only resample if index is DatetimeIndex
            if isinstance(df_merged.index, pd.DatetimeIndex):
                df_aligned = df_merged.resample('ME').mean() # Resample to monthly mean
                print("Data resampled to monthly frequency for alignment.", file=sys.stdout)
            else:
                print("Warning: Index is not DatetimeIndex, cannot resample. Correlations might be based on misaligned data if time steps differ.", file=sys.stderr)
                df_aligned = df_merged # Use as is if not datetime index
        except Exception as e_resample:
            print(f"Warning: Could not resample DataFrame to monthly: {e_resample}. Using original alignment.", file=sys.stderr)
            df_aligned = df_merged


        # Compute the correlation matrix (NaNs will be handled pairwise by default by .corr())
        corr_matrix = df_aligned.corr()

        # Save the correlation matrix to CSV
        try:
            output_csv_path = Path(output_csv_path_str)
            corr_matrix.to_csv(output_csv_path)
            print(f"Correlation matrix data saved to {output_csv_path}", file=sys.stdout)
        except Exception as e_csv:
            print(f"Error: Failed to save correlation data to CSV {output_csv_path_str}: {e_csv}", file=sys.stderr)
            # Continue to plot image even if CSV fails

        # Plot the correlation matrix
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', fmt=".2f", vmin=-1, vmax=1)
        plt.title("Correlation Matrix of Spatially Averaged Variables", fontsize=14)
        plt.xticks(rotation=45, ha='right', fontsize=10)
        plt.yticks(rotation=0, fontsize=10)
        plt.tight_layout()
        
        try:
            output_image_path = Path(output_image_path_str)
            plt.savefig(output_image_path, dpi=300, bbox_inches='tight')
            print(f"Correlation matrix plot saved to {output_image_path}", file=sys.stdout)
        except Exception as e_img:
            print(f"Error: Failed to save correlation plot to {output_image_path_str}: {e_img}", file=sys.stderr)
            sys.exit(1) # Exit if image saving fails, as it's the primary output
        finally:
            plt.close()

    except Exception as e:
        print(f"An unexpected error occurred in plot_correlation_matrix: {e}", file=sys.stderr)
        # import traceback # For debugging
        # traceback.print_exc(file=sys.stderr) # For debugging
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute and plot a correlation matrix for given variables.")
    parser.add_argument('--files', required=True, nargs='+', help='Paths to the unique input NetCDF files involved.')
    parser.add_argument('--varnames', required=True, nargs='+', help="NetCDF variable names corresponding to each internal key in --vars.")
    parser.add_argument('--vars', required=True, nargs='+', help="Internal variable keys (e.g., 'noaa_sst', 'ecmwf_msl') for all variables to correlate.")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image for the heatmap.')
    parser.add_argument('--plot_type', help='Name of the plot type (e.g., "Correlation Matrix Plot").')

    args = parser.parse_args()

    # --- THIS CHECK WAS INCORRECT AND IS REMOVED/COMMENTED ---
    # if not (len(args.files) == len(args.varnames) == len(args.vars)):
    #     print("Error: The number of --files, --varnames, and --vars must be the same.", file=sys.stderr)
    #     print(f"Got {len(args.files)} files, {len(args.varnames)} varnames, {len(args.vars)} vars.", file=sys.stderr)
    #     sys.exit(1)
    # --- END OF REMOVED CHECK ---

    # Ensure that the number of varnames matches the number of vars (internal_keys)
    if len(args.varnames) != len(args.vars):
        print("Error: The number of --varnames must match the number of --vars (internal keys).", file=sys.stderr)
        print(f"Got {len(args.varnames)} varnames and {len(args.vars)} vars.", file=sys.stderr)
        sys.exit(1)
    
    # The number of files can be less than vars/varnames if multiple vars come from the same file.

    try:
        # Reconstruct the input_files_dict needed by the plotting function
        reconstructed_input_files_dict = {}
        
        # This temporary map helps find the original filename for an internal_key.
        # In a real-world scenario, this mapping should ideally come from a shared config
        # or be passed more directly from the Streamlit app if it's dynamic.
        # For now, this mirrors the structure from your Streamlit app's VARIABLE_DETAILS_MAP.
        temp_var_to_file_basename_map = {
            'noaa_precip': "precip_konkan.nc",
            'noaa_sst':    "noaa_sst_masked.nc",
            'ecmwf_msl':   "ecmwf datas.nc",
            'ecmwf_tcc':   "ecmwf datas.nc",
            'ecmwf_u10':   "ecmwf datas.nc",
            'ecmwf_v10':   "ecmwf datas.nc",
            'ecmwf_t2m':   "ecmwf datas.nc",
        }
        # User-friendly names (can also be derived from internal_key if simple pattern)
        temp_internal_to_user_friendly_map = {
            'noaa_precip': 'precip',
            'noaa_sst':    'sst',
            'ecmwf_msl':   'msl',
            'ecmwf_tcc':   'tcc',
            'ecmwf_u10':   'u10',
            'ecmwf_v10':   'v10',
            'ecmwf_t2m':   't2m',
        }


        for i, internal_key in enumerate(args.vars):
            nc_var_name = args.varnames[i]
            user_friendly_name = temp_internal_to_user_friendly_map.get(internal_key, internal_key) # Fallback to internal_key
            
            original_file_basename = temp_var_to_file_basename_map.get(internal_key)
            actual_file_path_for_var = None
            
            if original_file_basename:
                for fp_str_from_args in args.files: # args.files contains unique, full paths
                    if Path(fp_str_from_args).name == original_file_basename:
                        actual_file_path_for_var = fp_str_from_args
                        break
            
            if not actual_file_path_for_var:
                print(f"Warning: Could not determine file path for internal key '{internal_key}' (NetCDF var: '{nc_var_name}') using basename '{original_file_basename}'. Searched in {args.files}. Skipping this variable.", file=sys.stderr)
                continue

            reconstructed_input_files_dict[internal_key] = {
                'file_path': actual_file_path_for_var,
                'var_name_in_nc': nc_var_name,
                'user_friendly_name': user_friendly_name
            }

        if not reconstructed_input_files_dict:
            print("Error: Could not reconstruct any variable details for processing. Exiting.", file=sys.stderr)
            sys.exit(1)

        output_csv_path = args.output.replace('.png', '_data.csv') # Convention
        plot_correlation_matrix(reconstructed_input_files_dict, args.output, output_csv_path)

    except Exception as e:
        print(f"Critical error in main execution of correlation matrix script: {e}", file=sys.stderr)
        # import traceback # For debugging
        # traceback.print_exc(file=sys.stderr) # For debugging
        sys.exit(1)