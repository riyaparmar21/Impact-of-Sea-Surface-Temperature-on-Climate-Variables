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


# # In[98]:


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
    
#     # For Total Cloud Cover (tcc)
#     # Seasonal cloud cover
#     tcc_mean = ds.tcc.mean(dim=['latitude', 'longitude'])
#     seasons = {'December-January-February': [12, 1, 2], 'March-April-May': [3, 4, 5],
#                'June-July-August': [6, 7, 8], 'September-October-November': [9, 10, 11]}
   
#     seasonal_tcc = {season: [] for season in seasons}
#     for i, t in enumerate(times):
#         for season, months in seasons.items():
#             if t.month in months:
#                 seasonal_tcc[season].append(tcc_mean.values[i])
   
#     seasonal_tcc_mean = {season: np.mean(vals) for season, vals in seasonal_tcc.items()}
   
#     plt.figure(figsize=(8, 5))
#     plt.bar(seasonal_tcc_mean.keys(), seasonal_tcc_mean.values(), color='skyblue')
#     plt.xlabel('Season')
#     plt.ylabel('Mean Total Cloud Cover')
#     plt.title('Seasonal Mean Cloud Cover')
#     plt.grid(axis='y', linestyle='--', alpha=0.7)
#     plt.tight_layout()
#     plt.show()
#     #plt.savefig('ecmwf_tcc_seasonal.png')
#     plt.close()


# # In[99]:


# plot_ecmwf_vars(r"ecmwf datas.nc")


# # In[ ]:




#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd # For pd.to_datetime if time conversion is tricky
# from matplotlib.dates import DateFormatter # Not used
# import scipy.stats as stats # Not used
# import seaborn as sns # Not directly used, can be for colors
# import cartopy.crs as ccrs # Not used
# import matplotlib.dates as mdates # Not used
# from matplotlib.colors import LinearSegmentedColormap # Not used
# from matplotlib.ticker import MaxNLocator # Not used

import argparse
import sys
from pathlib import Path

def plot_tcc_seasonal_mean_bar(ecmwf_file_path_str, nc_variable_name, output_path_str):
    """
    Plots the seasonal mean Total Cloud Cover (TCC) as a bar chart.
   
    Parameters:
    -----------
    ecmwf_file_path_str : str
        Path to the ECMWF data file.
    nc_variable_name : str
        The NetCDF variable name for TCC (e.g., 'tcc').
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
            print(f"Error: Variable '{nc_variable_name}' (expected for TCC) not found in dataset {ecmwf_file_path}", file=sys.stderr)
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
            # print(f"Using dimension '{dim_to_use_for_time}' for time operations on variable '{nc_variable_name}'.", file=sys.stdout)

        # Determine spatial coordinate names for averaging
        lat_coord_name = 'latitude' if 'latitude' in data_var.dims else 'lat'
        lon_coord_name = 'longitude' if 'longitude' in data_var.dims else 'lon'

        if lat_coord_name not in data_var.dims or lon_coord_name not in data_var.dims:
            print(f"Error: Could not determine lat/lon dimensions for '{nc_variable_name}'. Dims: {data_var.dims}", file=sys.stderr)
            sys.exit(1)
            
        # Area average first
        tcc_area_mean = data_var.mean(dim=[lat_coord_name, lon_coord_name], skipna=True)

        # Define seasons using standard meteorological definitions
        def get_season(time_coord_da, time_dim_for_da):
            month = time_coord_da.dt.month
            seasons_map = {
                1: 'DJF', 2: 'DJF', 3: 'MAM', 4: 'MAM', 5: 'MAM',
                6: 'JJA', 7: 'JJA', 8: 'JJA', 9: 'SON', 10: 'SON',
                11: 'SON', 12: 'DJF'
            }
            # Use the actual time dimension name for the new coordinate
            return xr.DataArray([seasons_map[m] for m in month.data], coords={time_dim_for_da: time_coord_da[time_dim_for_da]}, name='season')

        # Assign season coordinate to the area-averaged DataArray
        tcc_area_mean_with_season = tcc_area_mean.assign_coords(season=get_season(tcc_area_mean[actual_time_coord_name], dim_to_use_for_time))
        
        # Calculate mean for each season
        seasonal_tcc_mean = tcc_area_mean_with_season.groupby('season').mean(skipna=True) # No dim needed as it's 1D now

        # Convert TCC from fraction (0-1) to percentage (0-100) if applicable
        display_units = data_var.attrs.get('units', 'fraction')
        if seasonal_tcc_mean.max(skipna=True).item() <= 1.0 and seasonal_tcc_mean.min(skipna=True).item() >=0.0 and display_units == 'fraction':
            seasonal_tcc_mean = seasonal_tcc_mean * 100
            display_units = '%' # Update units for label
            print(f"TCC data converted from fraction to percentage.", file=sys.stdout)
        elif display_units != '%': # If not fraction and not already %, keep original
            print(f"TCC data presented in original units: {display_units}.", file=sys.stdout)


        # Ensure seasons are plotted in a standard order
        season_order = ['DJF', 'MAM', 'JJA', 'SON']
        # Filter and reorder the data according to season_order
        plot_data_values = []
        plot_season_labels = []
        for s_name in season_order:
            if s_name in seasonal_tcc_mean.season.values:
                plot_data_values.append(seasonal_tcc_mean.sel(season=s_name).item())
                plot_season_labels.append(s_name)
        
        if not plot_data_values:
            print("Error: No seasonal data available for plotting.", file=sys.stderr)
            # Create an empty plot with a message if you want
            fig, ax = plt.subplots(figsize=(8,5))
            ax.text(0.5, 0.5, "No seasonal data available", horizontalalignment='center', verticalalignment='center')
            ax.set_title('Seasonal Mean Total Cloud Cover')
        else:
            plt.figure(figsize=(8, 5))
            plt.bar(plot_season_labels, plot_data_values, color='skyblue')
            plt.xlabel('Season', fontsize=12)
            plt.ylabel(f'Mean Total Cloud Cover ({display_units})', fontsize=12)
            plt.title('Seasonal Mean Total Cloud Cover', fontsize=14)
            plt.grid(axis='y', linestyle='--', alpha=0.7)
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
    parser = argparse.ArgumentParser(description="Generate a seasonal mean bar plot for Total Cloud Cover (TCC).")
    parser.add_argument('--files', required=True, nargs=1, help='Path to the input ECMWF NetCDF file.')
    parser.add_argument('--varnames', required=True, nargs=1, help="NetCDF variable name for TCC (e.g., 'tcc').")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type (unused by this script).')
    parser.add_argument('--vars', nargs='+', help='Internal variable keys (unused by this script).')
    
    args = parser.parse_args()
        
    ecmwf_file = args.files[0]
    nc_var_name = args.varnames[0]

    plot_tcc_seasonal_mean_bar(ecmwf_file, nc_var_name, args.output)