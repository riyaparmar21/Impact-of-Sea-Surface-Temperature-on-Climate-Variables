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


# # In[51]:


# # Load your datasets
# ds_precip = xr.open_dataset( r"precip_konkan.nc")
# ds_ecmwf = xr.open_dataset( r"ecmwf datas.nc")
# ds_sst = xr.open_dataset( r"noaa_sst_masked.nc")


# # In[60]:


# # Seasonal Averages 
# def plot_seasonal_averages(ds, var_name, time_dim='valid_time', lat_dim='latitude', lon_dim='longitude'):
#     """
#     Plot seasonal averages for a given variable
#     """
#     # Convert dataset time to datetime if not already
#     if not np.issubdtype(ds[time_dim].dtype, np.datetime64):
#         ds[time_dim] = ds[time_dim].astype('datetime64[ns]')
   
#     # Create a copy of the dataset with a 'season' coordinate
#     ds_with_season = ds.copy()
   
#     # Add season information
#     month_to_season = {1: 'Dec-Jan-Feb', 2: 'Dec-Jan-Feb,', 3: 'Mar-Apr-May', 4: 'Mar-Apr-May', 5: 'Mar-Apr-May',
#                        6: 'Jun-Jul-Aug', 7: 'Jun-Jul-Aug', 8: 'Jun-Jul-Aug', 9: 'Sep-Oct-Nov', 10: 'Sep-Oct-Nov',
#                        11: 'Sep-Oct-Nov', 12: 'Dec-Jan-Feb,'}
   
#     # Extract month and map to season
#     months = ds[time_dim].dt.month.values
#     seasons = [month_to_season[m] for m in months]
#     ds_with_season = ds_with_season.assign_coords(season=('valid_time', seasons))
   
#     # Calculate seasonal averages
#     seasonal_avg = ds_with_season[var_name].groupby('season').mean(dim=time_dim)
   
#     # Create figure with 4 subplots (one for each season)
#     fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 14))
#     axes = axes.flatten()
   
#     # Define color maps based on variable
#     if 'precip' in var_name:
#         cmap = 'Blues'
#         label = 'Precipitation (mm)'
#     elif 'sst' in var_name:
#         cmap = 'RdYlBu_r'
#         label = 'Temperature (°C)'
#     elif 't2m' in var_name:
#         cmap = 'RdYlBu_r'
#         label = 'Temperature (K)'
#     else:
#         cmap = 'viridis'
#         label = var_name
   
#     # Plot each season
#     for i, season in enumerate(['Dec-Jan-Feb', 'Mar-Apr-May', 'Jun-Jul-Aug', 'Sep-Oct-Nov']):
#         if season in seasonal_avg.season.values:
#             season_data = seasonal_avg.sel(season=season)
           
#             # Get coordinate names for proper plotting
#             lon_coords = ds[lon_dim]
#             lat_coords = ds[lat_dim]
           
#             # Create heatmap
#             im = axes[i].pcolormesh(lon_coords, lat_coords, season_data,
#                                   cmap=cmap, shading='auto')
           
#             # Add colorbar
#             cbar = plt.colorbar(im, ax=axes[i], pad=0.05)
#             cbar.set_label(f'Average {label}', fontsize=10)
           
#             axes[i].set_title(f'{season} Average {var_name.upper()}', fontsize=12)
#             axes[i].set_xlabel('Longitude', fontsize=10)
#             axes[i].set_ylabel('Latitude', fontsize=10)
   
#     plt.suptitle(f'Seasonal Average {var_name.upper()}', fontsize=16)
#     plt.tight_layout(rect=[0, 0, 1, 0.96])
#     plt.show()
#     return fig


# # In[64]:


# fig_t2m_seasonal = plot_seasonal_averages(ds_ecmwf, 't2m')


# # In[ ]:




#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd # Used for pd.to_datetime if time conversion is tricky
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

def plot_2m_temp_seasonal_avg_heatmaps(ecmwf_file_path_str, nc_variable_name, output_path_str):
    """
    Plot seasonal average heatmaps for 2m temperature.
   
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
        actual_time_coord_name = None # To store the name of the time coordinate/dimension
        try:
            ds = xr.open_dataset(ecmwf_file_path, decode_times=True)
            # Try to identify the time coordinate name
            if 'time' in ds.coords: actual_time_coord_name = 'time'
            elif 'valid_time' in ds.coords: actual_time_coord_name = 'valid_time'
            # Add other common time coord names if necessary
            
            if not actual_time_coord_name:
                 print(f"Error: Could not identify a primary time coordinate ('time' or 'valid_time') after initial load in {ecmwf_file_path}", file=sys.stderr)
                 sys.exit(1)
            if actual_time_coord_name == 'valid_time': # If primary time is valid_time, rename it to 'time'
                ds = ds.rename({'valid_time': 'time'})
                actual_time_coord_name = 'time' # Now it's 'time'

        except Exception as e:
            print(f"Warning: Initial time decoding failed: {e}. Trying with decode_times=False.", file=sys.stderr)
            try:
                ds = xr.open_dataset(ecmwf_file_path, decode_times=False)
                
                if 'time' in ds.coords: actual_time_coord_name = 'time'
                elif 'valid_time' in ds.coords: actual_time_coord_name = 'valid_time'
                else:
                    print(f"Error: Could not find a recognizable time coordinate ('time', 'valid_time') in {ecmwf_file_path} during fallback.", file=sys.stderr)
                    sys.exit(1)
                
                # Create a new dataset with just the time coordinate for decode_cf
                time_ds_to_decode = xr.Dataset({actual_time_coord_name: ds[actual_time_coord_name]})
                decoded_time_ds = xr.decode_cf(time_ds_to_decode)
                ds[actual_time_coord_name] = decoded_time_ds[actual_time_coord_name] # Update the coordinate in the main ds

                if actual_time_coord_name == 'valid_time': # If it was 'valid_time', rename it to 'time'
                    ds = ds.rename({'valid_time': 'time'})
                    actual_time_coord_name = 'time' # The consistent name is now 'time'
            
            except Exception as e_fallback:
                print(f"Error: Failed to load dataset {ecmwf_file_path} even with fallback: {e_fallback}", file=sys.stderr)
                sys.exit(1)
        
        # At this point, actual_time_coord_name should be 'time' if renaming occurred,
        # or the original name if it was already 'time'.
        # The dataset ds should have this coordinate correctly decoded.

        if nc_variable_name not in ds:
            print(f"Error: Variable '{nc_variable_name}' (expected for 2m Temp) not found in dataset {ecmwf_file_path}", file=sys.stderr)
            print(f"Available variables: {list(ds.data_vars)}", file=sys.stderr)
            sys.exit(1)
        
        data_var = ds[nc_variable_name]

        # Ensure the identified time coordinate is a dimension of the data variable
        if actual_time_coord_name not in data_var.dims:
            # This can happen if 'time' is a coordinate but not a dimension of this specific variable
            # Or if the dimension has a different name even if the coordinate is 'time'
            # Let's try to find a dimension that looks like time
            potential_time_dims = [d for d in data_var.dims if actual_time_coord_name in d or 'time' in d.lower()]
            if not potential_time_dims:
                print(f"Error: The time coordinate '{actual_time_coord_name}' is not a dimension of variable '{nc_variable_name}'. Dims: {data_var.dims}", file=sys.stderr)
                sys.exit(1)
            # If multiple, pick the one that exactly matches actual_time_coord_name, or the first one.
            dim_to_use_for_time = actual_time_coord_name if actual_time_coord_name in potential_time_dims else potential_time_dims[0]
            print(f"Using dimension '{dim_to_use_for_time}' for time operations on variable '{nc_variable_name}'.", file=sys.stdout)
        else:
            dim_to_use_for_time = actual_time_coord_name


        # Define seasons based on month of the time coordinate
        def get_season(time_data_array, time_dim_name_for_data_array):
            month = time_data_array.dt.month
            seasons_map = {
                1: 'DJF', 2: 'DJF', 3: 'MAM', 4: 'MAM', 5: 'MAM',
                6: 'JJA', 7: 'JJA', 8: 'JJA', 9: 'SON', 10: 'SON',
                11: 'SON', 12: 'DJF'
            }
            # Ensure the coords for the new season DataArray match the input time_data_array's time dimension
            return xr.DataArray([seasons_map[m] for m in month.data], coords={time_dim_name_for_data_array: time_data_array[time_dim_name_for_data_array]}, name='season')

        # Use the actual_time_coord_name to access the coordinate from ds for season assignment
        # And dim_to_use_for_time for the dimension in that coordinate array.
        ds_with_season = ds.assign_coords(season=get_season(ds[actual_time_coord_name], dim_to_use_for_time))

        # Calculate seasonal averages using the correct time dimension name for the groupby
        seasonal_avg_data_var = ds_with_season[nc_variable_name].groupby('season').mean(dim=dim_to_use_for_time, skipna=True)
        
        # Unit conversion for 2m temperature
        display_units = 'K'
        cmap = 'RdYlBu_r' # Standard for temperature
        label_text = 'Temperature'

        if not seasonal_avg_data_var.season.size == 0 : 
            sample_mean_val = seasonal_avg_data_var.isel(season=0).mean(skipna=True).item()
            if np.isfinite(sample_mean_val) and sample_mean_val > 100: 
                seasonal_avg_data_var = seasonal_avg_data_var - 273.15
                display_units = '°C'
                print(f"Data for '{nc_variable_name}' appears to be in Kelvin, converted to Celsius for seasonal plots.", file=sys.stdout)
            else:
                original_units = data_var.attrs.get('units', 'K')
                display_units = original_units
                print(f"Data for '{nc_variable_name}' presented in original units: {display_units} (or already Celsius-like values).", file=sys.stdout)
        else:
            print(f"Warning: No seasons could be grouped for variable '{nc_variable_name}'. Plot may be empty or incorrect.", file=sys.stderr)
        
        cbar_label = f'Average {label_text} ({display_units})'

        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 12)) 
        axes = axes.flatten()
        season_order = ['DJF', 'MAM', 'JJA', 'SON']
        
        lat_coord_name = 'latitude' if 'latitude' in ds.coords else 'lat'
        lon_coord_name = 'longitude' if 'longitude' in ds.coords else 'lon'

        if lat_coord_name not in ds.coords or lon_coord_name not in ds.coords:
            print(f"Error: Could not determine lat/lon coordinates. Found: {list(ds.coords)}", file=sys.stderr)
            sys.exit(1)

        lon_coords = ds[lon_coord_name]
        lat_coords = ds[lat_coord_name]

        for i, season_name in enumerate(season_order):
            ax = axes[i]
            if season_name in seasonal_avg_data_var.season.values:
                season_data_to_plot = seasonal_avg_data_var.sel(season=season_name)
                im = ax.pcolormesh(lon_coords.data, lat_coords.data, season_data_to_plot.data,
                                      cmap=cmap, shading='auto')
                cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.05, aspect=25)
                cbar.set_label(cbar_label, fontsize=10)
                ax.set_title(f'{season_name} Average', fontsize=12)
            else:
                ax.set_title(f'{season_name} (No Data)', fontsize=12)
                ax.text(0.5, 0.5, 'No data for this season', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

            ax.set_xlabel('Longitude', fontsize=10)
            ax.set_ylabel('Latitude', fontsize=10)
            ax.tick_params(axis='both', which='major', labelsize=8)

        plt.suptitle(f'Seasonal Average 2m Temperature', fontsize=16) 
        plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
        
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
    parser = argparse.ArgumentParser(description="Generate seasonal average heatmaps for 2m temperature.")
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

    plot_2m_temp_seasonal_avg_heatmaps(ecmwf_file, nc_var_name, args.output)