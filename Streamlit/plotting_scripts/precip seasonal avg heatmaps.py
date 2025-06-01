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


# # In[25]:


# # Load your datasets
# ds_precip = xr.open_dataset( r"precip_konkan.nc")
# ds_ecmwf = xr.open_dataset( r"ecmwf datas.nc")
# ds_sst = xr.open_dataset( r"noaa_sst_masked.nc")


# # In[43]:


# # Seasonal Averages 
# def plot_seasonal_averages(ds, var_name, time_dim='time', lat_dim='lat', lon_dim='lon'):
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
#     ds_with_season = ds_with_season.assign_coords(season=('time', seasons))
   
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


# # In[44]:


# fig_precip_seasonal = plot_seasonal_averages(ds_precip, 'precip')

#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd # For pd.to_datetime if time conversion is tricky
# from matplotlib.dates import DateFormatter # Not used
# import scipy.stats as stats # Not used
# import seaborn as sns # Not directly used, standard cmaps are fine
# import cartopy.crs as ccrs # Not used
# import matplotlib.dates as mdates # Not used
from matplotlib.colors import LinearSegmentedColormap # If using custom cmap
# from matplotlib.ticker import MaxNLocator # Not used

import argparse
import sys
from pathlib import Path

def plot_precip_seasonal_avg_heatmaps(noaa_precip_file_path_str, nc_variable_name, output_path_str):
    """
    Plot seasonal average heatmaps for NOAA precipitation.
   
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
            print(f"Using dimension '{dim_to_use_for_time}' for time operations on variable '{nc_variable_name}'.", file=sys.stdout)

        # Define seasons
        def get_season(time_coord_da, time_dim_for_da):
            month = time_coord_da.dt.month
            seasons_map = {
                1: 'DJF', 2: 'DJF', 3: 'MAM', 4: 'MAM', 5: 'MAM',
                6: 'JJA', 7: 'JJA', 8: 'JJA', 9: 'SON', 10: 'SON',
                11: 'SON', 12: 'DJF'
            }
            return xr.DataArray([seasons_map[m] for m in month.data], coords={time_dim_for_da: time_coord_da[time_dim_for_da]}, name='season')

        ds_with_season = ds_precip.assign_coords(season=get_season(ds_precip[actual_time_coord_name], dim_to_use_for_time))
        seasonal_avg_data_var = ds_with_season[nc_variable_name].groupby('season').mean(dim=dim_to_use_for_time, skipna=True)
        
        # Units and Colormap for Precipitation
        display_units = data_var.attrs.get('units', 'mm') # Default to mm
        cmap = 'Blues' # Standard sequential colormap for precipitation
        # Or your custom one:
        # cmap = LinearSegmentedColormap.from_list('precip_blues',
        #                                         ['#FFFFFF', '#E0F2F7', '#B3E0F2', '#80CEEB', 
        #                                          '#4DBCE6', '#1AA9E0', '#008BC7', '#006A9B', '#004B70'])
        label_text = 'Precipitation'
        cbar_label = f'Average {label_text} ({display_units})'

        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 12))
        axes = axes.flatten()
        season_order = ['DJF', 'MAM', 'JJA', 'SON']
        
        lat_coord_name = 'lat' if 'lat' in ds_precip.coords else 'latitude'
        lon_coord_name = 'lon' if 'lon' in ds_precip.coords else 'longitude'

        if lat_coord_name not in ds_precip.coords or lon_coord_name not in ds_precip.coords:
            print(f"Error: Could not determine lat/lon coordinates. Found: {list(ds_precip.coords)}", file=sys.stderr)
            sys.exit(1)

        lon_coords = ds_precip[lon_coord_name]
        lat_coords = ds_precip[lat_coord_name]

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

        plt.suptitle(f'Seasonal Average Precipitation', fontsize=16)
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
    parser = argparse.ArgumentParser(description="Generate seasonal average heatmaps for NOAA precipitation.")
    parser.add_argument('--files', required=True, nargs=1, help='Path to the input NOAA precipitation NetCDF file.')
    parser.add_argument('--varnames', required=True, nargs=1, help="NetCDF variable name for precipitation (e.g., 'precip').")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type (unused by this script).')
    parser.add_argument('--vars', nargs='+', help='Internal variable keys (unused by this script).')
    
    args = parser.parse_args()
        
    noaa_file = args.files[0]
    nc_var_name = args.varnames[0]

    plot_precip_seasonal_avg_heatmaps(noaa_file, nc_var_name, args.output)