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


# # In[26]:


# # ECMWF Dataset - Each variable separately (averaged over time)
# def plot_ecmwf_variable_average(ds_ecmwf, var_name, cmap=None):
#     """
#     Plot average ECMWF variables as heatmaps over entire time range
#     var_name options: 'msl', 'tcc', 'u10', 'v10', 't2m'
#     """
#     # Set default colormap if not provided
#     if cmap is None:
#         cmaps = {
#             'msl': 'viridis',
#             'tcc': 'Blues',
#             'u10': 'RdBu_r',
#             'v10': 'RdBu_r',
#             't2m': 'RdYlBu_r'
#         }
#         cmap = cmaps.get(var_name, 'viridis')
   
#     # Calculate time average
#     var_avg = ds_ecmwf[var_name].mean(dim='valid_time')
   
#     # Get variable-specific title and colorbar label
#     var_labels = {
#         'msl': ('Mean Sea Level Pressure', 'Pressure (hPa)'),
#         'tcc': ('Total Cloud Cover', 'Cloud Cover (%)'),
#         'u10': ('U-Component of Wind at 10m', 'Wind Speed (m/s)'),
#         'v10': ('V-Component of Wind at 10m', 'Wind Speed (m/s)'),
#         't2m': ('Temperature at 2m', 'Temperature (K)')
#     }
#     title, cbar_label = var_labels.get(var_name, (var_name, var_name))
   
#     # Create plot
#     fig, ax = plt.subplots(figsize=(10, 8))
   
#     # Create heatmap
#     var_plot = ax.pcolormesh(ds_ecmwf.longitude, ds_ecmwf.latitude, var_avg,
#                             cmap=cmap, shading='auto')
   
#     # Add colorbar and labels
#     cbar = plt.colorbar(var_plot, ax=ax, pad=0.05)
#     cbar.set_label(f'Average {cbar_label}', fontsize=12)
   
#     ax.set_xlabel('Longitude', fontsize=12)
#     ax.set_ylabel('Latitude', fontsize=12)
   
#     # Add time range info to title
#     start_date = np.datetime_as_string(ds_ecmwf.valid_time[0].values, unit="D")
#     end_date = np.datetime_as_string(ds_ecmwf.valid_time[-1].values, unit="D")
#     ax.set_title(f'Average {title}\n{start_date} to {end_date}', fontsize=14)
   
#     plt.tight_layout()
#     plt.show()
#     return fig


# # In[30]:


# for var_name in ['v10']:
#     fig_var_avg = plot_ecmwf_variable_average(ds_ecmwf, var_name)


# # In[ ]:




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
# from matplotlib.colors import LinearSegmentedColormap # Not used
# from matplotlib.ticker import MaxNLocator # Not used

import argparse
import sys
from pathlib import Path

def plot_v10_avg_heatmap(ecmwf_file_path_str, nc_variable_name, output_path_str):
    """
    Plot the time-averaged 10m V-Wind Component (V10) as a heatmap.
   
    Parameters:
    -----------
    ecmwf_file_path_str : str
        Path to the ECMWF data file.
    nc_variable_name : str
        The NetCDF variable name for V10 (e.g., 'v10').
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
            ds_ecmwf = xr.open_dataset(ecmwf_file_path, decode_times=True)
            if 'time' in ds_ecmwf.coords: actual_time_coord_name = 'time'
            elif 'valid_time' in ds_ecmwf.coords: actual_time_coord_name = 'valid_time'
            if not actual_time_coord_name:
                 print(f"Error: Could not identify 'time' or 'valid_time' coordinate after initial load in {ecmwf_file_path}", file=sys.stderr)
                 sys.exit(1)
            if actual_time_coord_name == 'valid_time': 
                ds_ecmwf = ds_ecmwf.rename({'valid_time': 'time'})
                actual_time_coord_name = 'time'
        except Exception as e: # Fallback
            print(f"Warning: Initial time decoding failed: {e}. Trying with decode_times=False.", file=sys.stderr)
            try:
                ds_ecmwf = xr.open_dataset(ecmwf_file_path, decode_times=False)
                if 'time' in ds_ecmwf.coords: actual_time_coord_name = 'time'
                elif 'valid_time' in ds_ecmwf.coords: actual_time_coord_name = 'valid_time'
                else: print(f"Error: No time coord in fallback for {ecmwf_file_path}.", file=sys.stderr); sys.exit(1)
                
                time_ds_to_decode = xr.Dataset({actual_time_coord_name: ds_ecmwf[actual_time_coord_name]})
                decoded_time_ds = xr.decode_cf(time_ds_to_decode)
                ds_ecmwf[actual_time_coord_name] = decoded_time_ds[actual_time_coord_name]
                if actual_time_coord_name == 'valid_time':
                    ds_ecmwf = ds_ecmwf.rename({'valid_time': 'time'})
                    actual_time_coord_name = 'time'
            except Exception as e_fallback: print(f"Error: Failed to load {ecmwf_file_path} with fallback: {e_fallback}", file=sys.stderr); sys.exit(1)

        if nc_variable_name not in ds_ecmwf:
            print(f"Error: Variable '{nc_variable_name}' (expected for V10) not found in dataset {ecmwf_file_path}", file=sys.stderr)
            print(f"Available variables: {list(ds_ecmwf.data_vars)}", file=sys.stderr)
            sys.exit(1)
        
        data_var = ds_ecmwf[nc_variable_name]
        
        dim_to_use_for_time = actual_time_coord_name
        if actual_time_coord_name not in data_var.dims:
            potential_time_dims = [d for d in data_var.dims if actual_time_coord_name in d or 'time' in d.lower()]
            if not potential_time_dims:
                print(f"Error: Time coordinate '{actual_time_coord_name}' is not a dimension of '{nc_variable_name}'. Dims: {data_var.dims}", file=sys.stderr)
                sys.exit(1)
            dim_to_use_for_time = potential_time_dims[0]

        var_avg = data_var.mean(dim=dim_to_use_for_time, skipna=True)
        
        display_units = data_var.attrs.get('units', 'm/s') 
        cbar_label = f'V-Wind ({display_units})'
        plot_title = f'Average 10m V-Wind Component'
        cmap = 'RdBu_r' 

        lat_coord_name = 'latitude' if 'latitude' in ds_ecmwf.coords else 'lat'
        lon_coord_name = 'longitude' if 'longitude' in ds_ecmwf.coords else 'lon'

        if lat_coord_name not in ds_ecmwf.coords or lon_coord_name not in ds_ecmwf.coords:
            print(f"Error: Could not determine lat/lon coordinates. Found coords: {list(ds_ecmwf.coords)}", file=sys.stderr)
            sys.exit(1)
            
        longitude_coords = ds_ecmwf[lon_coord_name].values
        latitude_coords = ds_ecmwf[lat_coord_name].values

        fig, ax = plt.subplots(figsize=(10, 8))
        
        max_abs_val = np.nanmax(np.abs(var_avg.data))
        vmin, vmax = -max_abs_val, max_abs_val

        var_plot = ax.pcolormesh(longitude_coords, latitude_coords, var_avg.data,
                                cmap=cmap, shading='auto', vmin=vmin, vmax=vmax)
       
        cbar = plt.colorbar(var_plot, ax=ax, orientation='vertical', pad=0.05, aspect=30)
        cbar.set_label(f'Average {cbar_label}', fontsize=12)
       
        ax.set_xlabel('Longitude', fontsize=12)
        ax.set_ylabel('Latitude', fontsize=12)
       
        try:
            start_date = pd.to_datetime(str(ds_ecmwf[actual_time_coord_name].min().values)).strftime('%Y-%m-%d')
            end_date = pd.to_datetime(str(ds_ecmwf[actual_time_coord_name].max().values)).strftime('%Y-%m-%d')
            ax.set_title(f'{plot_title}\n({start_date} to {end_date})', fontsize=14)
        except Exception as e_time_title:
            print(f"Warning: Could not determine time range for title: {e_time_title}", file=sys.stderr)
            ax.set_title(plot_title, fontsize=14)

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
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a time-averaged heatmap for 10m V-Wind (V10).")
    parser.add_argument('--files', required=True, nargs=1, help='Path to the input ECMWF NetCDF file.')
    parser.add_argument('--varnames', required=True, nargs=1, help="NetCDF variable name for V10 (e.g., 'v10').")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type (unused by this script).')
    parser.add_argument('--vars', nargs='+', help='Internal variable keys (unused by this script).')
    
    args = parser.parse_args()
        
    ecmwf_file = args.files[0]
    nc_var_name = args.varnames[0]

    plot_v10_avg_heatmap(ecmwf_file, nc_var_name, args.output)