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


# # In[33]:


# # ECMWF Wind Vector Plot (Average U10 and V10)
# def plot_ecmwf_wind_vectors_average(ds_ecmwf):
#     """
#     Plot average wind vectors using u10 and v10 components over entire time range
#     """
#     # Calculate time average
#     u10_avg = ds_ecmwf.u10.mean(dim='valid_time')
#     v10_avg = ds_ecmwf.v10.mean(dim='valid_time')
#     wind_speed_avg = np.sqrt(u10_avg**2 + v10_avg**2)
   
#     # Create plot with wind vectors and wind speed as background
#     fig, ax = plt.subplots(figsize=(12, 10))
   
#     # Plot wind speed as background
#     speed_plot = ax.pcolormesh(ds_ecmwf.longitude, ds_ecmwf.latitude, wind_speed_avg,
#                               cmap='viridis', shading='auto')
   
#     # Add wind vectors (subsample for clarity)
#     quiver_stride = 1  # adjust for desired density (smaller stride for higher density)
#     lon_mesh, lat_mesh = np.meshgrid(ds_ecmwf.longitude[::quiver_stride], ds_ecmwf.latitude[::quiver_stride])
#     ax.quiver(lon_mesh, lat_mesh,
#              u10_avg[::quiver_stride, ::quiver_stride],
#              v10_avg[::quiver_stride, ::quiver_stride],
#              scale=30, color='white', alpha=0.8)
   
#     # Add colorbar and labels
#     cbar = plt.colorbar(speed_plot, ax=ax, pad=0.05)
#     cbar.set_label('Average Wind Speed (m/s)', fontsize=12)
   
#     ax.set_xlabel('Longitude', fontsize=12)
#     ax.set_ylabel('Latitude', fontsize=12)
   
#     # Add time range info to title
#     start_date = np.datetime_as_string(ds_ecmwf.valid_time[0].values, unit="D")
#     end_date = np.datetime_as_string(ds_ecmwf.valid_time[-1].values, unit="D")
#     ax.set_title(f'Average Wind Vectors\n{start_date} to {end_date}', fontsize=14)
   
#     plt.tight_layout()
#     plt.show()
#     return fig


# # In[34]:


# fig_wind_avg = plot_ecmwf_wind_vectors_average(ds_ecmwf)


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

def plot_wind_vectors_avg_heatmap(ecmwf_file_path_str, u10_nc_var_name, v10_nc_var_name, output_path_str):
    """
    Plot time-averaged wind vectors (from u10, v10) with wind speed as background heatmap.
   
    Parameters:
    -----------
    ecmwf_file_path_str : str
        Path to the ECMWF data file.
    u10_nc_var_name : str
        NetCDF variable name for the 10m U-component of wind.
    v10_nc_var_name : str
        NetCDF variable name for the 10m V-component of wind.
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

        if u10_nc_var_name not in ds_ecmwf:
            print(f"Error: U10 variable '{u10_nc_var_name}' not found in {ecmwf_file_path}", file=sys.stderr); sys.exit(1)
        if v10_nc_var_name not in ds_ecmwf:
            print(f"Error: V10 variable '{v10_nc_var_name}' not found in {ecmwf_file_path}", file=sys.stderr); sys.exit(1)
            
        u10_data_var = ds_ecmwf[u10_nc_var_name]
        v10_data_var = ds_ecmwf[v10_nc_var_name]
        
        # Ensure the identified time coordinate is a dimension of the data variables
        dim_to_use_for_time_u = actual_time_coord_name if actual_time_coord_name in u10_data_var.dims else [d for d in u10_data_var.dims if 'time' in d.lower()][0]
        dim_to_use_for_time_v = actual_time_coord_name if actual_time_coord_name in v10_data_var.dims else [d for d in v10_data_var.dims if 'time' in d.lower()][0]
        if not dim_to_use_for_time_u or not dim_to_use_for_time_v :
             print(f"Error: Could not determine time dimension for U10 or V10.", file=sys.stderr); sys.exit(1)


        # Calculate time average
        u10_avg = u10_data_var.mean(dim=dim_to_use_for_time_u, skipna=True)
        v10_avg = v10_data_var.mean(dim=dim_to_use_for_time_v, skipna=True)
        wind_speed_avg = np.sqrt(u10_avg**2 + v10_avg**2)
       
        # Determine spatial coordinate names
        lat_coord_name = 'latitude' if 'latitude' in ds_ecmwf.coords else 'lat'
        lon_coord_name = 'longitude' if 'longitude' in ds_ecmwf.coords else 'lon'

        if lat_coord_name not in ds_ecmwf.coords or lon_coord_name not in ds_ecmwf.coords:
            print(f"Error: Could not determine lat/lon coordinates. Found: {list(ds_ecmwf.coords)}", file=sys.stderr)
            sys.exit(1)
            
        longitude_coords = ds_ecmwf[lon_coord_name].values
        latitude_coords = ds_ecmwf[lat_coord_name].values

        fig, ax = plt.subplots(figsize=(12, 10))
       
        # Plot wind speed as background
        speed_cmap = 'GnBu' # Good for wind speed magnitude
        speed_plot = ax.pcolormesh(longitude_coords, latitude_coords, wind_speed_avg.data,
                                  cmap=speed_cmap, shading='auto')
       
        # Add wind vectors (subsample for clarity)
        # Dynamic stride to avoid overcrowding, ensuring at least 1 if grid is too small
        quiver_stride_lon = max(1, len(longitude_coords) // 20) 
        quiver_stride_lat = max(1, len(latitude_coords) // 20)
        
        # Meshgrid for quiver, using the original full coordinates for indexing
        lon_mesh_full, lat_mesh_full = np.meshgrid(longitude_coords, latitude_coords)
        
        ax.quiver(lon_mesh_full[::quiver_stride_lat, ::quiver_stride_lon], 
                  lat_mesh_full[::quiver_stride_lat, ::quiver_stride_lon],
                  u10_avg.data[::quiver_stride_lat, ::quiver_stride_lon],
                  v10_avg.data[::quiver_stride_lat, ::quiver_stride_lon],
                  scale=None, # Auto-scale based on vector magnitudes
                  scale_units='inches', # Scale relative to plot size in inches
                  headwidth=4, headlength=5, width=0.004, # Adjust arrow appearance
                  color='black', alpha=0.8, pivot='middle')
       
        cbar = plt.colorbar(speed_plot, ax=ax, orientation='vertical', pad=0.05, aspect=30)
        cbar.set_label('Average Wind Speed (m/s)', fontsize=12)
       
        ax.set_xlabel('Longitude', fontsize=12)
        ax.set_ylabel('Latitude', fontsize=12)
       
        try:
            start_date = pd.to_datetime(str(ds_ecmwf[actual_time_coord_name].min().values)).strftime('%Y-%m-%d')
            end_date = pd.to_datetime(str(ds_ecmwf[actual_time_coord_name].max().values)).strftime('%Y-%m-%d')
            ax.set_title(f'Average Wind Speed and Vectors\n({start_date} to {end_date})', fontsize=14)
        except Exception as e_time_title:
            print(f"Warning: Could not determine time range for title: {e_time_title}", file=sys.stderr)
            ax.set_title(f'Average Wind Speed and Vectors', fontsize=14)

        plt.tight_layout(pad=1.5)
        
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
    parser = argparse.ArgumentParser(description="Generate a time-averaged wind vector plot with speed heatmap.")
    parser.add_argument('--files', required=True, nargs=1, 
                        help='Path to the input ECMWF NetCDF file (containing u10 and v10).')
    parser.add_argument('--varnames', required=True, nargs=2, 
                        help="NetCDF variable names for U10 and V10 (e.g., 'u10' 'v10'). Order: U then V.")
    parser.add_argument('--vars', required=True, nargs=2, 
                        help="Internal variable keys (e.g., 'ecmwf_u10' 'ecmwf_v10'). Order should match --varnames.")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type (unused by this script).')
    
    args = parser.parse_args()

    if len(args.vars) != 2 or len(args.varnames) != 2:
        print("Error: Requires exactly two internal keys (--vars) and two NetCDF variable names (--varnames) for u10 and v10.", file=sys.stderr)
        sys.exit(1)

    # Determine u10 and v10 nc_names based on internal keys from --vars
    u10_nc_name_arg, v10_nc_name_arg = None, None
    if 'u10' in args.vars[0].lower() and 'v10' in args.vars[1].lower():
        u10_nc_name_arg = args.varnames[0]
        v10_nc_name_arg = args.varnames[1]
    elif 'v10' in args.vars[0].lower() and 'u10' in args.vars[1].lower():
        v10_nc_name_arg = args.varnames[0] # This was u10's varname
        u10_nc_name_arg = args.varnames[1] # This was v10's varname
    else:
        print("Error: Could not reliably determine U10 and V10 from --vars. Ensure one key contains 'u10' and the other 'v10'.", file=sys.stderr)
        sys.exit(1)
        
    ecmwf_file = args.files[0] # Expects a single file containing both u10 and v10

    plot_wind_vectors_avg_heatmap(ecmwf_file, u10_nc_name_arg, v10_nc_name_arg, args.output)