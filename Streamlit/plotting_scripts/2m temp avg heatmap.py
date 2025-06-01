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


# # In[28]:


# for var_name in ['t2m']:
#     fig_var_avg = plot_ecmwf_variable_average(ds_ecmwf, var_name)


# # In[ ]:




#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd # Not directly used
# from matplotlib.dates import DateFormatter # Not used
# import scipy.stats as stats # Not used
import seaborn as sns # For colormaps if needed, though standard matplotlib cmaps are used
# import cartopy.crs as ccrs # Not used for this specific plot type
# import matplotlib.dates as mdates # Not used
# from matplotlib.colors import LinearSegmentedColormap # Not used
# from matplotlib.ticker import MaxNLocator # Not used
import pandas as pd
import argparse
import sys
from pathlib import Path

def plot_2m_temp_avg_heatmap(ecmwf_file_path_str, nc_variable_name, output_path_str):
    """
    Plot the time-averaged 2m temperature as a heatmap.
   
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
        try:
            ds_ecmwf = xr.open_dataset(ecmwf_file_path, decode_times=True)
        except Exception as e:
            print(f"Warning: Initial time decoding failed: {e}. Trying with decode_times=False.", file=sys.stderr)
            try:
                ds_ecmwf = xr.open_dataset(ecmwf_file_path, decode_times=False)
                # Attempt to decode/use common time coordinates
                time_coord_name = None
                if 'time' in ds_ecmwf.coords: time_coord_name = 'time'
                elif 'valid_time' in ds_ecmwf.coords: time_coord_name = 'valid_time'
                
                if time_coord_name:
                    # Create a new dataset with just the time coordinate for decode_cf
                    time_ds = xr.Dataset({time_coord_name: ds_ecmwf[time_coord_name]})
                    decoded_time_ds = xr.decode_cf(time_ds)
                    ds_ecmwf[time_coord_name] = decoded_time_ds[time_coord_name]
                    if time_coord_name == 'valid_time': # Rename for consistency if needed
                         ds_ecmwf = ds_ecmwf.rename({'valid_time': 'time'})
                    time_coord_name = 'time' # Ensure it's 'time' now
                else:
                    print(f"Error: Could not find or decode a recognizable time coordinate in {ecmwf_file_path}", file=sys.stderr)
                    sys.exit(1)
            except Exception as e_fallback:
                print(f"Error: Failed to load dataset {ecmwf_file_path} even with fallback: {e_fallback}", file=sys.stderr)
                sys.exit(1)

        if nc_variable_name not in ds_ecmwf:
            print(f"Error: Variable '{nc_variable_name}' (expected for 2m Temp) not found in dataset {ecmwf_file_path}", file=sys.stderr)
            print(f"Available variables: {list(ds_ecmwf.data_vars)}", file=sys.stderr)
            sys.exit(1)
        
        data_var = ds_ecmwf[nc_variable_name]
        time_dim_name = 'time' # Assumed after potential rename above

        if time_dim_name not in data_var.dims:
             # Fallback if 'time' is not a dimension of the data variable itself,
             # but present in ds_ecmwf.coords (e.g. scalar time or other complex structure)
             # This would be unusual for typical gridded data.
            potential_time_dims = [d for d in data_var.dims if 'time' in d.lower()]
            if not potential_time_dims:
                print(f"Error: Time dimension ('{time_dim_name}' or similar) not found in variable '{nc_variable_name}'. Dims: {data_var.dims}", file=sys.stderr)
                sys.exit(1)
            time_dim_name = potential_time_dims[0]


        # Calculate time average
        var_avg = data_var.mean(dim=time_dim_name, skipna=True)
        
        # Unit handling for 2m temperature
        current_units = data_var.attrs.get('units', 'K').lower() # Default to K, convert to lower for comparison
        display_units = 'K'
        cbar_label_val = 'Temperature'

        if 'k' in current_units and var_avg.mean(skipna=True).item() > 100: # If units contain 'k' (Kelvin) and mean is high
            var_avg = var_avg - 273.15  # Convert to Celsius
            display_units = '°C'
            print(f"Data for '{nc_variable_name}' appears to be in Kelvin, converted to Celsius.", file=sys.stdout)
        elif 'c' in current_units: # If units explicitly Celsius
             display_units = '°C'
        else: # Otherwise, assume Kelvin or keep original if not explicitly K/C and not converted
            display_units = data_var.attrs.get('units', 'K') # Use original attribute if available
            print(f"Data for '{nc_variable_name}' presented in original units: {display_units} (or already Celsius-like values).", file=sys.stdout)

        cbar_label = f'{cbar_label_val} ({display_units})'
        plot_title = f'Average 2m Temperature'
        cmap = 'RdYlBu_r' # Common for temperature

        # Determine spatial coordinate names
        lat_coord_name = 'latitude' if 'latitude' in ds_ecmwf.coords else 'lat'
        lon_coord_name = 'longitude' if 'longitude' in ds_ecmwf.coords else 'lon'

        if lat_coord_name not in ds_ecmwf.coords or lon_coord_name not in ds_ecmwf.coords:
            print(f"Error: Could not determine lat/lon coordinates. Looked for 'latitude'/'lat' and 'longitude'/'lon'. Found coords: {list(ds_ecmwf.coords)}", file=sys.stderr)
            sys.exit(1)
            
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 8))
       
        # Create heatmap
        # Ensure longitude and latitude are 1D for pcolormesh if they come from the var_avg itself
        # Or use ds_ecmwf.longitude, ds_ecmwf.latitude if they are guaranteed to be 1D grid coords
        
        # Assuming ds_ecmwf.longitude and ds_ecmwf.latitude are the 1D coordinate arrays for the grid
        longitude_coords = ds_ecmwf[lon_coord_name].values
        latitude_coords = ds_ecmwf[lat_coord_name].values

        var_plot = ax.pcolormesh(longitude_coords, latitude_coords, var_avg.data,
                                cmap=cmap, shading='auto') # Use .data to pass numpy array
       
        cbar = plt.colorbar(var_plot, ax=ax, orientation='vertical', pad=0.05, aspect=30)
        cbar.set_label(f'Average {cbar_label}', fontsize=12)
       
        ax.set_xlabel('Longitude', fontsize=12)
        ax.set_ylabel('Latitude', fontsize=12)
       
        try:
            # Add time range info to title if time coordinate is valid
            start_date = pd.to_datetime(str(ds_ecmwf[time_dim_name].min().values)).strftime('%Y-%m-%d')
            end_date = pd.to_datetime(str(ds_ecmwf[time_dim_name].max().values)).strftime('%Y-%m-%d')
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
    parser = argparse.ArgumentParser(description="Generate a time-averaged heatmap for 2m temperature.")
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

    plot_2m_temp_avg_heatmap(ecmwf_file, nc_var_name, args.output)