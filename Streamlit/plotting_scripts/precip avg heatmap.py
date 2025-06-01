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


# # In[16]:


# # Set plotting style
# plt.style.use('seaborn-whitegrid')

# # 1. NOAA Precipitation Dataset - Average over time
# def plot_noaa_precip_average(ds_precip):
#     """
#     Plot average NOAA precipitation heatmap over entire time range
#     """
#     # Calculate time average
#     precip_avg = ds_precip.precip.mean(dim='time')
   
#     # Create plot
#     fig, ax = plt.subplots(figsize=(10, 8))
   
#     # Custom blue colormap for precipitation
#     cmap_precip = LinearSegmentedColormap.from_list('precip_blues',
#                                                     ['#FFFFFF', '#D1EEFC', '#75C6EF', '#1E90FF', '#0066CC', '#003366'])
   
#     # Create heatmap
#     precip_plot = ax.pcolormesh(ds_precip.lon, ds_precip.lat, precip_avg,
#                                cmap=cmap_precip, shading='auto')
   
#     # Add colorbar and labels
#     cbar = plt.colorbar(precip_plot, ax=ax, pad=0.05)
#     cbar.set_label('Average Precipitation (mm)', fontsize=12)
   
#     ax.set_xlabel('Longitude', fontsize=12)
#     ax.set_ylabel('Latitude', fontsize=12)
   
#     # Add time range info to title
#     start_date = np.datetime_as_string(ds_precip.time[0].values, unit="D")
#     end_date = np.datetime_as_string(ds_precip.time[-1].values, unit="D")
#     ax.set_title(f'Average Precipitation\n{start_date} to {end_date}', fontsize=14)
   
#     plt.tight_layout()
#     plt.show()
#     return fig


# # In[17]:


# # Load your datasets
# ds_precip = xr.open_dataset( r"precip_konkan.nc")
# ds_ecmwf = xr.open_dataset( r"ecmwf datas.nc")
# ds_sst = xr.open_dataset( r"noaa_sst_masked.nc")


# # In[18]:


# fig_precip_avg = plot_noaa_precip_average(ds_precip)


# # In[ ]:




#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd # For pd.to_datetime if time conversion is tricky
# from matplotlib.dates import DateFormatter # Not used
# import scipy.stats as stats # Not used
# import seaborn as sns # Not directly used, custom cmap is defined
# import cartopy.crs as ccrs # Not used
# import matplotlib.dates as mdates # Not used
from matplotlib.colors import LinearSegmentedColormap # For custom cmap
# from matplotlib.ticker import MaxNLocator # Not used

import argparse
import sys
from pathlib import Path

def plot_precip_avg_heatmap(noaa_precip_file_path_str, nc_variable_name, output_path_str):
    """
    Plot the time-averaged NOAA precipitation as a heatmap.
   
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
        actual_time_coord_name = None # To store the name of the time coordinate/dimension
        try:
            ds_precip = xr.open_dataset(noaa_precip_file_path, decode_times=True)
            # Try to identify the time coordinate name
            if 'time' in ds_precip.coords: actual_time_coord_name = 'time'
            # Add other common time coord names if necessary (e.g., 'valid_time' if it could be from another source)
            
            if not actual_time_coord_name:
                 print(f"Error: Could not identify a primary time coordinate ('time') after initial load in {noaa_precip_file_path}", file=sys.stderr)
                 sys.exit(1)
            # No rename needed if 'time' is standard for NOAA precip files

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

        # Calculate time average
        precip_avg = data_var.mean(dim=dim_to_use_for_time, skipna=True)
        
        # Define units and colormap
        display_units = data_var.attrs.get('units', 'mm') # Default to mm if not specified
        cbar_label = f'Average Precipitation ({display_units})'
        plot_title = f'Average Precipitation'
        
        # Custom blue colormap for precipitation
        cmap_precip = LinearSegmentedColormap.from_list('precip_blues',
                                                        ['#FFFFFF', '#E0F2F7', '#B3E0F2', '#80CEEB', 
                                                         '#4DBCE6', '#1AA9E0', '#008BC7', '#006A9B', '#004B70'])


        # Determine spatial coordinate names (typically 'lat', 'lon' for NOAA CPC data)
        lat_coord_name = 'lat' if 'lat' in ds_precip.coords else 'latitude'
        lon_coord_name = 'lon' if 'lon' in ds_precip.coords else 'longitude'

        if lat_coord_name not in ds_precip.coords or lon_coord_name not in ds_precip.coords:
            print(f"Error: Could not determine lat/lon coordinates. Looked for 'lat'/'latitude' and 'lon'/'longitude'. Found coords: {list(ds_precip.coords)}", file=sys.stderr)
            sys.exit(1)
            
        longitude_coords = ds_precip[lon_coord_name].values
        latitude_coords = ds_precip[lat_coord_name].values

        fig, ax = plt.subplots(figsize=(10, 8))
        precip_plot = ax.pcolormesh(longitude_coords, latitude_coords, precip_avg.data,
                                   cmap=cmap_precip, shading='auto')
       
        cbar = plt.colorbar(precip_plot, ax=ax, orientation='vertical', pad=0.05, aspect=30)
        cbar.set_label(cbar_label, fontsize=12)
       
        ax.set_xlabel('Longitude', fontsize=12)
        ax.set_ylabel('Latitude', fontsize=12)
       
        try:
            start_date = pd.to_datetime(str(ds_precip[actual_time_coord_name].min().values)).strftime('%Y-%m-%d')
            end_date = pd.to_datetime(str(ds_precip[actual_time_coord_name].max().values)).strftime('%Y-%m-%d')
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
        # import traceback # For debugging
        # traceback.print_exc(file=sys.stderr) # For debugging
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a time-averaged heatmap for NOAA precipitation.")
    parser.add_argument('--files', required=True, nargs=1, help='Path to the input NOAA precipitation NetCDF file.')
    parser.add_argument('--varnames', required=True, nargs=1, help="NetCDF variable name for precipitation (e.g., 'precip').")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type (unused by this script).')
    parser.add_argument('--vars', nargs='+', help='Internal variable keys (unused by this script).')
    
    args = parser.parse_args()
        
    noaa_file = args.files[0]
    nc_var_name = args.varnames[0]

    plot_precip_avg_heatmap(noaa_file, nc_var_name, args.output)