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


# # In[11]:


# def plot_noaa_precip_grid_points(precip_file):
#     """
#     Create timeseries plot for NOAA Precipitation data - Grid points only
#     """
#     # Load the dataset
#     ds = xr.open_dataset(precip_file)
   
#     # Create a figure for grid points precipitation
#     plt.figure(figsize=(12, 5))
   
#     # Plot: Precipitation at specific points
#     for i, lat in enumerate(ds.lat.values):
#         for j, lon in enumerate(ds.lon.values):
#             if i == 0 and j == 0:  # Only plot a selection of grid points
#                 label = f"lat={lat:.2f}, lon={lon:.2f}"
#                 ds.precip.sel(lat=lat, lon=lon).plot(label=label)
   
#     plt.title('NOAA Precipitation at Selected Grid Points', fontsize=14)
#     plt.xlabel('Time', fontsize=12)
#     plt.ylabel('Precipitation (mm/day)', fontsize=12)
#     plt.legend(loc='upper right', fontsize=10)
#     plt.gca().xaxis.set_major_locator(mdates.YearLocator())
#     plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
   
#     plt.tight_layout()
#     #plt.savefig('noaa_precip_grid_points.png', dpi=300)
#     plt.show()
#     plt.close()
   


# # In[12]:


# if __name__ == "__main__":
#     # Define your file paths here
#     noaa_precip_file = r"precip_konkan.nc"

#     # Generate individual dataset plots
#     plot_noaa_precip_grid_points(noaa_precip_file)


# # In[ ]:




#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd # For pd.to_datetime if time conversion is tricky
# from matplotlib.dates import DateFormatter # mdates.DateFormatter is used
# import scipy.stats as stats # Not used
# import seaborn as sns # Not directly used for this plot
# import cartopy.crs as ccrs # Not used
import matplotlib.dates as mdates # For x-axis formatting
# from matplotlib.colors import LinearSegmentedColormap # Not used
# from matplotlib.ticker import MaxNLocator # Not used

import argparse
import sys
from pathlib import Path

def plot_precip_at_selected_grid_points(noaa_precip_file_path_str, nc_variable_name, output_path_str, num_points_to_plot=1):
    """
    Create a timeseries plot for NOAA Precipitation at a few selected grid points.
   
    Parameters:
    -----------
    noaa_precip_file_path_str : str
        Path to the NOAA precipitation data file.
    nc_variable_name : str
        The NetCDF variable name for precipitation (e.g., 'precip').
    output_path_str : str
        Path to save the output plot image.
    num_points_to_plot : int
        Number of distinct grid points to plot (selected from the start of lat/lon arrays).
        To keep the plot clean, this is usually a small number.
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
        
        # Determine spatial coordinate names (typically 'lat', 'lon' for NOAA CPC data)
        lat_coord_name = 'lat' if 'lat' in data_var.dims else 'latitude'
        lon_coord_name = 'lon' if 'lon' in data_var.dims else 'longitude'

        if lat_coord_name not in data_var.dims or lon_coord_name not in data_var.dims:
            print(f"Error: Could not determine lat/lon dimensions for '{nc_variable_name}'. Dims: {data_var.dims}", file=sys.stderr)
            sys.exit(1)

        # Select grid points to plot
        # For simplicity, select the first few points or points from the middle of the grid.
        # Your original script selected only the first point (0,0).
        # Let's select `num_points_to_plot` distinct points, spread out if possible or just the first few.
        
        lat_values = ds_precip[lat_coord_name].values
        lon_values = ds_precip[lon_coord_name].values
        
        points_plotted = 0
        
        plt.figure(figsize=(12, 6)) # Adjusted figure size slightly

        # Plot a few selected points to avoid overcrowding
        # Example: plot num_points_to_plot, taking them somewhat spread out
        # or just the first num_points_to_plot if the grid is small
        
        # Simple selection: first `num_points_to_plot` points along one dimension, or a diagonal
        # For a cleaner plot, let's pick up to `num_points_to_plot` distinct points
        # by iterating through lat and then lon.
        
        # If num_points_to_plot = 1, it plots the very first grid cell ds.lat[0], ds.lon[0]
        # If num_points_to_plot = 2, it plots ds.lat[0],ds.lon[0] and ds.lat[0],ds.lon[1] (if lon has >1 point)
        # etc. This tries to take points from the "corner" of the grid.

        # A slightly better strategy for a few points might be to take corners or center.
        # For now, let's stick to a simple configurable number of points from the start.
        
        # Select up to `num_points_to_plot` grid cells
        # Example: take the first point, a point from middle lat / first lon, first lat / middle lon
        
        indices_to_plot = []
        if len(lat_values) > 0 and len(lon_values) > 0:
            indices_to_plot.append({'lat_idx': 0, 'lon_idx': 0, 'label_suffix': ' (Corner 1)'}) # First point
            if num_points_to_plot > 1 and len(lat_values) > 1 and len(lon_values) > 1:
                indices_to_plot.append({'lat_idx': len(lat_values)-1, 'lon_idx': len(lon_values)-1, 'label_suffix': ' (Corner 2)'}) # Opposite corner
            if num_points_to_plot > 2 and len(lat_values) > 2 and len(lon_values) > 2: # Approx middle
                 indices_to_plot.append({'lat_idx': len(lat_values)//2, 'lon_idx': len(lon_values)//2, 'label_suffix': ' (Center)'})
        
        # Ensure we don't try to plot more than requested or available points through this selection
        actual_indices_to_plot = indices_to_plot[:min(num_points_to_plot, len(indices_to_plot))]
        if not actual_indices_to_plot and len(lat_values) > 0 and len(lon_values) > 0 : # Fallback if num_points_to_plot was 0 or logic failed
            actual_indices_to_plot = [{'lat_idx': 0, 'lon_idx': 0, 'label_suffix': ''}]


        for point_info in actual_indices_to_plot:
            lat_val = lat_values[point_info['lat_idx']]
            lon_val = lon_values[point_info['lon_idx']]
            label = f"Lat={lat_val:.2f}, Lon={lon_val:.2f}{point_info.get('label_suffix','')}"
            
            # Use .sel for selection by coordinate values
            data_var.sel({lat_coord_name: lat_val, lon_coord_name: lon_val}).plot(ax=plt.gca(), label=label, linewidth=1.2)
            points_plotted += 1

        if points_plotted == 0:
            print("Warning: No grid points were selected for plotting. The plot will be empty.", file=sys.stderr)
            # Optionally, plot area average as a fallback or just show an empty plot with a message
            plt.text(0.5, 0.5, "No grid points selected/plotted", horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)

        display_units = data_var.attrs.get('units', 'mm/day')
        plot_title = f'Precipitation at Selected Grid Points'
        
        plt.title(plot_title, fontsize=14)
        plt.xlabel('Time', fontsize=12)
        plt.ylabel(f'Precipitation ({display_units})', fontsize=12)
        
        if points_plotted > 0:
            plt.legend(loc='best', fontsize=9) # 'best' might be better than 'upper right'
        
        plt.gca().xaxis.set_major_locator(mdates.YearLocator())
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        plt.gca().tick_params(axis='x', rotation=45)
       
        plt.grid(True, linestyle='--', alpha=0.7)
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
    parser = argparse.ArgumentParser(description="Generate a timeseries plot for NOAA precipitation at selected grid points.")
    parser.add_argument('--files', required=True, nargs=1, help='Path to the input NOAA precipitation NetCDF file.')
    parser.add_argument('--varnames', required=True, nargs=1, help="NetCDF variable name for precipitation (e.g., 'precip').")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type (unused by this script).')
    parser.add_argument('--vars', nargs='+', help='Internal variable keys (unused by this script).')
    parser.add_argument('--num_points', type=int, default=3, help='Number of distinct grid points to attempt to plot (default: 3). Max 3 from corners/center.')


    args = parser.parse_args()
        
    noaa_file = args.files[0]
    nc_var_name = args.varnames[0]

    plot_precip_at_selected_grid_points(noaa_file, nc_var_name, args.output, num_points_to_plot=args.num_points)