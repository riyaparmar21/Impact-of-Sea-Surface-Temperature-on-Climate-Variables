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


# # In[41]:


# # Multi-Variable ECMWF Panel (all in one figure) - Averaged
# def plot_ecmwf_panel_average(ds_ecmwf):
#     """
#     Plot all averaged ECMWF variables in a single panel
#     """
#     # Define variable names and their associated cmaps
#     vars_and_cmaps = {
#         'msl': 'viridis',
#         'tcc': 'Blues',
#         't2m': 'RdYlBu_r',
#     }
   
#     # Create figure with subplots
#     fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 14))
#     axes = axes.flatten()
   
#     # Plot each variable
#     for i, (var_name, cmap) in enumerate(vars_and_cmaps.items()):
#         if i < len(axes) - 1:  # One spot reserved for wind vectors
#             var_avg = ds_ecmwf[var_name].mean(dim='valid_time')
#             var_plot = axes[i].pcolormesh(ds_ecmwf.longitude, ds_ecmwf.latitude, var_avg,
#                                        cmap=cmap, shading='auto')
#             cbar = plt.colorbar(var_plot, ax=axes[i], pad=0.05)
           
#             # Get variable labels
#             var_labels = {
#                 'msl': ('Mean Sea Level Pressure', 'Pressure (hPa)'),
#                 'tcc': ('Total Cloud Cover', 'Cloud Cover (%)'),
#                 't2m': ('Temperature at 2m', 'Temperature (K)')
#             }
#             title, cbar_label = var_labels.get(var_name, (var_name, var_name))
           
#             cbar.set_label(f'Average {cbar_label}', fontsize=10)
#             axes[i].set_title(f'Average {title}', fontsize=12)
#             axes[i].set_xlabel('Longitude', fontsize=10)
#             axes[i].set_ylabel('Latitude', fontsize=10)
   
#     # Plot average wind vectors in the last panel
#     u10_avg = ds_ecmwf.u10.mean(dim='valid_time')
#     v10_avg = ds_ecmwf.v10.mean(dim='valid_time')
#     wind_speed_avg = np.sqrt(u10_avg**2 + v10_avg**2)
   
#     # Plot wind speed as background
#     speed_plot = axes[-1].pcolormesh(ds_ecmwf.longitude, ds_ecmwf.latitude, wind_speed_avg,
#                                   cmap='viridis', shading='auto')
   
#     # Add wind vectors
#     quiver_stride = 1  # Smaller stride for higher vector density
#     lon_mesh, lat_mesh = np.meshgrid(ds_ecmwf.longitude[::quiver_stride], ds_ecmwf.latitude[::quiver_stride])
#     axes[-1].quiver(lon_mesh, lat_mesh,
#                  u10_avg[::quiver_stride, ::quiver_stride],
#                  v10_avg[::quiver_stride, ::quiver_stride],
#                  scale=30, color='white', alpha=0.8)
   
#     cbar = plt.colorbar(speed_plot, ax=axes[-1], pad=0.05)
#     cbar.set_label('Average Wind Speed (m/s)', fontsize=10)
#     axes[-1].set_title('Average Wind Vectors (u10, v10)', fontsize=12)
#     axes[-1].set_xlabel('Longitude', fontsize=10)
#     axes[-1].set_ylabel('Latitude', fontsize=10)
   
#     # Add time range to main title
#     start_date = np.datetime_as_string(ds_ecmwf.valid_time[0].values, unit="D")
#     end_date = np.datetime_as_string(ds_ecmwf.valid_time[-1].values, unit="D")
#     plt.suptitle(f'ECMWF Average Variables\n{start_date} to {end_date}', fontsize=16)
   
#     plt.tight_layout(rect=[0, 0, 1, 0.96])
#     plt.show()
#     return fig


# # In[42]:


# fig_panel_avg = plot_ecmwf_panel_average(ds_ecmwf)


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

def plot_ecmwf_multipanel_average(ecmwf_file_path_str, 
                                  nc_var_names_map, # Dict: {'msl_key': 'msl_nc', 't2m_key': 't2m_nc', ...}
                                  output_path_str):
    """
    Plot time-averaged MSL, TCC, T2M, and Wind Vectors in a single multi-panel figure.
   
    Parameters:
    -----------
    ecmwf_file_path_str : str
        Path to the ECMWF data file.
    nc_var_names_map : dict
        A dictionary mapping internal variable keys (e.g., 'ecmwf_msl', 'ecmwf_t2m', 
        'ecmwf_u10', 'ecmwf_v10', 'ecmwf_tcc') to their actual NetCDF variable names.
    output_path_str : str
        Path to save the output plot image.
    """
    try:
        ecmwf_file_path = Path(ecmwf_file_path_str)
        output_path = Path(output_path_str)

        if not ecmwf_file_path.exists():
            print(f"Error: ECMWF data file not found at {ecmwf_file_path}", file=sys.stderr)
            sys.exit(1)

        # Define which variables are expected for this multipanel plot and their properties
        # Keys here are simplified (e.g., 'msl') but will be looked up in nc_var_names_map
        vars_to_plot_info = {
            'msl': {'cmap': 'viridis', 'title': 'Mean Sea Level Pressure', 'cbar_label_base': 'Pressure', 'expected_units': 'hPa', 'panel_idx': 0, 'internal_key_suffix': 'msl'},
            'tcc': {'cmap': 'Blues', 'title': 'Total Cloud Cover', 'cbar_label_base': 'Cloud Cover', 'expected_units': '%', 'panel_idx': 1, 'internal_key_suffix': 'tcc'},
            't2m': {'cmap': 'RdYlBu_r', 'title': '2m Temperature', 'cbar_label_base': 'Temperature', 'expected_units': '°C', 'panel_idx': 2, 'internal_key_suffix': 't2m'},
            'wind': {'panel_idx': 3, 'internal_key_suffix_u': 'u10', 'internal_key_suffix_v': 'v10'} # Special case for wind
        }
        
        # Find the actual NetCDF names using the map provided from Streamlit's --vars and --varnames
        actual_nc_names = {}
        missing_vars = []
        for panel_key, info in vars_to_plot_info.items():
            if panel_key == 'wind':
                u_key_to_find = next((k for k in nc_var_names_map if k.endswith(info['internal_key_suffix_u'])), None)
                v_key_to_find = next((k for k in nc_var_names_map if k.endswith(info['internal_key_suffix_v'])), None)
                if u_key_to_find and v_key_to_find:
                    actual_nc_names['u10'] = nc_var_names_map[u_key_to_find]
                    actual_nc_names['v10'] = nc_var_names_map[v_key_to_find]
                else:
                    if not u_key_to_find: missing_vars.append(f"u10 (key suffix {info['internal_key_suffix_u']})")
                    if not v_key_to_find: missing_vars.append(f"v10 (key suffix {info['internal_key_suffix_v']})")
            else:
                key_to_find = next((k for k in nc_var_names_map if k.endswith(info['internal_key_suffix'])), None)
                if key_to_find:
                    actual_nc_names[panel_key] = nc_var_names_map[key_to_find]
                else:
                    missing_vars.append(f"{panel_key} (key suffix {info['internal_key_suffix']})")
        
        if missing_vars:
            print(f"Error: Could not find mappings for the following essential variables in the provided --vars/--varnames: {', '.join(missing_vars)}", file=sys.stderr)
            print(f"Provided mappings (internal_key -> nc_name): {nc_var_names_map}", file=sys.stderr)
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

        # Check if all required NetCDF variables are in the dataset
        for role, nc_name in actual_nc_names.items():
            if nc_name not in ds_ecmwf:
                print(f"Error: Expected NetCDF variable '{nc_name}' (for role '{role}') not found in {ecmwf_file_path}.", file=sys.stderr)
                sys.exit(1)

        # Determine spatial coordinate names
        lat_coord_name = 'latitude' if 'latitude' in ds_ecmwf.coords else 'lat'
        lon_coord_name = 'longitude' if 'longitude' in ds_ecmwf.coords else 'lon'
        if lat_coord_name not in ds_ecmwf.coords or lon_coord_name not in ds_ecmwf.coords:
            print(f"Error: Could not determine lat/lon coordinates. Found: {list(ds_ecmwf.coords)}", file=sys.stderr)
            sys.exit(1)
        longitude_coords = ds_ecmwf[lon_coord_name].values
        latitude_coords = ds_ecmwf[lat_coord_name].values

        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 12))
        axes_flat = axes.flatten()
        
        # Plot each variable (MSL, TCC, T2M)
        for var_key, info in vars_to_plot_info.items():
            if var_key == 'wind': continue # Handle wind separately

            ax = axes_flat[info['panel_idx']]
            nc_name = actual_nc_names[var_key]
            data_var = ds_ecmwf[nc_name]
            
            # Ensure time dimension is correctly identified for this variable
            dim_to_use_for_time = actual_time_coord_name
            if actual_time_coord_name not in data_var.dims:
                potential_time_dims = [d for d in data_var.dims if actual_time_coord_name in d or 'time' in d.lower()]
                if not potential_time_dims: print(f"Error: Time dim for '{nc_name}' not found.", file=sys.stderr); sys.exit(1)
                dim_to_use_for_time = potential_time_dims[0]

            var_avg = data_var.mean(dim=dim_to_use_for_time, skipna=True)
            display_units_val = info['expected_units']

            if var_key == 'msl':
                current_units_attr = data_var.attrs.get('units', 'Pa').lower()
                if var_avg.max(skipna=True).item() > 10000 and 'pa' in current_units_attr and 'hpa' not in current_units_attr:
                    var_avg = var_avg / 100
                    print(f"MSL data converted to hPa.", file=sys.stdout)
            elif var_key == 't2m':
                if var_avg.mean(skipna=True).item() > 100: # Heuristic for Kelvin
                    var_avg = var_avg - 273.15
                    print(f"T2M data converted to °C.", file=sys.stdout)
            elif var_key == 'tcc': # TCC is often 0-1, convert to %
                if var_avg.max(skipna=True).item() <= 1.0 and var_avg.min(skipna=True).item() >=0.0 :
                    var_avg = var_avg * 100
                    print(f"TCC data converted to %.", file=sys.stdout)
            
            var_plot = ax.pcolormesh(longitude_coords, latitude_coords, var_avg.data,
                                          cmap=info['cmap'], shading='auto')
            cbar = plt.colorbar(var_plot, ax=ax, orientation='vertical', pad=0.05, aspect=25)
            cbar.set_label(f'Average {info["cbar_label_base"]} ({display_units_val})', fontsize=10)
            ax.set_title(f'Average {info["title"]}', fontsize=12)
            ax.set_xlabel('Longitude', fontsize=10); ax.set_ylabel('Latitude', fontsize=10)
            ax.tick_params(axis='both', which='major', labelsize=8)

        # Plot average wind vectors in the last panel
        wind_info = vars_to_plot_info['wind']
        ax_wind = axes_flat[wind_info['panel_idx']]
        u10_nc_name = actual_nc_names['u10']
        v10_nc_name = actual_nc_names['v10']

        u10_data_var = ds_ecmwf[u10_nc_name]
        v10_data_var = ds_ecmwf[v10_nc_name]
        
        dim_u_time = actual_time_coord_name if actual_time_coord_name in u10_data_var.dims else [d for d in u10_data_var.dims if 'time' in d.lower()][0]
        dim_v_time = actual_time_coord_name if actual_time_coord_name in v10_data_var.dims else [d for d in v10_data_var.dims if 'time' in d.lower()][0]

        u10_avg = u10_data_var.mean(dim=dim_u_time, skipna=True)
        v10_avg = v10_data_var.mean(dim=dim_v_time, skipna=True)
        wind_speed_avg = np.sqrt(u10_avg**2 + v10_avg**2)
   
        speed_plot = ax_wind.pcolormesh(longitude_coords, latitude_coords, wind_speed_avg.data,
                                      cmap='GnBu', shading='auto') # Changed cmap for wind speed
   
        quiver_stride = max(1, len(longitude_coords) // 20, len(latitude_coords) // 20) # Dynamic stride
        
        # Ensure lon/lat for quiver match the potentially strided u10_avg/v10_avg if they are also strided
        # If u10_avg/v10_avg are full res, then lon_mesh/lat_mesh should also be full res before striding for quiver
        lon_mesh, lat_mesh = np.meshgrid(longitude_coords, latitude_coords)

        ax_wind.quiver(lon_mesh[::quiver_stride, ::quiver_stride], 
                       lat_mesh[::quiver_stride, ::quiver_stride],
                       u10_avg.data[::quiver_stride, ::quiver_stride],
                       v10_avg.data[::quiver_stride, ::quiver_stride],
                       scale=None, scale_units='inches', headwidth=4, headlength=5, width=0.004, color='black', alpha=0.7)
   
        cbar_wind = plt.colorbar(speed_plot, ax=ax_wind, orientation='vertical', pad=0.05, aspect=25)
        cbar_wind.set_label('Average Wind Speed (m/s)', fontsize=10)
        ax_wind.set_title('Average Wind Speed & Vectors', fontsize=12)
        ax_wind.set_xlabel('Longitude', fontsize=10); ax_wind.set_ylabel('Latitude', fontsize=10)
        ax_wind.tick_params(axis='both', which='major', labelsize=8)
   
        try:
            start_date = pd.to_datetime(str(ds_ecmwf[actual_time_coord_name].min().values)).strftime('%Y-%m-%d')
            end_date = pd.to_datetime(str(ds_ecmwf[actual_time_coord_name].max().values)).strftime('%Y-%m-%d')
            plt.suptitle(f'ECMWF Average Variables ({start_date} to {end_date})', fontsize=16)
        except Exception: plt.suptitle(f'ECMWF Average Variables', fontsize=16)
   
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        try:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"Multipanel plot saved successfully to {output_path}", file=sys.stdout)
        except Exception as e:
            print(f"Error: Failed to save multipanel plot to {output_path}: {e}", file=sys.stderr)
            sys.exit(1)
        finally:
            plt.close()

    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        # import traceback # For debugging
        # traceback.print_exc(file=sys.stderr) # For debugging
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a multipanel plot of time-averaged ECMWF variables.")
    parser.add_argument('--files', required=True, nargs=1, help='Path to the input ECMWF NetCDF file.')
    parser.add_argument('--varnames', required=True, nargs='+', help="NetCDF variable names for [msl, tcc, t2m, u10, v10] in order or as mapped by --vars.")
    parser.add_argument('--vars', required=True, nargs='+', help="Internal variable keys (e.g., ecmwf_msl, ecmwf_tcc, etc.) corresponding to varnames.")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type (unused by this script).')
    
    args = parser.parse_args()

    # We need to map the internal keys from --vars to the nc_names from --varnames
    # The Streamlit app sends these lists. `args.vars[i]` corresponds to `args.varnames[i]`.
    if len(args.vars) != len(args.varnames):
        print("Error: Number of --vars must match number of --varnames.", file=sys.stderr)
        sys.exit(1)
    
    # Check if all required internal key suffixes are present in args.vars
    expected_suffixes = ['msl', 'tcc', 't2m', 'u10', 'v10']
    provided_suffixes_in_keys = {key.split('_')[-1] for key in args.vars} # e.g. from 'ecmwf_msl' extracts 'msl'
    
    all_expected_found = True
    for suffix in expected_suffixes:
        if suffix not in provided_suffixes_in_keys:
            print(f"Error: Expected variable type '{suffix}' (based on internal key suffix) not found in --vars list: {args.vars}", file=sys.stderr)
            all_expected_found = False
    if not all_expected_found:
        sys.exit(1)

    # Create the nc_var_names_map for the plotting function
    # Maps internal_key (from --vars) to its nc_name (from --varnames)
    _nc_var_names_map = {args.vars[i]: args.varnames[i] for i in range(len(args.vars))}
        
    plot_ecmwf_multipanel_average(args.files[0], _nc_var_names_map, args.output)