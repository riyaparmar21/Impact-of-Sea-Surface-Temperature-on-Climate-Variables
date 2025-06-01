# #!/usr/bin/env python
# # coding: utf-8

# # In[4]:


# import xarray as xr
# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd
# from matplotlib.dates import DateFormatter
# import scipy.stats as stats


# # In[13]:


# sst_data = xr.open_dataset(r"noaa_sst_masked.nc")
# ecmwf_data = xr.open_dataset(r"ecmwf datas.nc")

# # Ensure same time period
# sst = sst_data.sst.sel(time=slice('2019-01-01', '2024-12-01'))
# tcc = ecmwf_data.tcc.sel(valid_time=slice('2019-01-01', '2024-12-01'))

# # Interpolate to common grid
# tcc_interp = tcc.interp(latitude=sst.lat, longitude=sst.lon)
# tcc_interp = tcc_interp.rename({'valid_time': 'time'})

# # Compute correlation at each grid point
# correlation = xr.corr(sst, tcc_interp, dim='time')

# # Plot
# plt.figure(figsize=(12, 10))
# im = plt.pcolormesh(sst.lon, sst.lat, correlation, cmap='RdBu_r', vmin=-1, vmax=1)
# plt.colorbar(im, label='Correlation Coefficient')
# plt.title('Spatial Correlation: SST vs Total Cloud Cover')

# # Add correlation values as text on the plot
# for i in range(len(sst.lat)):
#     for j in range(len(sst.lon)):
#         # Only show values at some grid points to avoid overcrowding
#         if i % 2 == 0 and j % 2 == 0:  # Show every other point
#             plt.text(sst.lon[j], sst.lat[i], 
#                      f'{correlation.values[i, j]:.2f}',  # Format to 2 decimal places
#                      ha='center', va='center', 
#                      color='black', fontsize=8,
#                      bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', pad=1))

# plt.xlabel('Longitude')
# plt.ylabel('Latitude')
# plt.tight_layout()
# plt.show()


# # In[ ]:



#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import sys
from pathlib import Path

def plot_sst_tcc_spatial_correlation(sst_file_path_str, sst_nc_var_name, 
                                     tcc_file_path_str, tcc_nc_var_name, 
                                     output_path_str,
                                     time_slice_start='2019-01-01', 
                                     time_slice_end='2024-12-01'):
    """
    Computes and plots the spatial correlation between SST and TCC.
    Overlays correlation values on the heatmap.
    """
    try:
        sst_file_path = Path(sst_file_path_str)
        tcc_file_path = Path(tcc_file_path_str)
        output_path = Path(output_path_str)

        if not sst_file_path.exists():
            print(f"Error: SST data file not found at {sst_file_path}", file=sys.stderr)
            sys.exit(1)
        if not tcc_file_path.exists():
            print(f"Error: TCC data file not found at {tcc_file_path}", file=sys.stderr)
            sys.exit(1)

        # Load SST data
        try:
            sst_ds = xr.open_dataset(sst_file_path, decode_times=True)
            if 'time' not in sst_ds.coords:
                print(f"Error: 'time' coordinate not found in SST dataset {sst_file_path}", file=sys.stderr)
                sys.exit(1)
        except Exception as e:
            print(f"Error loading SST dataset {sst_file_path}: {e}", file=sys.stderr)
            sys.exit(1)
        
        if sst_nc_var_name not in sst_ds:
            print(f"Error: SST variable '{sst_nc_var_name}' not found in {sst_file_path}", file=sys.stderr)
            sys.exit(1)
        sst_data_var = sst_ds[sst_nc_var_name]

        # Load TCC data
        tcc_time_coord = None
        try:
            tcc_ds = xr.open_dataset(tcc_file_path, decode_times=True)
            if 'valid_time' in tcc_ds.coords:
                tcc_time_coord = 'valid_time'
            elif 'time' in tcc_ds.coords:
                tcc_time_coord = 'time'
            else:
                print(f"Error: 'time' or 'valid_time' coord not in TCC dataset {tcc_file_path}", file=sys.stderr)
                sys.exit(1)
            if tcc_time_coord != 'time':
                tcc_ds = tcc_ds.rename({tcc_time_coord: 'time'})
            tcc_time_coord = 'time'
        except Exception as e:
            print(f"Error loading TCC dataset {tcc_file_path}: {e}", file=sys.stderr)
            sys.exit(1)

        if tcc_nc_var_name not in tcc_ds:
            print(f"Error: TCC variable '{tcc_nc_var_name}' not found in {tcc_file_path}", file=sys.stderr)
            sys.exit(1)
        tcc_data_var = tcc_ds[tcc_nc_var_name]

        # Ensure common time period
        try:
            sst_sliced = sst_data_var.sel(time=slice(time_slice_start, time_slice_end))
            tcc_sliced = tcc_data_var.sel(time=slice(time_slice_start, time_slice_end))
        except Exception as e_slice:
            print(f"Error slicing time for SST or TCC: {e_slice}", file=sys.stderr)
            print(f"SST time range: {sst_data_var.time.min().values} to {sst_data_var.time.max().values}", file=sys.stderr)
            print(f"TCC time range: {tcc_data_var.time.min().values} to {tcc_data_var.time.max().values}", file=sys.stderr)
            sys.exit(1)
        
        if sst_sliced.time.size < 2 or tcc_sliced.time.size < 2:
            print(f"Error: Insufficient time points (SST: {sst_sliced.time.size}, TCC: {tcc_sliced.time.size}) for correlation.", file=sys.stderr)
            sys.exit(1)

        # Interpolate TCC to SST grid
        sst_lat_name = 'lat' if 'lat' in sst_sliced.coords else 'latitude'
        sst_lon_name = 'lon' if 'lon' in sst_sliced.coords else 'longitude'
        if sst_lat_name not in sst_sliced.coords or sst_lon_name not in sst_sliced.coords:
            print(f"Error: Could not find lat/lon for SST grid. Coords: {list(sst_sliced.coords)}", file=sys.stderr)
            sys.exit(1)

        tcc_lat_name = 'latitude' if 'latitude' in tcc_sliced.coords else 'lat'
        tcc_lon_name = 'longitude' if 'longitude' in tcc_sliced.coords else 'lon'
        if tcc_lat_name not in tcc_sliced.coords or tcc_lon_name not in tcc_sliced.coords:
            print(f"Error: Could not find lat/lon for TCC grid. Coords: {list(tcc_sliced.coords)}", file=sys.stderr)
            sys.exit(1)

        try:
            tcc_interp = tcc_sliced.interp({tcc_lat_name: sst_sliced[sst_lat_name], 
                                            tcc_lon_name: sst_sliced[sst_lon_name]})
        except Exception as e_interp:
            print(f"Error during interpolation of TCC to SST grid: {e_interp}", file=sys.stderr)
            sys.exit(1)

        # Align time steps
        try:
            sst_aligned, tcc_interp_aligned = xr.align(sst_sliced, tcc_interp, join='inner')
        except Exception as e_align:
            print(f"Error aligning SST and interpolated TCC time steps: {e_align}", file=sys.stderr)
            sys.exit(1)

        if sst_aligned.time.size < 2:
            print(f"Error: Less than 2 common time points after alignment ({sst_aligned.time.size}). Cannot compute correlation.", file=sys.stderr)
            sys.exit(1)

        # Check valid data points per grid point
        valid_counts = (~sst_aligned.isnull() & ~tcc_interp_aligned.isnull()).sum(dim='time')
        if (valid_counts < 2).any():
            print(f"Warning: Some grid points have < 2 valid time points for correlation. These will be masked.", file=sys.stdout)

        # Compute correlation, masking invalid points
        correlation = xr.corr(sst_aligned, tcc_interp_aligned, dim='time')
        correlation = correlation.where(valid_counts >= 2)  # Mask points with < 2 valid data points

        # Plot
        plt.figure(figsize=(12, 10))
        im = plt.pcolormesh(correlation[sst_lon_name], correlation[sst_lat_name], correlation.data, 
                            cmap='RdBu_r', vmin=-1, vmax=1, shading='auto')
        plt.colorbar(im, label='Correlation Coefficient', fraction=0.046, pad=0.04)
        plt.title('Spatial Correlation: SST vs Total Cloud Cover', fontsize=14)

        # Add correlation values as text
        lat_stride = max(1, len(correlation[sst_lat_name]) // 10)
        lon_stride = max(1, len(correlation[sst_lon_name]) // 10)

        for i in range(0, len(correlation[sst_lat_name]), lat_stride):
            for j in range(0, len(correlation[sst_lon_name]), lon_stride):
                if not np.isnan(correlation.data[i, j]):
                    plt.text(correlation[sst_lon_name][j].item(), correlation[sst_lat_name][i].item(), 
                             f'{correlation.data[i, j]:.2f}',
                             ha='center', va='center', 
                             color='black', fontsize=7,
                             bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=0.5))

        plt.xlabel('Longitude', fontsize=12)
        plt.ylabel('Latitude', fontsize=12)
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
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute and plot spatial correlation between SST and TCC.")
    parser.add_argument('--files', required=True, nargs='+', 
                        help='Paths to the input NetCDF files (SST and TCC sources).')
    parser.add_argument('--varnames', required=True, nargs='+',
                        help="NetCDF variable names. Order should match --vars.")
    parser.add_argument('--vars', required=True, nargs=2,
                        help="Two internal variable keys (e.g., 'noaa_sst' 'ecmwf_tcc').")
    parser.add_argument('--output', required=True, help='Path to save the output PNG image.')
    parser.add_argument('--plot_type', help='Name of the plot type (unused by this script).')
    parser.add_argument('--time_start', default='2019-01-01', help='Start date for time slicing (YYYY-MM-DD).')
    parser.add_argument('--time_end', default='2024-12-01', help='End date for time slicing (YYYY-MM-DD).')

    args = parser.parse_args()

    if len(args.varnames) != len(args.vars) or len(args.vars) != 2:
        print("Error: This script expects exactly two --vars (internal keys) and a corresponding number of --varnames.", file=sys.stderr)
        sys.exit(1)
    
    SST_KEY_IDENTIFIER_SUBSTRING = 'sst'
    TCC_KEY_IDENTIFIER_SUBSTRING = 'tcc'

    SST_FILE_BASENAME_EXPECTED = "noaa_sst_masked.nc"
    TCC_FILE_BASENAME_EXPECTED = "ecmwf datas.nc"

    sst_file_arg, tcc_file_arg = None, None
    sst_varname_arg, tcc_varname_arg = None, None
    
    key1, key2 = args.vars[0], args.vars[1]
    varname1, varname2 = args.varnames[0], args.varnames[1]

    if SST_KEY_IDENTIFIER_SUBSTRING in key1.lower() and TCC_KEY_IDENTIFIER_SUBSTRING in key2.lower():
        sst_varname_arg, tcc_varname_arg = varname1, varname2
        sst_file_arg = next((f for f in args.files if Path(f).name == SST_FILE_BASENAME_EXPECTED), None)
        tcc_file_arg = next((f for f in args.files if Path(f).name == TCC_FILE_BASENAME_EXPECTED), None)
    elif TCC_KEY_IDENTIFIER_SUBSTRING in key1.lower() and SST_KEY_IDENTIFIER_SUBSTRING in key2.lower():
        tcc_varname_arg, sst_varname_arg = varname1, varname2
        tcc_file_arg = next((f for f in args.files if Path(f).name == TCC_FILE_BASENAME_EXPECTED), None)
        sst_file_arg = next((f for f in args.files if Path(f).name == SST_FILE_BASENAME_EXPECTED), None)
    else:
        print(f"Error: Could not identify SST and TCC roles from --vars: {args.vars}.", file=sys.stderr)
        sys.exit(1)

    if not sst_file_arg or not tcc_file_arg:
        print("Error: Could not determine file paths for SST and TCC.", file=sys.stderr)
        sys.exit(1)

    print(f"Plotting SST-TCC Spatial Correlation: SST from '{sst_file_arg}' (var: '{sst_varname_arg}') and TCC from '{tcc_file_arg}' (var: '{tcc_varname_arg}')", file=sys.stdout)

    plot_sst_tcc_spatial_correlation(sst_file_arg, sst_varname_arg, 
                                     tcc_file_arg, tcc_varname_arg, 
                                     args.output,
                                     args.time_start, args.time_end)
