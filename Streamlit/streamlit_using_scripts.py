#!/usr/bin/env python
# coding: utf-8

import streamlit as st
import xarray as xr
import os
import subprocess # For running external scripts
import tempfile   # For temporary directories
from pathlib import Path # Ensure Path is imported
from PIL import Image
import pandas as pd
import sys # For sys.executable in run_plotting_script

# --- Page Configuration ---
st.set_page_config(page_title="SST Impact Analysis-demo", layout="wide")

# --- File Definitions ---
# Use pathlib.Path for all path manipulations
SCRIPT_DIRECTORY = Path(os.path.dirname(os.path.abspath(__file__))) # Convert to Path object

NOAA_PRECIP_FILE = SCRIPT_DIRECTORY / 'precip_konkan.nc'
NOAA_SST_FILE = SCRIPT_DIRECTORY / 'noaa_sst_masked.nc'
ECMWF_FILE =  SCRIPT_DIRECTORY / 'ecmwf datas.nc'
PLOTTING_SCRIPTS_DIR = SCRIPT_DIRECTORY / 'plotting_scripts'
CYCLONE_IMAGES_DIR = SCRIPT_DIRECTORY / 'cyclone images'

if not PLOTTING_SCRIPTS_DIR.exists():
    st.error(f"CRITICAL: Plotting scripts directory not found at '{PLOTTING_SCRIPTS_DIR}'. Application cannot generate plots.")
    st.stop()

# --- Plot Script Mapping ---
PLOT_SCRIPT_MAP = {
    "SST-PRECIP Combined Analysis": PLOTTING_SCRIPTS_DIR / "sst-precip combined analyis.py",
    "SST-2M TEMP Combined Analysis": PLOTTING_SCRIPTS_DIR / "sst - t2m combined analysis.py",
    "SST-MSL Combined Analysis": PLOTTING_SCRIPTS_DIR / "sst - msl combined analysis.py",
    "SST-TCC Combined Analysis": PLOTTING_SCRIPTS_DIR / "sst - tcc combined analysis.py",
    "SST-U10 Combined Analysis": PLOTTING_SCRIPTS_DIR / "sst - u10 combined analysis.py",
    "SST-V10 Combined Analysis": PLOTTING_SCRIPTS_DIR / "sst - v10 combined analysis.py",
    "SST-TCC Spatial Correlation": PLOTTING_SCRIPTS_DIR / "sst - tcc spatial correlation with values.py",
    "Correlation Matrix Plot": PLOTTING_SCRIPTS_DIR / "correlation matrix plot.py",
    "ECMWF Multipanel Overview": PLOTTING_SCRIPTS_DIR / "multipanel ecmwf variables avg heatmaps.py",
    "PRECIP Average Heatmap": PLOTTING_SCRIPTS_DIR / "precip avg heatmap.py",
    "PRECIP seasonal time series": PLOTTING_SCRIPTS_DIR / "precipitation seasonal time series.py",
    "PRECIP Seasonal Heatmaps": PLOTTING_SCRIPTS_DIR / "precip seasonal avg heatmaps.py",
    "PRECIP Timeseries": PLOTTING_SCRIPTS_DIR / "precipitation area avg time series.py",
    "PRECIP At Selected Grid Points": PLOTTING_SCRIPTS_DIR / "precipitation at selected grid points.py",
    "SST Average Heatmap": PLOTTING_SCRIPTS_DIR / "SST avg heatmap.py",
    "SST Anomalies Timeseries": PLOTTING_SCRIPTS_DIR / "SST anomalies time series.py",
    "SST Area Average Timeseries": PLOTTING_SCRIPTS_DIR / "SST area avg time series.py",
    "SST Seasonal Heatmaps": PLOTTING_SCRIPTS_DIR / "SST seasonal avg heatmaps.py",
    "SST Seasonal Cycle Timeseries": PLOTTING_SCRIPTS_DIR / "SST seasonal cycle time series.py",
    "T2M Area Average Timeseries": PLOTTING_SCRIPTS_DIR / "2m temp area avg time series.py",
    "T2M Monthly Mean Barplot": PLOTTING_SCRIPTS_DIR / "2m temp monthly mean bar plot.py",
    "T2M Average Heatmap": PLOTTING_SCRIPTS_DIR / "2m temp avg heatmap.py",
    "T2M Seasonal Heatmaps": PLOTTING_SCRIPTS_DIR / "2m temp seasonal avg heatmaps.py",
    "T2M Seasonal Timeseries": PLOTTING_SCRIPTS_DIR / "2m temp seasonal time series.py",
    "MSL Average Heatmap": PLOTTING_SCRIPTS_DIR / "msl avg heatmap.py",
    "MSL Area Average Timeseries": PLOTTING_SCRIPTS_DIR / "msl area avg time series.py",
    "MSL Seasonal Heatmaps": PLOTTING_SCRIPTS_DIR / "msl seasonal avg heatmaps.py",
    "MSL Seasonal Timeseries": PLOTTING_SCRIPTS_DIR / "msl seasonal time series.py",
    "TCC Average Heatmap": PLOTTING_SCRIPTS_DIR / "total cloud cover avg heatmap.py",
    "TCC Area Average Timeseries": PLOTTING_SCRIPTS_DIR / "total cloud cover area avg time series.py",
    "TCC Seasonal Heatmaps": PLOTTING_SCRIPTS_DIR / "total cloud cover seasonal avg heatmaps.py",
    "TCC Seasonal Timeseries": PLOTTING_SCRIPTS_DIR / "total cloud cover seasonal time series.py",
    "TCC Seasonal Mean Barplot": PLOTTING_SCRIPTS_DIR / "seasonal mean cloud cover bar plot.py",
    "U10 Average Heatmap": PLOTTING_SCRIPTS_DIR / "u-10 wind avg heatmap.py",
    "U10 Area Average Timeseries": PLOTTING_SCRIPTS_DIR / "u wind area avg time series.py",
    "U10 Seasonal Heatmaps": PLOTTING_SCRIPTS_DIR / "u-10 wind seasonal avg heatmaps.py",
    "U10 Seasonal Timeseries": PLOTTING_SCRIPTS_DIR / "u wind seasonal time series.py",
    "V10 Average Heatmap": PLOTTING_SCRIPTS_DIR / "v-10 wind avg heatmap.py",
    "V10 Area Average Timeseries": PLOTTING_SCRIPTS_DIR / "v wind area avg time series.py",
    "V10 Seasonal Heatmaps": PLOTTING_SCRIPTS_DIR / "v-10 wind seasonal avg heatmaps.py",
    "V10 Seasonal Timeseries": PLOTTING_SCRIPTS_DIR / "v wind seasonal time series.py",
    "Wind Speed Area Average": PLOTTING_SCRIPTS_DIR / "10m wind speed area avg.py",
    "Wind Vectors Average Heatmap": PLOTTING_SCRIPTS_DIR / "wind vectors avg heatmap.py",
    "Generic Two Variable Timeseries": PLOTTING_SCRIPTS_DIR / "plot_combined_timeseries.py",
}
PLOT_SCRIPT_MAP = {k.strip(): v for k, v in PLOT_SCRIPT_MAP.items()}

# --- Variable Details Mapping ---
VARIABLE_DETAILS_MAP = {
    'precip': {'file_key': 'noaa_precip_data', 'nc_var_name': 'precip', 'display_name': 'Precipitation', 'data_file': NOAA_PRECIP_FILE, 'time_coord': 'time'},
    'sst':    {'file_key': 'noaa_sst_data',    'nc_var_name': 'sst',    'display_name': 'Sea Surface Temperature', 'data_file': NOAA_SST_FILE,    'time_coord': 'time'},
    'msl':    {'file_key': 'ecmwf_msl_data',   'nc_var_name': 'msl',    'display_name': 'Mean Sea Level Pressure', 'data_file': ECMWF_FILE,     'time_coord': 'valid_time'},
    'tcc':    {'file_key': 'ecmwf_tcc_data',   'nc_var_name': 'tcc',    'display_name': 'Total Cloud Cover',       'data_file': ECMWF_FILE,     'time_coord': 'valid_time'},
    't2m':    {'file_key': 'ecmwf_t2m_data',   'nc_var_name': 't2m',    'display_name': '2 meter Temperature',   'data_file': ECMWF_FILE,     'time_coord': 'valid_time'},
    'v10':    {'file_key': 'ecmwf_v10_data',   'nc_var_name': 'v10',    'display_name': 'V Wind (10m)',          'data_file': ECMWF_FILE,     'time_coord': 'valid_time'},
    'u10':    {'file_key': 'ecmwf_u10_data',   'nc_var_name': 'u10',    'display_name': 'U Wind (10m)',          'data_file': ECMWF_FILE,     'time_coord': 'valid_time'},
}
USER_FRIENDLY_VARS = sorted(VARIABLE_DETAILS_MAP.keys())
VARS_FOR_SST_COMBINED = sorted([v for v in USER_FRIENDLY_VARS if v != 'sst'])

# --- Master Config from successfully loaded files ---
master_config = {}
for uf_key, details in VARIABLE_DETAILS_MAP.items():
    if details['data_file'].exists(): # This now correctly uses Path.exists()
        master_config[details['file_key']] = {
            'data_file': details['data_file'], # This is a Path object
            'nc_var_name': details['nc_var_name'],
            'user_friendly_name': uf_key,
            'display_name': details['display_name'],
            'time_coord': details['time_coord']
        }
    else:
        st.sidebar.warning(f"Data file for '{details['display_name']}' not found: {details['data_file']}. This variable will be unavailable.")

if not master_config:
    st.error("CRITICAL: No usable data files found. Application cannot proceed with plotting.")
    # st.stop() # Consider stopping if no data is available

USER_FRIENDLY_VARS_FILTERED = sorted(list(set(details['user_friendly_name'] for fk, details in master_config.items())))
VARS_FOR_SST_COMBINED_FILTERED = sorted([v for v in USER_FRIENDLY_VARS_FILTERED if v != 'sst'])

# --- Analysis Plot Mapping ---
ANALYSIS_PLOT_MAPPING = {
    "Single Variable Analysis": {
        'precip': ["PRECIP Average Heatmap", "PRECIP Seasonal Heatmaps", "PRECIP Timeseries", "PRECIP At Selected Grid Points","PRECIP seasonal time series"],
        'sst': ["SST Average Heatmap", "SST Anomalies Timeseries", "SST Area Average Timeseries", "SST Seasonal Heatmaps", "SST Seasonal Cycle Timeseries"],
        't2m': ["T2M Area Average Timeseries", "T2M Monthly Mean Barplot", "T2M Average Heatmap", "T2M Seasonal Heatmaps", "T2M Seasonal Timeseries"],
        'msl': ["MSL Average Heatmap", "MSL Area Average Timeseries", "MSL Seasonal Heatmaps", "MSL Seasonal Timeseries"],
        'tcc': ["TCC Average Heatmap", "TCC Area Average Timeseries", "TCC Seasonal Heatmaps", "TCC Seasonal Timeseries", "TCC Seasonal Mean Barplot"],
        'u10': ["U10 Average Heatmap", "U10 Area Average Timeseries", "U10 Seasonal Heatmaps", "U10 Seasonal Timeseries"],
        'v10': ["V10 Average Heatmap", "V10 Area Average Timeseries", "V10 Seasonal Heatmaps", "V10 Seasonal Timeseries"],
    },
    "SST Combined Analysis": {
        'precip': ["SST-PRECIP Combined Analysis"],
        't2m': ["SST-2M TEMP Combined Analysis"],
        'msl': ["SST-MSL Combined Analysis"],
        'tcc': ["SST-TCC Combined Analysis", "SST-TCC Spatial Correlation"],
        'u10': ["SST-U10 Combined Analysis"],
        'v10': ["SST-V10 Combined Analysis"],
    },
    "Two Variable Combined Timeseries": ["Generic Two Variable Timeseries"],
    "Correlation Matrix (All Variables)": ["Correlation Matrix Plot"],
    "ECMWF All Variables Overview": ["ECMWF Multipanel Overview", "Wind Speed Area Average", "Wind Vectors Average Heatmap"],
}

# --- Plot Explanations ---
def create_two_var_key(var1_ufk, var2_ufk):
    if not var1_ufk or not var2_ufk: return None
    return "_vs_".join(sorted((str(var1_ufk), str(var2_ufk))))

DEFAULT_PLOT_EXPLANATION = "Detailed explanation for this plot is pending."
PLOT_EXPLANATIONS = { 
    # ... (Your extensive PLOT_EXPLANATIONS dictionary remains unchanged) ...
    create_two_var_key('precip', 'v10'): """**Precipitation vs. 10 m Wind Speed**\n\nOver 2019–2024 the monsoon (Jun–Sep) rainfall peaks coincide with elevated 10 m wind speeds, illustrating the seasonal coupling of moisture and momentum. A sharp joint spike in June 2020 reflects Cyclone Nisarga’s heavy rains and storm‐force winds. Smaller simultaneous upticks in May 2021 and June 2023 mark Tauktae and Biparjoy passages. Outside storm periods, wind and rain retain distinct but overlapping seasonal phasing.""",
    # (Include all your explanations here)
    "Wind Vectors Average Heatmap": """Overlays average wind vectors (combining u10 and v10) atop mean SST contours. The arrows illustrate prevailing wind directions and magnitudes, revealing coastal upwelling zones where winds run parallel to the shore. Linking this with SST patterns helps interpret how wind-driven mixing affects regional sea surface temperatures."""
}


# --- Cyclone Data ---
CYCLONES_DATA = [
    {"name": "Cyclone Kyarr (2019)", "region": "Western India, Oman, Yemen, Somalia", "image": CYCLONE_IMAGES_DIR / "kyaar.jpeg", "explanation": "Cyclone Kyarr (October 2019) skirted..."},
    {"name": "Cyclone Nisarga (2020)", "region": "Maharashtra (Alibag),Maharashtra,", "image": CYCLONE_IMAGES_DIR / "nisarga.jpg", "explanation": "Cyclone Nisarga (June 2020) made landfall..."},
    {"name": "Cyclone Tauktae (2021)", "region": "Gujarat, Maharashtra, Karnataka, Goa, Kerala", "image": CYCLONE_IMAGES_DIR / "tauktae.png", "explanation": "Cyclone Tauktae (May 2021) tracked parallel..."},
    {"name": "Cyclone Biparjoy (2023)", "region": "Gujarat (Did not make landfall as expected, skirted coast)", "image": CYCLONE_IMAGES_DIR / "Biparjoy.jpg", "explanation": "Cyclone Biparjoy (June 2023) passed north..."}
]
DATASETS_TO_DOWNLOAD = {
    "NOAA Precipitation Data": NOAA_PRECIP_FILE, 
    "NOAA SST Data": NOAA_SST_FILE, 
    "ECMWF Combined Data (All Vars)": ECMWF_FILE,
}

# --- Helper Function to Run Plotting Scripts ---
def run_plotting_script(script_path: Path, vars_info_for_script: dict, output_image_path: Path, plot_type_name: str):
    if not script_path.exists():
        st.error(f"Plotting script not found: {script_path}")
        return False, None, f"Plotting script not found: {script_path}", None
    
    # Check for Generic Two Variable Timeseries specifically, as it has unique argument needs.
    is_generic_ts_plot = (plot_type_name == "Generic Two Variable Timeseries")

    if not vars_info_for_script and is_generic_ts_plot:
        # This specific plot type *always* needs vars_info. Other scripts might be self-contained.
        st.error(f"Script '{script_path.name}' for plot type '{plot_type_name}' requires variable information, but none was provided.")
        return False, None, "Missing variable arguments for generic script.", None
    
    details_list_for_cmd = list(vars_info_for_script.values()) if vars_info_for_script else []
    script_internal_keys_cmd = list(vars_info_for_script.keys()) if vars_info_for_script else []

    unique_files_for_cmd = list(dict.fromkeys([str(d['data_file']) for d in details_list_for_cmd]))
    script_nc_varnames_for_cmd = [d['nc_var_name'] for d in details_list_for_cmd]
    script_display_names_for_cmd = [d['display_name'] for d in details_list_for_cmd]
    script_time_coords_for_cmd = [d['time_coord'] for d in details_list_for_cmd]
    
    # Convert Path objects to strings for subprocess command
    cmd = [ sys.executable, str(script_path), "--output", str(output_image_path), "--plot_type", str(plot_type_name)]
    
    if script_internal_keys_cmd: cmd.extend(["--vars"] + script_internal_keys_cmd)
    if unique_files_for_cmd: cmd.extend(["--files"] + unique_files_for_cmd)
    if script_nc_varnames_for_cmd: cmd.extend(["--varnames"] + script_nc_varnames_for_cmd)
    
    if is_generic_ts_plot:
        if script_display_names_for_cmd: cmd.extend(["--displaynames"] + script_display_names_for_cmd)
        if script_time_coords_for_cmd: cmd.extend(["--timecoords"] + script_time_coords_for_cmd)
        # Argument count validation for Generic Two Variable Timeseries
        if not (len(unique_files_for_cmd) in [1,2] and \
                len(script_nc_varnames_for_cmd) == 2 and \
                len(script_display_names_for_cmd) == 2 and \
                len(script_time_coords_for_cmd) == 2):
            err_msg = f"Argument mismatch for '{plot_type_name}'. Expected 2 sets of var, displayname, timecoord, and 1 or 2 files. Got: {len(script_nc_varnames_for_cmd)} vars, {len(unique_files_for_cmd)} files."
            st.error(err_msg); return False, None, err_msg, None
            
    st.info(f"Running: `{' '.join(cmd)}`")
    try:
        process = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=300)
        stdout, stderr = process.stdout, process.stderr
        if process.returncode != 0: return False, stdout, stderr, None
        
        generated_csv_path = None
        if plot_type_name == "Correlation Matrix Plot":
            # output_image_path is a Path object, so with_name works.
            potential_csv_path = output_image_path.with_name(output_image_path.stem + "_data.csv")
            if potential_csv_path.exists(): generated_csv_path = potential_csv_path
            elif process.returncode == 0 : st.caption(f"Note: Corr Matrix script ran, but CSV not found: {potential_csv_path}")
        return True, stdout, stderr, generated_csv_path
    except FileNotFoundError: 
        err_msg = f"Error: Python interpreter ({sys.executable}) or script path ({script_path}) incorrect or not found."; 
        st.error(err_msg); return False, None, err_msg, None
    except subprocess.TimeoutExpired: 
        err_msg = f"Error: Script {script_path.name} timed out after 300 seconds."; 
        st.error(err_msg); return False, "Script timed out.", err_msg, None
    except Exception as e: 
        err_msg = f"Unexpected error running {script_path.name}: {e}"; 
        st.error(err_msg); return False, str(e), err_msg, None

# --- Streamlit App UI and Logic ---
st.title("SST Impact Analysis on climate variables")
st.markdown("""**Region:** Arabian Sea near Konkan Coast, Maharashtra (Approx. 16°N - 21°N, 69°E - 74°E)""")
st.markdown("---")

with st.sidebar:
    st.header("Tools & Info")
    if st.button("🌀 Cyclone Information (2019-2023)", key="cyclone_info_button_sidebar"):
        st.session_state.show_cyclone_dialog = not st.session_state.get('show_cyclone_dialog', False)
    
    with st.expander("📥 Download Datasets", expanded=False):
        for dn, fp in DATASETS_TO_DOWNLOAD.items():
            if fp.exists(): # fp is a Path object
                with open(fp, "rb") as file_content: # open() works with Path objects
                    st.download_button(f"Download {dn} ({fp.name})", file_content, fp.name, "application/x-netcdf", key=f"dl_{dn.replace(' ', '_')}")
            else: 
                st.caption(f"{dn} ({fp.name}) not found at {fp}.") # Provide full path for debugging
    
    st.markdown("---"); st.header("Plotting Controls")

    if 'view_all_plots_mode' not in st.session_state: st.session_state.view_all_plots_mode = False
    cb_view_all = st.checkbox("🔓 View All Available Plot Scripts", value=st.session_state.view_all_plots_mode, key="view_all_cb")
    if cb_view_all != st.session_state.view_all_plots_mode:
        st.session_state.view_all_plots_mode = cb_view_all; st.rerun()

    analysis_type_selected_sb = None
    selected_primary_var_ufk_sb = None
    selected_secondary_var_ufk_for_sst_sb = None
    selected_var1_ufk_for_combined_ts_sb = None
    selected_var2_ufk_for_combined_ts_sb = None

    if st.session_state.view_all_plots_mode:
        st.subheader("All Available Plot Scripts")
        analysis_type_selected_sb = "View All Plots Mode"
        # PLOT_SCRIPT_MAP values are Path objects, so .exists() works
        available_plot_scripts_sb = sorted([name for name, path_obj in PLOT_SCRIPT_MAP.items() if path_obj.exists()])
        if len(available_plot_scripts_sb) != len(PLOT_SCRIPT_MAP): 
            st.caption("Note: Some plot scripts in config not found on disk.")
        unavailable_plot_scripts_info_sb = [] # In this mode, we only list what's available
    else: # Guided Mode
        analysis_type_selected_sb = st.radio( "Select Analysis Type:",
            ["Single Variable Analysis", "SST Combined Analysis", "Two Variable Combined Timeseries", "Correlation Matrix (All Variables)", "ECMWF All Variables Overview"],
            key="analysis_type_radio_guided", horizontal=True )
        
        if analysis_type_selected_sb == "Single Variable Analysis":
            if not USER_FRIENDLY_VARS_FILTERED: st.warning("No variables available for single analysis.")
            else: selected_primary_var_ufk_sb = st.selectbox("Select Variable:", USER_FRIENDLY_VARS_FILTERED, format_func=lambda x: VARIABLE_DETAILS_MAP.get(x, {}).get('display_name', x.upper()), key="single_var_sel_sb")
        
        elif analysis_type_selected_sb == "SST Combined Analysis":
            if 'sst' not in USER_FRIENDLY_VARS_FILTERED: st.warning("SST data is not available for combined analysis.")
            elif not VARS_FOR_SST_COMBINED_FILTERED: st.warning("No other variables available to combine with SST.")
            else: selected_secondary_var_ufk_for_sst_sb = st.selectbox("Combine SST with:", VARS_FOR_SST_COMBINED_FILTERED, format_func=lambda x: VARIABLE_DETAILS_MAP.get(x, {}).get('display_name', x.upper()), key="sst_comb_var_sel_sb")
        
        elif analysis_type_selected_sb == "Two Variable Combined Timeseries":
            st.caption("Select two different variables for combined timeseries plot.")
            if len(USER_FRIENDLY_VARS_FILTERED) < 2: st.warning("At least two variables must be available for this plot type.")
            else:
                c1, ca, c2 = st.columns([5,1,5])
                with c1: selected_var1_ufk_for_combined_ts_sb = st.selectbox("Variable 1:", USER_FRIENDLY_VARS_FILTERED, format_func=lambda x: VARIABLE_DETAILS_MAP.get(x, {}).get('display_name', x.upper()), key="comb_ts_v1_sb")
                with ca: st.markdown("<br><p style='text-align:center;font-size:1.1em'>&</p>", unsafe_allow_html=True)
                with c2:
                    opts_v2 = [v for v in USER_FRIENDLY_VARS_FILTERED if v != selected_var1_ufk_for_combined_ts_sb] if selected_var1_ufk_for_combined_ts_sb else []
                    if selected_var1_ufk_for_combined_ts_sb and not opts_v2: st.warning("No other variable available for Variable 2.")
                    elif opts_v2: selected_var2_ufk_for_combined_ts_sb = st.selectbox("Variable 2:", opts_v2, format_func=lambda x: VARIABLE_DETAILS_MAP.get(x, {}).get('display_name', x.upper()), key="comb_ts_v2_sb")
        
        plot_candidates_guided = []
        current_mapping = ANALYSIS_PLOT_MAPPING.get(analysis_type_selected_sb)
        if current_mapping:
            if analysis_type_selected_sb == "Single Variable Analysis" and selected_primary_var_ufk_sb: plot_candidates_guided = current_mapping.get(selected_primary_var_ufk_sb, [])
            elif analysis_type_selected_sb == "SST Combined Analysis" and selected_secondary_var_ufk_for_sst_sb: plot_candidates_guided = current_mapping.get(selected_secondary_var_ufk_for_sst_sb, [])
            elif analysis_type_selected_sb == "Two Variable Combined Timeseries":
                if selected_var1_ufk_for_combined_ts_sb and selected_var2_ufk_for_combined_ts_sb and selected_var1_ufk_for_combined_ts_sb != selected_var2_ufk_for_combined_ts_sb:
                    plot_candidates_guided = ["Generic Two Variable Timeseries"] # This is the specific plot type name
            elif isinstance(current_mapping, list): plot_candidates_guided = current_mapping
        
        # Filter candidates by actual script existence
        available_plot_scripts_sb = [p_name for p_name in plot_candidates_guided if p_name in PLOT_SCRIPT_MAP and PLOT_SCRIPT_MAP[p_name].exists()]
        unavailable_plot_scripts_info_sb = [p_name for p_name in plot_candidates_guided if p_name not in available_plot_scripts_sb]


    st.subheader("Select Plot Script(s) to Run")
    selected_plot_script_names_final = []
    if analysis_type_selected_sb == "Two Variable Combined Timeseries":
        if selected_var1_ufk_for_combined_ts_sb and selected_var2_ufk_for_combined_ts_sb and selected_var1_ufk_for_combined_ts_sb != selected_var2_ufk_for_combined_ts_sb:
            if "Generic Two Variable Timeseries" in PLOT_SCRIPT_MAP and PLOT_SCRIPT_MAP["Generic Two Variable Timeseries"].exists():
                selected_plot_script_names_final = ["Generic Two Variable Timeseries"]
            else: st.warning("Script for 'Generic Two Variable Timeseries' missing or not configured.")
    elif available_plot_scripts_sb:
        selected_plot_script_names_final = st.multiselect("Choose from available plot scripts:", available_plot_scripts_sb, default=[], key="plot_script_ms_common")
    elif st.session_state.view_all_plots_mode: st.warning("No plot scripts found or configured in 'View All Plots Mode'.")
    elif analysis_type_selected_sb: st.info("No plot scripts available for the current guided selection. Check if data files are loaded and scripts exist.")

    if unavailable_plot_scripts_info_sb and not st.session_state.view_all_plots_mode and analysis_type_selected_sb != "Two Variable Combined Timeseries":
        with st.expander("Unavailable/missing scripts for current selection", expanded=False):
            for plot_name_unavail in sorted(list(set(unavailable_plot_scripts_info_sb))): st.caption(f"- {plot_name_unavail} (Script file missing or not in PLOT_SCRIPT_MAP)")

# --- Cyclone Dialog ---
if st.session_state.get('show_cyclone_dialog', False):
    # Call st.dialog directly. Subsequent st commands will populate it.
    st.dialog(title="Major Cyclones Affecting the Region (2019-2023)")
    
    # Content of the dialog starts here
    st.header("Notable Cyclones") 
    if CYCLONES_DATA:
        cyclone_tab_names = [cyc['name'] for cyc in CYCLONES_DATA]
        if cyclone_tab_names:
            cyclone_display_tabs = st.tabs(cyclone_tab_names) 
            for i, tab_item in enumerate(cyclone_display_tabs):
                with tab_item: # st.tabs returns context managers for each tab
                    cyc_data_item = CYCLONES_DATA[i]
                    st.write(f"**Region of Impact:** {cyc_data_item['region']}")
                    if cyc_data_item['image'].exists():
                        try: 
                            st.image(Image.open(cyc_data_item['image']), caption=cyc_data_item['name'], use_container_width=True) 
                        except Exception as e_img: st.warning(f"Could not load image for {cyc_data_item['name']}: {e_img}") 
                    else: st.warning(f"Image file not found: {cyc_data_item['image']}")
                    st.markdown("---"); st.subheader("Details & Impact")
                    st.markdown(cyc_data_item.get("explanation", "Detailed explanation pending."))
    else: st.write("No cyclone data configured.") 
    
    if st.button("Close Cyclone Info", key="close_cyclone_dialog_button_standard"): 
        st.session_state.show_cyclone_dialog = False; st.rerun()

# --- Main Area: Plot Display ---
st.header("Data Visualization")

if not selected_plot_script_names_final:
    if st.session_state.view_all_plots_mode: st.info("⬅️ In 'View All Plots Mode'. Select plot scripts from the sidebar to generate visualizations.")
    else: st.info("⬅️ Complete selections in the sidebar to generate visualizations.")
else:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir) # tmp_path is a Path object
        
        def determine_vars_for_script_call(plot_script_name, analysis_type, 
                                           primary_var_ufk, sst_secondary_ufk, 
                                           ts_var1_ufk, ts_var2_ufk, 
                                           is_view_all, master_cfg, var_details_map, user_vars_filt):
            vars_to_pass_dict = {} 
            
            # Ensure all keys exist before trying to access them
            def _add_var_to_dict(ufk_key):
                if ufk_key in var_details_map and var_details_map[ufk_key]['file_key'] in master_cfg:
                    vars_to_pass_dict[var_details_map[ufk_key]['file_key']] = master_cfg[var_details_map[ufk_key]['file_key']]
                else:
                    st.caption(f"Warning: Details or master_config entry missing for user friendly key: {ufk_key}")

            if is_view_all:
                if plot_script_name == "Correlation Matrix Plot":
                    for ufk_ in user_vars_filt: _add_var_to_dict(ufk_)
                elif plot_script_name == "ECMWF Multipanel Overview":
                    for ufk_ in user_vars_filt: 
                        if ufk_ in var_details_map and var_details_map[ufk_]['data_file'] == ECMWF_FILE: _add_var_to_dict(ufk_)
                elif plot_script_name == "Generic Two Variable Timeseries":
                    if len(user_vars_filt) >= 2:
                        _add_var_to_dict(user_vars_filt[0])
                        _add_var_to_dict(user_vars_filt[1])
                    else: st.warning("Generic Timeseries needs at least 2 available variables in View All mode."); return {}
                else: 
                    found = False
                    for at_map_key, mapping_val in ANALYSIS_PLOT_MAPPING.items():
                        if isinstance(mapping_val, dict):
                            for p_var_key, p_list_val in mapping_val.items():
                                if plot_script_name in p_list_val:
                                    if at_map_key == "Single Variable Analysis": _add_var_to_dict(p_var_key)
                                    elif at_map_key == "SST Combined Analysis":
                                        _add_var_to_dict('sst')
                                        _add_var_to_dict(p_var_key)
                                    found = True; break
                        if found: break
                    if not found and not vars_to_pass_dict: # If specific inference fails and dict is still empty
                        st.caption(f"For '{plot_script_name}' in View All Mode: Could not auto-determine variables. The script might be self-contained or require manual argument setup if it fails.")
                        # For truly general scripts or those that load all data, this might be okay.
                        # If it's expected to take specific vars, this could lead to failure.

            elif analysis_type == "Single Variable Analysis" and primary_var_ufk: _add_var_to_dict(primary_var_ufk)
            elif analysis_type == "SST Combined Analysis" and sst_secondary_ufk:
                _add_var_to_dict('sst')
                _add_var_to_dict(sst_secondary_ufk)
            elif analysis_type == "Two Variable Combined Timeseries" and ts_var1_ufk and ts_var2_ufk:
                _add_var_to_dict(ts_var1_ufk)
                _add_var_to_dict(ts_var2_ufk)
            elif analysis_type == "Correlation Matrix (All Variables)":
                for ufk_ in user_vars_filt: _add_var_to_dict(ufk_)
            elif analysis_type == "ECMWF All Variables Overview":
                for ufk_ in user_vars_filt:
                    if ufk_ in var_details_map and var_details_map[ufk_]['data_file'] == ECMWF_FILE: _add_var_to_dict(ufk_)
            
            if not vars_to_pass_dict and plot_script_name != "Generic Two Variable Timeseries" and plot_script_name not in ["Correlation Matrix Plot", "ECMWF Multipanel Overview"]:
                 # Only warn if vars_to_pass_dict is empty AND it's not one of the known general plots.
                 # Or if it's the generic timeseries, which explicitly needs vars.
                if not (is_view_all and plot_script_name not in ["Generic Two Variable Timeseries"]): # Allow some view_all to proceed
                    st.warning(f"No variables determined to pass to script '{plot_script_name}'. It might fail if it expects specific inputs.")


            return vars_to_pass_dict

        # Logic for single plot or multiple plots in tabs
        plot_execution_failed = False # Flag to track if any plot fails

        if len(selected_plot_script_names_final) == 1:
            plot_script_name_to_run_single = selected_plot_script_names_final[0]
            title_context_single = ""
            # ... (determine title_context_single - same as your existing logic) ...
            if st.session_state.view_all_plots_mode: title_context_single = "(from All Plots list)"
            elif analysis_type_selected_sb == "Single Variable Analysis" and selected_primary_var_ufk_sb: title_context_single = f"for {VARIABLE_DETAILS_MAP.get(selected_primary_var_ufk_sb, {}).get('display_name', selected_primary_var_ufk_sb.upper())}"
            # ... (add other conditions for title_context_single)

            st.subheader(f"Displaying: {plot_script_name_to_run_single} {title_context_single}")
            script_to_run_single = PLOT_SCRIPT_MAP.get(plot_script_name_to_run_single) # Path object

            if not script_to_run_single or not script_to_run_single.exists(): 
                st.error(f"Script for '{plot_script_name_to_run_single}' not found/configured at {script_to_run_single}.")
                plot_execution_failed = True
            else:
                vars_for_script_call_single = determine_vars_for_script_call(
                    plot_script_name_to_run_single, analysis_type_selected_sb,
                    selected_primary_var_ufk_sb, selected_secondary_var_ufk_for_sst_sb,
                    selected_var1_ufk_for_combined_ts_sb, selected_var2_ufk_for_combined_ts_sb,
                    st.session_state.view_all_plots_mode, master_config, VARIABLE_DETAILS_MAP, USER_FRIENDLY_VARS_FILTERED
                )
                
                # Specific check for generic timeseries needing vars
                if plot_script_name_to_run_single == "Generic Two Variable Timeseries" and not vars_for_script_call_single:
                    st.error(f"Critical: Variables for '{plot_script_name_to_run_single}' could not be determined. Cannot generate plot.")
                    plot_execution_failed = True
                else:
                    output_img_name_single_str = f"{plot_script_name_to_run_single.replace(' ', '_').lower()}_plot.png"
                    if plot_script_name_to_run_single == "Generic Two Variable Timeseries":
                        v1k, v2k = None, None
                        if analysis_type_selected_sb == "Two Variable Combined Timeseries" and selected_var1_ufk_for_combined_ts_sb and selected_var2_ufk_for_combined_ts_sb:
                            v1k, v2k = selected_var1_ufk_for_combined_ts_sb, selected_var2_ufk_for_combined_ts_sb
                        elif st.session_state.view_all_plots_mode and len(USER_FRIENDLY_VARS_FILTERED) >= 2:
                            v1k, v2k = USER_FRIENDLY_VARS_FILTERED[0], USER_FRIENDLY_VARS_FILTERED[1]
                        if v1k and v2k: output_img_name_single_str = f"timeseries_{v1k}_vs_{v2k}.png"
                    
                    output_img_path_single = tmp_path / output_img_name_single_str # Path object
                    
                    success_s, stdout_s, stderr_s, csv_path_s = run_plotting_script(script_to_run_single, vars_for_script_call_single, output_img_path_single, plot_script_name_to_run_single)
                    
                    if stdout_s: 
                        with st.expander(f"Output from script: {script_to_run_single.name}", expanded=False): st.code(stdout_s, language='text')
                    if stderr_s: 
                        with st.expander(f"Errors/Warnings from script: {script_to_run_single.name}", expanded=not success_s): st.code(stderr_s, language='bash')

                    if success_s and output_img_path_single.exists():
                        st.image(str(output_img_path_single), caption=f"Generated: {plot_script_name_to_run_single}", use_container_width=True) # str() for st.image
                        with open(output_img_path_single, "rb") as f_img_s: st.download_button(f"Download Plot", f_img_s, output_img_name_single_str, "image/png", key=f"dl_s_{plot_script_name_to_run_single.replace(' ','_')}")
                        if csv_path_s and csv_path_s.exists(): # csv_path_s is Path object
                            with open(csv_path_s, "rb") as f_csv_s: st.download_button(f"Download CSV Data", f_csv_s, csv_path_s.name, "text/csv", key=f"dl_csv_s_{plot_script_name_to_run_single.replace(' ','_')}")
                        
                        st.markdown("---"); st.subheader("Plot Explanation")
                        exp_key_s = plot_script_name_to_run_single
                        if plot_script_name_to_run_single == "Generic Two Variable Timeseries":
                            # ... (logic to determine exp_key_s for generic timeseries, same as before) ...
                            v1fk, v2fk = None,None
                            if analysis_type_selected_sb == "Two Variable Combined Timeseries": v1fk, v2fk = selected_var1_ufk_for_combined_ts_sb, selected_var2_ufk_for_combined_ts_sb
                            elif st.session_state.view_all_plots_mode and len(USER_FRIENDLY_VARS_FILTERED) >=2 : v1fk, v2fk = USER_FRIENDLY_VARS_FILTERED[0], USER_FRIENDLY_VARS_FILTERED[1]
                            if v1fk and v2fk: exp_key_s = create_two_var_key(v1fk, v2fk)
                        st.markdown(PLOT_EXPLANATIONS.get(exp_key_s, DEFAULT_PLOT_EXPLANATION))
                    elif success_s: 
                        st.error(f"Plotting script '{script_to_run_single.name}' ran successfully but the output image was not found at '{output_img_path_single}'.")
                        plot_execution_failed = True
                    elif not success_s:
                        st.error(f"Failed to generate plot for '{plot_script_name_to_run_single}'. Check error messages above.")
                        plot_execution_failed = True


        else: # Multiple plots in tabs
            tabs_objects = st.tabs([f"📊 {name}" for name in selected_plot_script_names_final])
            for i, tab_item_disp in enumerate(tabs_objects):
                with tab_item_disp:
                    current_plot_script_name_tab = selected_plot_script_names_final[i]
                    title_context_tab = ""
                     # ... (determine title_context_tab - same as your existing logic) ...
                    if st.session_state.view_all_plots_mode: title_context_tab = "(from All Plots list)"
                    elif analysis_type_selected_sb: title_context_tab = f"({analysis_type_selected_sb})" # Simplified for tabs

                    st.markdown(f"#### {current_plot_script_name_tab} {title_context_tab}")
                    script_to_run_tab = PLOT_SCRIPT_MAP.get(current_plot_script_name_tab) # Path object

                    if not script_to_run_tab or not script_to_run_tab.exists(): 
                        st.error(f"Script for '{current_plot_script_name_tab}' not found/configured at {script_to_run_tab}.")
                        plot_execution_failed = True; continue # Skip to next tab
                    
                    vars_for_script_call_tab = determine_vars_for_script_call(
                        current_plot_script_name_tab, analysis_type_selected_sb,
                        selected_primary_var_ufk_sb, selected_secondary_var_ufk_for_sst_sb,
                        selected_var1_ufk_for_combined_ts_sb, selected_var2_ufk_for_combined_ts_sb,
                        st.session_state.view_all_plots_mode, master_config, VARIABLE_DETAILS_MAP, USER_FRIENDLY_VARS_FILTERED
                    )

                    if current_plot_script_name_tab == "Generic Two Variable Timeseries" and not vars_for_script_call_tab:
                        st.error(f"Critical: Variables for '{current_plot_script_name_tab}' could not be determined. Cannot generate plot in this tab.")
                        plot_execution_failed = True; continue
                    
                    output_img_name_tab_str = f"{current_plot_script_name_tab.replace(' ', '_').lower()}_plot.png"
                    if current_plot_script_name_tab == "Generic Two Variable Timeseries":
                        v1k_t, v2k_t = None, None
                        if analysis_type_selected_sb == "Two Variable Combined Timeseries" and selected_var1_ufk_for_combined_ts_sb and selected_var2_ufk_for_combined_ts_sb:
                            v1k_t, v2k_t = selected_var1_ufk_for_combined_ts_sb, selected_var2_ufk_for_combined_ts_sb
                        elif st.session_state.view_all_plots_mode and len(USER_FRIENDLY_VARS_FILTERED) >= 2:
                             v1k_t, v2k_t = USER_FRIENDLY_VARS_FILTERED[0], USER_FRIENDLY_VARS_FILTERED[1]
                        if v1k_t and v2k_t: output_img_name_tab_str = f"timeseries_{v1k_t}_vs_{v2k_t}.png"

                    output_img_path_tab_tmp = tmp_path / output_img_name_tab_str # Path object
                    
                    success_t, stdout_t, stderr_t, csv_path_t = run_plotting_script(script_to_run_tab, vars_for_script_call_tab, output_img_path_tab_tmp, current_plot_script_name_tab)

                    if stdout_t:
                        with st.expander(f"Output from script: {script_to_run_tab.name}", expanded=False): st.code(stdout_t, language='text')
                    if stderr_t:
                        with st.expander(f"Errors/Warnings from script: {script_to_run_tab.name}", expanded=not success_t): st.code(stderr_t, language='bash')

                    if success_t and output_img_path_tab_tmp.exists():
                        st.image(str(output_img_path_tab_tmp), caption=f"Generated: {current_plot_script_name_tab}", use_container_width=True) # str() for st.image
                        with open(output_img_path_tab_tmp, "rb") as f_img_t: st.download_button(f"Download Plot", f_img_t, output_img_name_tab_str, "image/png", key=f"dl_t_{current_plot_script_name_tab.replace(' ','_')}_{i}")
                        if csv_path_t and csv_path_t.exists(): # csv_path_t is Path object
                            with open(csv_path_t, "rb") as f_csv_t: st.download_button(f"Download CSV Data", f_csv_t, csv_path_t.name, "text/csv", key=f"dl_csv_t_{current_plot_script_name_tab.replace(' ','_')}_{i}")
                        
                        st.markdown("---"); st.subheader("Plot Explanation")
                        exp_key_t = current_plot_script_name_tab
                        if current_plot_script_name_tab == "Generic Two Variable Timeseries":
                            # ... (logic to determine exp_key_t for generic timeseries, same as before) ...
                            v1fk_t, v2fk_t = None,None
                            if analysis_type_selected_sb == "Two Variable Combined Timeseries": v1fk_t, v2fk_t = selected_var1_ufk_for_combined_ts_sb, selected_var2_ufk_for_combined_ts_sb
                            elif st.session_state.view_all_plots_mode and len(USER_FRIENDLY_VARS_FILTERED) >=2 : v1fk_t, v2fk_t = USER_FRIENDLY_VARS_FILTERED[0], USER_FRIENDLY_VARS_FILTERED[1]
                            if v1fk_t and v2fk_t: exp_key_t = create_two_var_key(v1fk_t, v2fk_t)
                        st.markdown(PLOT_EXPLANATIONS.get(exp_key_t, DEFAULT_PLOT_EXPLANATION))
                    elif success_t: 
                        st.error(f"Plotting script '{script_to_run_tab.name}' ran successfully but the output image was not found at '{output_img_path_tab_tmp}'.")
                        plot_execution_failed = True
                    elif not success_t:
                        st.error(f"Failed to generate plot for '{current_plot_script_name_tab}' in this tab. Check error messages above.")
                        plot_execution_failed = True
        
        if plot_execution_failed:
            st.sidebar.error("One or more plots failed to generate. Check messages in the main area.")


st.markdown("---")
st.markdown("""
<style>.footer-left {text-align: left; font-size: 16px; margin-top: 50px;} .footer-left img {height: 16px; vertical-align: middle; margin-left: 5px;}</style>
<div class="footer-left">Developed by: Riya Parmar <a href="https://www.linkedin.com/in/riya-parmar-44890a26b/" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/c/ca/LinkedIn_logo_initials.png" alt="LinkedIn"></a></div>
""", unsafe_allow_html=True)
