#!/usr/bin/env python
# coding: utf-8

from Streamlit.plot_helpers import my_plot_function

import streamlit as st
import xarray as xr
import os
import subprocess # For running external scripts
import tempfile   # For temporary directories
from pathlib import Path
from PIL import Image
import pandas as pd
import sys # For sys.executable in run_plotting_script

# --- Page Configuration ---
st.set_page_config(page_title="SST Impact Analysis-demo", layout="wide")

# --- File Definitions ---
BASE_DIR = Path(".")
NOAA_PRECIP_FILE = BASE_DIR / "precip_konkan.nc"
ECMWF_FILE = BASE_DIR / "ecmwf datas.nc"
NOAA_SST_FILE = BASE_DIR / "noaa_sst_masked.nc"
CYCLONE_IMAGES_DIR = BASE_DIR / "cyclone images"
PLOTTING_SCRIPTS_DIR = BASE_DIR / "plotting_scripts"

if not PLOTTING_SCRIPTS_DIR.exists():
    st.error(f"CRITICAL: Plotting scripts directory not found at '{PLOTTING_SCRIPTS_DIR}'. Application cannot generate plots.")

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
    if details['data_file'].exists():
        master_config[details['file_key']] = {
            'data_file': details['data_file'],
            'nc_var_name': details['nc_var_name'],
            'user_friendly_name': uf_key,
            'display_name': details['display_name'],
            'time_coord': details['time_coord']
        }
    else:
        st.sidebar.warning(f"Data file for '{details['display_name']}' not found: {details['data_file']}. This variable will be unavailable.")

if not master_config:
    st.error("CRITICAL: No usable data files found. Application cannot proceed with plotting.")
    # st.stop()

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
# Config Check (Optional but good for debugging)
# ... (config check logic can be kept here) ...

# --- Plot Explanations ---
def create_two_var_key(var1_ufk, var2_ufk): # Takes user_friendly_keys
    if not var1_ufk or not var2_ufk: return None
    return "_vs_".join(sorted((str(var1_ufk), str(var2_ufk))))

DEFAULT_PLOT_EXPLANATION = "Detailed explanation for this plot is pending."
PLOT_EXPLANATIONS = { 
    create_two_var_key('precip', 'v10'): """**Precipitation vs. 10 m Wind Speed**\n\nOver 2019–2024 the monsoon (Jun–Sep) rainfall peaks coincide with elevated 10 m wind speeds, illustrating the seasonal coupling of moisture and momentum. A sharp joint spike in June 2020 reflects Cyclone Nisarga’s heavy rains and storm‐force winds. Smaller simultaneous upticks in May 2021 and June 2023 mark Tauktae and Biparjoy passages. Outside storm periods, wind and rain retain distinct but overlapping seasonal phasing.""",

    create_two_var_key('sst', 'v10'): """**SST vs. 10 m Wind Speed**\n\nSea‐surface temperature rises steadily into the pre-monsoon maximum (Mar–May), then cools under monsoon cloud cover, while surface winds strengthen with monsoon onset. Cyclone Kyarr in Oct 2019 produces a brief SST dip and wind surge, evidencing storm‐driven mixing. Later events (Tauktae, Biparjoy) induce similar but muted wind peaks and minor SST cooling. Overall, wind bursts during cyclones punctuate the smooth seasonal SST cycle.""",

    create_two_var_key('t2m', 'u10'): """**2 m Temperature vs. U-Wind Component**\n\nNear-surface air temperature peaks in April–May each year then falls through the monsoon, whereas the zonal (u10) wind component intensifies with the southwest flow. Cyclone Tauktae (May 2021) shows a transient u10 boost and slight t2m dip as the storm’s westerly winds push ashore. Minor temperature fluctuations accompany other cyclones but recoup quickly post‐event. The underlying pattern remains a clear pre-monsoon thermal maximum and monsoon wind strengthening.""",

    create_two_var_key('msl', 'precip'): """**MSL Pressure vs. Precipitation**\n\nMean sea level pressure falls sharply at monsoon onset, inversely tracking rainfall peaks. During Nisarga (Jun 2020) and Biparjoy (Jun 2023), MSLP plunges coincide with torrential precipitation spikes, marking classic cyclone signatures. Outside those storms, monsoon trough dynamics yield broad low‐pressure periods with sustained rain. The strong anticorrelation underscores how pressure drops drive convective rainfall in this region.""",

    create_two_var_key('msl', 'sst'): """**MSL Pressure vs. SST**\n\nSeasonal MSLP minima in peak summer align with SST maxima just beforehand, as land/sea temperature gradients deepen. Kyarr (Oct 2019) briefly forces an anomalous pressure drop with a concomitant slight SST dip, reflecting enhanced mixing. Similar but smaller MSLP dips and SST coolings appear around Tauktae and Biparjoy. Otherwise, both variables follow their own seasonal rhythms—thermal buildup then monsoon‐driven cooling.""",

    create_two_var_key('msl', 'v10'): """**MSL Pressure vs. 10 m Wind Speed**\n\nLow‐pressure monsoon months bring stronger surface winds; this inverse relationship holds across 2019–2024. Notably, Tauktae and Biparjoy yield sharp MSLP troughs with concurrent wind‐speed peaks, illustrating storm intensity. Between cyclones, wind variations trace gradual pressure changes as the monsoon advances and retreats. The tight coupling during extreme events contrasts with smoother seasonal cycles otherwise.""",

    create_two_var_key('sst', 'tcc'): """**SST vs. Total Cloud Cover**\n\nSST reaches annual highs in pre-monsoon (Mar–May), while total cloud cover surges during the monsoon (Jun–Sep). Cyclone Nisarga stands out as a massive TCC spike in Jun 2020, alongside a brief SST drop from upwelling and mixing. Biparjoy in Jun 2023 shows a smaller cloud burst and surface cooling. Beyond storms, the eye-pleasing seesaw between clear, warm pre-monsoon skies and overcast monsoon persists.""",

    create_two_var_key('precip', 'tcc'): """**Precipitation vs. Total Cloud Cover**\n\nMonsoon months exhibit tightly coupled peaks of rainfall and cloud cover, with cloudiness arcing just ahead of rain maxima. Storm events like Nisarga and Biparjoy appear as abrupt concurrent surges—extreme TCC and intense precipitation within days. Outside cyclones, the smooth monsoon build-up and retreat maintain a strong precipitation–cloud relationship. The regional convective regime is clearly manifest.""",

    create_two_var_key('msl', 't2m'): """**MSL Pressure vs. 2 m Temperature**\n\nMean pressure lows during the monsoon coincide with moderate dips in 2 m temperature, driven by cloud shading and evaporative cooling. Cyclones Kyarr and Nisarga cause pronounced pressure drops with accompanying t2m dips of 1–2 °C. These thermal anomalies quickly rebound post-event. Overall, seasonal temperature minima align with monsoon low‐pressure periods.""",

    create_two_var_key('tcc', 'u10'): """**Cloud Cover vs. U-Wind Component**\n\nDuring Jun–Sep the u10 component strengthens alongside peak cloud cover, as westerly monsoon flow feeds convection. Cyclone Tauktae produces a clear transient jump in both u10 and TCC in May 2021, highlighting storm‐enhanced cloudiness. Smaller bumps accompany Nisarga and Biparjoy. Otherwise, cloudiness and zonal winds rise and fall together with the monsoon’s advance and withdrawal.""",

    create_two_var_key('sst', 't2m'): """**SST vs. 2 m Temperature**\n\nSST and air temperature both rise from winter into pre-monsoon, reflecting strong seasonal insolation. During June–September, SST levels off or drops slightly while 2 m temperature decreases more steeply under cloud cover. Cyclones like Nisarga and Biparjoy show momentary t2m cooling, while SST lags in its response. The pairing highlights surface–air heat exchange patterns disturbed only briefly by storms.""",

    create_two_var_key('tcc', 'v10'): """**Cloud Cover vs. V-Wind Component**\n\nV10 winds (southerly/northerly) intensify during the monsoon months in tandem with higher cloud fractions. Cyclones like Tauktae and Biparjoy cause sharp TCC spikes and increased v10 variability due to turbulent vertical transport and shear. Outside cyclone events, the seasonal march of monsoon cloud decks tracks well with elevated v10 readings, signifying vertical moisture movement and convective depth.""",

    create_two_var_key('sst', 'u10'): """**SST vs. U-Wind Component**\n\nSurface westerlies (u10) rise during the monsoon as SST begins to decline under cloud cover and evaporative cooling. Cyclone Kyarr’s late‐season burst in Oct 2019 generates a wind surge with subtle SST response. More prominent wind peaks emerge during Tauktae and Biparjoy in May–June. Overall, SST responds slowly to wind‐driven mixing, but storm pulses interrupt this gradual seasonal cooling arc.""",

    create_two_var_key('t2m', 'v10'): """**2 m Temperature vs. V-Wind Component**\n\nFrom pre-monsoon into mid-summer, rising v10 winds (south-north flow) accompany falling 2 m temperatures, as overcast skies block solar input. Cyclone Nisarga’s June 2020 signature includes a temporary v10 spike and t2m cooling. Similar patterns occur in May 2021 (Tauktae) and June 2023 (Biparjoy). The broader pattern reveals cooling of the land surface under strong monsoonal vertical winds.""",

    create_two_var_key('msl', 'u10'): """**MSL Pressure vs. U-Wind Component**\n\nPressure drops through the monsoon season coincide with strengthening westerly winds (u10), a hallmark of the southwest monsoon. Tauktae’s May 2021 landfall brings a steep drop in MSL and a matching u10 spike. Biparjoy echoes this pattern in 2023. These outliers punctuate the otherwise smooth seasonal descent in pressure and rise in zonal wind strength.""",

    create_two_var_key('tcc', 't2m'): """**Cloud Cover vs. 2 m Temperature**\n\nHigh total cloud cover during the monsoon limits surface heating, suppressing 2 m temperatures. Each year this inverse relationship is clear, especially in Jul–Aug. Cyclones further deepen this contrast: Tauktae and Biparjoy cause dense cloud spikes with brief temperature dips. The data captures how insolation control by clouds dominates regional near-surface temperatures.""",

    create_two_var_key('precip', 'sst'): """**Precipitation vs. SST**\n\nRising SST pre-monsoon precedes rain onset by 1–2 months, then dips during peak monsoon due to oceanic cooling under cloud and rain. Cyclone Nisarga produces intense rain and slight SST cooling in Jun 2020; Tauktae and Biparjoy also show SST dips during peak rainfall. The curve emphasizes how upper-ocean heat content modulates rainfall, and how cyclones briefly tip this balance.""",

    create_two_var_key('precip', 't2m'): """**Precipitation vs. 2 m Temperature**\n\nPrecipitation ramps up as 2 m temperatures fall with monsoon onset in June. Cyclone Biparjoy in Jun 2023 produces a dramatic rain spike and accompanying t2m dip. Nisarga and Tauktae show similar but smaller signatures. Outside these events, the pre-monsoon to monsoon shift displays a consistent pattern of declining temperature as rainfall surges.""",

    create_two_var_key('t2m', 'sst'): """**2 m Temperature vs. SST**\n\nBoth SST and 2 m temperature rise steadily through March–May, then diverge: t2m drops faster with monsoon cloud cover, while SST lags. During cyclones (e.g., Kyarr, Nisarga), t2m responds with short-lived cooling, while SST exhibits a delayed, milder dip. This pairing highlights land–sea thermal gradients during summer and storm-driven disruptions.""",

    create_two_var_key('u10', 'v10'): """**U vs. V Wind Components**\n\nU10 and V10 wind components surge in synchrony during the monsoon, reflecting southwest flow dominance. Cyclones create temporary divergence—e.g., Tauktae (May 2021) and Biparjoy (Jun 2023) show offset wind spikes from directional shifts. Outside storm windows, the seasonal co‐variation suggests broad monsoon circulation patterns dominate both wind directions.""",
    "SST-PRECIP Combined Analysis": """
Here SST vs. precipitation correlations are displayed spatially. Moderate positive values (0.3–0.6) along the coast indicate that warmer seas boost evaporation and moisture convergence, driving heavier rainfall in the Konkan belt.
""",

    "SST-2M TEMP Combined Analysis": """
This focused correlation heatmap isolates the relationship between 2 m air temperature and SST at each grid point. High coefficients (= 0.9) along the coastal cells confirm that near-surface air heating is tightly coupled to ocean warming, especially during summer months. 
""",

    "SST-MSL Combined Analysis": """
A scatter-density plot with a trend line quantifying SST–MSL correlation across all time steps and grid cells. The negative slope visually confirms the inverse relationship: each 1 °C increase in SST corresponds to a notable drop in sea level pressure. This scatter view complements the spatial map by emphasizing the overall strength and linearity of the linkage.

""",

    "SST-TCC Combined Analysis": """
Annotates SST vs. total cloud cover correlations with exact values (≈ 0.4–0.6 offshore). Warmer SSTs drive enhanced convection and cloud formation, a key feedback for radiation and precipitation processes that you capture.
""",

    "SST-U10 Combined Analysis": """
This figure depicts the correlation between SST and the zonal wind (u10) component. Positive correlations (0.4–0.7) near the shore suggest that warm SSTs bolster easterly winds by altering pressure gradients.
""",

    "SST-V10 Combined Analysis": """
A map of SST vs. meridional wind (v10) correlations, showing values up to \~0.9 in some coastal cells. Strong positive coupling indicates that SST-driven thermal contrasts significantly influence north–south wind patterns, shaping monsoon incursions inland.
""",

    "SST-TCC Spatial Correlation": """
This annotated heatmap shows SST vs. total cloud cover (TCC) correlation with exact numerical labels in each cell. Values around 0.4–0.6 in offshore regions highlight moderate positive coupling—warmer water surfaces lead to more cloud formation.
""",

    "Correlation Matrix Plot": """
A full correlation matrix heatmap showing Pearson correlation coefficients among NOAA precipitation, NOAA SST, ECMWF MSL, TCC, u10, v10, and 2 m temperature. Warm colors highlight strong positive relationships (e.g., TCC vs. precipitation), while cool colors mark negative links (e.g., MSL vs. SST). 
""",

    "ECMWF Multipanel Overview": """
A multi-panel figure displaying the climatological mean of all key ECMWF variables (MSL, TCC, u10, v10, T2m) alongside NOAA SST. Each subplot uses a unified color scale to facilitate cross-variable comparisons. This panel synthesizes baseline conditions, showing—for instance—how high SST zones coincide with low-pressure troughs and enhanced cloud cover.
""",

    "PRECIP Average Heatmap": """
Shows the 2019–2024 mean precipitation at each grid cell. Heaviest rainfall concentrates along the central Konkan coast, tapering inland. This baseline clarifies where SST-induced moisture convergence most strongly translates into precipitation, informing spatial weighting in your model inputs.

""",

    "PRECIP Seasonal Heatmaps": """
Four-panel maps of seasonal precipitation ('Dec-Jan-Feb', 'Mar-Apr-May', 'Jun-Jul-Aug', 'Sep-Oct-Nov'). Jun-Jul-Aug exhibits a dramatic coastal rainfall spike with values often exceeding 400 mm, reflecting monsoon peak. In contrast, Dec-Jan-Feb is nearly dry. These seasonal snapshots help you learn distinct precipitation regimes tied to SST cycles.
""",

    "PRECIP Timeseries": """
Monthly precipitation averaged over the Konkan grid, with clear annual peaks during monsoon months and lows in winter. This temporal profile establishes the strong seasonality your LSTM must model.
""",

    "PRECIP At Selected Grid Points": """
Overlays precipitation time series from specific coastal and inland grid points. Contrasting profiles reveal that coastal cells have sharper monsoon peaks and larger interannual variability than inland cells.
""",

    "PRECIP Distribution Over Time": """
A violin (or box-whisker) plot showing the distribution of monthly precipitation values across years. It visualizes variability and outliers, indicating that July and August have high medians.
""",

    "PRECIP Anomalies Timeseries": """
Time series of precipitation anomalies (deviation from monthly climatological mean). Positive spikes denote wetter-than-normal months, while negatives mark drought phases. 
""",

    "PRECIP 6-Point Moving Average": """
Monthly precipitation with a 3-month moving average overlay. The smoothing emphasizes multi-month trends, revealing anomalous monsoon years (e.g., 2021 dip or 2023 surge). Incorporating moving averages can help you capture longer-term persistence in rainfall anomalies.
""",

    "SST Average Heatmap": """
Displays the mean sea surface temperature (SST) from 2019 to 2024 at each grid cell. The warmest waters (˜29–30 °C) concentrate along the coastal strip, tapering offshore. This baseline highlights the persistent heat reservoirs that modulate regional convection and drive monsoon dynamics in your climate-impact analysis.
""",

    "SST Anomalies Timeseries": """
Time series of SST anomalies—deviations from the multi-year monthly mean. Positive anomalies in 2020 and 2023 indicate warmer-than-normal years. Such anomaly tracking is crucial for identifying the impact of unusual ocean warming on rainfall and pressure patterns.
""",

    "SST Area Average Timeseries": """
Monthly area‐averaged SST for the Konkan domain, smoothing spatial variability to reveal broad seasonal trends and interannual differences. Key for correlating SST patterns with regional climate anomalies.
""",

    "SST Seasonal Heatmaps": """
Four-panel seasonal SST maps ('Dec-Jan-Feb', 'Mar-Apr-May', 'Jun-Jul-Aug', 'Sep-Oct-Nov'). Pre-monsoon (Mar-Apr-May) and peak monsoon (Jun-Jul-Aug) panels show maximum SSTs, with coastal temperatures peaking near 31 °C. Winter (Dec-Jan-Feb) depicts cooling to around 26 °C. These snapshots emphasize how seasonal SST swings set the stage for atmospheric responses.
""",

    "SST Seasonal Cycle Timeseries": """
A line plot of the monthly climatological SST cycle. It reveals a gradual rise from January’s \~27 °C to a June–July peak near 30 °C, followed by a decline into December. This smooth cycle aids by pinpointing lags between SST peaks and downstream atmospheric effects.
""",

    "T2M Area Average Timeseries": """
Area-averaged time series of 2 m temperature for the Konkan grid. It shows a clear seasonal cycle with maxima in May–June and minima in December.
""",

    "T2M Monthly Mean Barplot": """
Bar chart of monthly mean 2 m temperatures. The bars quantify the amplitude of seasonal warming, peaking pre-monsoon. This visualization aids in communicating baseline thermal loads when SST anomalies occur.
""",

    "T2M Seasonal Heatmaps": """
Four-panel seasonal maps of 2 m temperature. Summer (Sep-Oct-Nov, Jun-Jul-Aug) panels show the strongest coastal warming, while winter (Dec-Jan-Feb) has more uniform cool temperatures. These spatial patterns guide the CNN component to learn season-specific temperature gradients.
""",

    "T2M Seasonal Timeseries": """
Line plot of the seasonal cycle (monthly climatology) of 2 m temperature. Smooth curves highlight inter-month transitions, emphasizing how air temperature lags or leads SST—critical information.
""",

    "MSL Average Heatmap": """
Heatmap of mean sea level pressure (MSL) averaged over 2019–2024. It highlights low-pressure belts offshore where warmer SSTs drive enhanced convection. Coastal pressure minima align with known monsoon trough locations, critical for modeling storm genesis in your LSTM–CNN framework.
""",

    "MSL Area Average Timeseries": """
Spatial map of MSL at a representative time (e.g., July 2022). It visualizes pressure gradients across the Konkan coast, revealing zones where low pressure intensifies. Comparing this snapshot with concurrent SST maps can pinpoint specific events when ocean warmth triggered pressure drops.
""",

    "MSL Seasonal Heatmaps": """
Seasonal MSL maps ('Dec-Jan-Feb', 'Mar-Apr-May', 'Jun-Jul-Aug', 'Sep-Oct-Nov') illustrating the deepest Jun-Jul-Aug trough mirroring peak SSTs and transitional patterns in Mar-Apr-May/Sep-Oct-Nov. Informs how SST anomalies precede pressure changes.
""",

    "MSL Seasonal Timeseries": """
A time series plot showing the monthly climatological cycle of MSL. Noticeable dips in June–September reflect the monsoon low-pressure system, while peaks in December–February indicate the winter high-pressure regime. This seasonal signature is foundational for interpreting SST–MSL coupling.
""",

    "TCC Average Heatmap": """
Climatological mean total cloud cover (TCC) heatmap. Areas with frequent cloudiness align with warm SST and low-pressure zones. This baseline helps assess how SST anomalies shift cloud patterns, affecting radiation and precipitation.
""",

    "TCC Area Average Timeseries": """
Time series of area-averaged TCC over the Arabian Sea - Konkan grid. It peaks sharply during the monsoon months and dips in winter, mirroring MSL and SST cycles. 
""",

    "TCC Seasonal Heatmaps": """
Seasonal cloud cover maps showing maximum JJA cloudiness driven by strong SST‐induced convection and sparse DJF cover. Supports season-specific model training.
""",

    "TCC Seasonal Timeseries": """
Line plot of monthly mean TCC. The smooth seasonal curve illustrates the build-up and decline of cloudiness, providing a temporal context for when SST anomalies have maximum impact on cloud formation.
""",

    "TCC Seasonal Mean Barplot": """
Panel of seasonal mean TCC maps. The monsoon season (Jun-Jul-Aug) exhibits maximum cloud cover, reflecting high SST-driven convection. Winter (Dec-Jan-Feb) shows sparse clouds, reinforcing the seasonality captured in your correlation analysis.
""",

    "U10 Average Heatmap": """
Mean ECMWF 10 m zonal wind (u10) heatmap. Positive values indicate easterly winds dominating along the coast. Comparing this with SST reveals how thermal gradients drive zonal wind shifts relevant for moisture transport.

""",

    "U10 Area Average Timeseries": """
Area-averaged u10 time series, highlighting seasonal peaks during pre-monsoon and early monsoon, when SSTs are rising. These wind bursts enhance ocean mixing, impacting SST anomalies.
""",

    "U10 Seasonal Heatmaps": """
Four-panel seasonal maps of u10. The Jun-Jul-Aug panel shows strongest easterly flow; 
Dec-Jan-Feb shows weak westerlies or neutral winds.
""",

    "U10 Seasonal Timeseries": """
Line graph of the seasonal climatology of u10. It shows a clear increase in easterlies from March to June, followed by a decline post-monsoon.
""",

    "V10 Average Heatmap": """
Mean meridional wind (v10) averaged over 2019–2024. Positive (northward) flows dominate, with pockets of southward winds aligning with cooler SST zones. This map sets the baseline for v10–SST correlation studies.
""",

    "V10 Area Average Timeseries": """
Area-averaged v10 time series, showing strong northward winds in the monsoon months. Peaks in July correspond with peak SSTs.

""",

    "V10 Seasonal Heatmaps": """
Four-panel seasonal distribution of v10. Jun-Jul-Aug shows the strongest uniform northward wind band, while Dec-Jan-Feb  reveals scattered weaker flows.
""",

    "V10 Seasonal Timeseries": """
Seasonal cycle of v10 plotted as a smooth monthly climatology. It highlights how northward flow intensifies sharply in June and wanes by October.
""",

    "Wind Speed Area Average": """
Shows the spatial average of ECMWF 10-meter wind speed over the Konkan region from 2019 to 2024. By smoothing out local variability, it highlights seasonal peaks—stronger winds during monsoon onset (June–July) and calmer conditions in winter. Useful for benchmarking how wind energy potential and mixing in coastal waters relate to SST-driven changes.
""",

    "Wind Vectors Average Heatmap": """
Overlays average wind vectors (combining u10 and v10) atop mean SST contours. The arrows illustrate prevailing wind directions and magnitudes, revealing coastal upwelling zones where winds run parallel to the shore. Linking this with SST patterns helps interpret how wind-driven mixing affects regional sea surface temperatures.
"""
}

# --- Cyclone Data ---
CYCLONES_DATA = [
    {"name": "Cyclone Kyarr (2019)", "region": "Western India, Oman, Yemen, Somalia", "image": CYCLONE_IMAGES_DIR / "kyaar.jpeg", "explanation": "Cyclone Kyarr (October 2019) skirted..."},
    {"name": "Cyclone Nisarga (2020)", "region": "Maharashtra (Alibag),Maharashtra,", "image": CYCLONE_IMAGES_DIR / "nisarga.jpg", "explanation": "Cyclone Nisarga (June 2020) made landfall..."},
    {"name": "Cyclone Tauktae (2021)", "region": "Gujarat, Maharashtra, Karnataka, Goa, Kerala", "image": CYCLONE_IMAGES_DIR / "tauktae.png", "explanation": "Cyclone Tauktae (May 2021) tracked parallel..."},
    {"name": "Cyclone Biparjoy (2023)", "region": "Gujarat (Did not make landfall as expected, skirted coast)", "image": CYCLONE_IMAGES_DIR / "Biparjoy.jpg", "explanation": "Cyclone Biparjoy (June 2023) passed north..."}
]
DATASETS_TO_DOWNLOAD = {
    "NOAA Precipitation Data": NOAA_PRECIP_FILE, "NOAA SST Data": NOAA_SST_FILE, "ECMWF Combined Data (All Vars)": ECMWF_FILE,
}

# --- Helper Function to Run Plotting Scripts ---
def run_plotting_script(script_path, vars_info_for_script: dict, output_image_path, plot_type_name):
    if not script_path.exists():
        st.error(f"Plotting script not found: {script_path}")
        return False, None, f"Plotting script not found: {script_path}", None
    if not vars_info_for_script and plot_type_name == "Generic Two Variable Timeseries":
        st.error(f"Script '{script_path.name}' requires variable information, but none was provided.")
        return False, None, "Missing variable arguments for generic script.", None

    details_list_for_cmd = list(vars_info_for_script.values())
    script_internal_keys_cmd = list(vars_info_for_script.keys()) # These are file_keys from master_config

    unique_files_for_cmd = list(dict.fromkeys([str(d['data_file']) for d in details_list_for_cmd]))
    script_nc_varnames_for_cmd = [d['nc_var_name'] for d in details_list_for_cmd]
    script_display_names_for_cmd = [d['display_name'] for d in details_list_for_cmd]
    script_time_coords_for_cmd = [d['time_coord'] for d in details_list_for_cmd]
    
    cmd = [ sys.executable, str(script_path), "--output", str(output_image_path), "--plot_type", str(plot_type_name)]
    
    if script_internal_keys_cmd: cmd.extend(["--vars"] + script_internal_keys_cmd) # Contextual internal keys
    if unique_files_for_cmd: cmd.extend(["--files"] + unique_files_for_cmd)
    if script_nc_varnames_for_cmd: cmd.extend(["--varnames"] + script_nc_varnames_for_cmd)
    
    if plot_type_name == "Generic Two Variable Timeseries":
        if script_display_names_for_cmd: cmd.extend(["--displaynames"] + script_display_names_for_cmd)
        if script_time_coords_for_cmd: cmd.extend(["--timecoords"] + script_time_coords_for_cmd)
        if not (len(unique_files_for_cmd) in [1,2] and \
                len(script_nc_varnames_for_cmd) == 2 and \
                len(script_display_names_for_cmd) == 2 and \
                len(script_time_coords_for_cmd) == 2):
            err_msg = f"Argument mismatch for '{plot_type_name}'. Check selections."
            st.error(err_msg); return False, None, err_msg, None
            
    st.info(f"Running: `{' '.join(cmd)}`")
    try:
        process = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=300)
        stdout, stderr = process.stdout, process.stderr
        if process.returncode != 0: return False, stdout, stderr, None
        generated_csv_path = None
        if plot_type_name == "Correlation Matrix Plot":
            potential_csv_path = output_image_path.with_name(output_image_path.stem + "_data.csv")
            if potential_csv_path.exists(): generated_csv_path = potential_csv_path
            elif process.returncode == 0 : st.caption(f"Note: Corr Matrix script ran, but CSV not found: {potential_csv_path}")
        return True, stdout, stderr, generated_csv_path
    except FileNotFoundError: err_msg = f"Error: Python ({sys.executable}) or script path incorrect: {script_path}"; st.error(err_msg); return False, None, err_msg, None
    except subprocess.TimeoutExpired: err_msg = f"Error: Script {script_path.name} timed out."; st.error(err_msg); return False, "Script timed out.", err_msg, None
    except Exception as e: err_msg = f"Unexpected error running {script_path.name}: {e}"; st.error(err_msg); return False, str(e), err_msg, None

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
            if fp.exists():
                with open(fp, "rb") as file_content:
                    st.download_button(f"Download {dn} ({fp.name})", file_content, fp.name, "application/x-netcdf", key=f"dl_{dn.replace(' ', '_')}")
            else: st.caption(f"{dn} ({fp.name}) not found.")
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
        available_plot_scripts_sb = sorted([name for name, path in PLOT_SCRIPT_MAP.items() if path.exists()])
        if len(available_plot_scripts_sb) != len(PLOT_SCRIPT_MAP): st.caption("Note: Some plot scripts in config not found.")
        unavailable_plot_scripts_info_sb = []
    else: # Guided Mode
        analysis_type_selected_sb = st.radio( "Select Analysis Type:",
            ["Single Variable Analysis", "SST Combined Analysis", "Two Variable Combined Timeseries", "Correlation Matrix (All Variables)", "ECMWF All Variables Overview"],
            key="analysis_type_radio_guided", horizontal=True )
        if analysis_type_selected_sb == "Single Variable Analysis":
            if not USER_FRIENDLY_VARS_FILTERED: st.warning("No variables available.")
            else: selected_primary_var_ufk_sb = st.selectbox("Select Variable:", USER_FRIENDLY_VARS_FILTERED, format_func=lambda x: VARIABLE_DETAILS_MAP.get(x, {}).get('display_name', x.upper()), key="single_var_sel_sb")
        elif analysis_type_selected_sb == "SST Combined Analysis":
            if 'sst' not in USER_FRIENDLY_VARS_FILTERED: st.warning("SST data needed.")
            elif not VARS_FOR_SST_COMBINED_FILTERED: st.warning("No other vars for SST combination.")
            else: selected_secondary_var_ufk_for_sst_sb = st.selectbox("Combine SST with:", VARS_FOR_SST_COMBINED_FILTERED, format_func=lambda x: VARIABLE_DETAILS_MAP.get(x, {}).get('display_name', x.upper()), key="sst_comb_var_sel_sb")
        elif analysis_type_selected_sb == "Two Variable Combined Timeseries":
            st.caption("Select two different variables.")
            if len(USER_FRIENDLY_VARS_FILTERED) < 2: st.warning("Need at least two variables.")
            else:
                c1, ca, c2 = st.columns([5,1,5])
                with c1: selected_var1_ufk_for_combined_ts_sb = st.selectbox("Variable 1:", USER_FRIENDLY_VARS_FILTERED, format_func=lambda x: VARIABLE_DETAILS_MAP.get(x, {}).get('display_name', x.upper()), key="comb_ts_v1_sb")
                with ca: st.markdown("<br><p style='text-align:center;font-size:1.1em'>&</p>", unsafe_allow_html=True)
                with c2:
                    opts_v2 = [v for v in USER_FRIENDLY_VARS_FILTERED if v != selected_var1_ufk_for_combined_ts_sb] if selected_var1_ufk_for_combined_ts_sb else []
                    if selected_var1_ufk_for_combined_ts_sb and not opts_v2: st.warning("No other var for Variable 2.")
                    elif opts_v2: selected_var2_ufk_for_combined_ts_sb = st.selectbox("Variable 2:", opts_v2, format_func=lambda x: VARIABLE_DETAILS_MAP.get(x, {}).get('display_name', x.upper()), key="comb_ts_v2_sb")
        
        plot_candidates_guided = []
        current_mapping = ANALYSIS_PLOT_MAPPING.get(analysis_type_selected_sb)
        if current_mapping:
            if analysis_type_selected_sb == "Single Variable Analysis" and selected_primary_var_ufk_sb: plot_candidates_guided = current_mapping.get(selected_primary_var_ufk_sb, [])
            elif analysis_type_selected_sb == "SST Combined Analysis" and selected_secondary_var_ufk_for_sst_sb: plot_candidates_guided = current_mapping.get(selected_secondary_var_ufk_for_sst_sb, [])
            elif analysis_type_selected_sb == "Two Variable Combined Timeseries":
                if selected_var1_ufk_for_combined_ts_sb and selected_var2_ufk_for_combined_ts_sb and selected_var1_ufk_for_combined_ts_sb != selected_var2_ufk_for_combined_ts_sb:
                    plot_candidates_guided = ["Generic Two Variable Timeseries"]
            elif isinstance(current_mapping, list): plot_candidates_guided = current_mapping
        available_plot_scripts_sb = [p for p in plot_candidates_guided if p in PLOT_SCRIPT_MAP and PLOT_SCRIPT_MAP[p].exists()]
        unavailable_plot_scripts_info_sb = [p for p in plot_candidates_guided if p not in available_plot_scripts_sb]


    st.subheader("Select Plot Script(s) to Run")
    selected_plot_script_names_final = []
    if analysis_type_selected_sb == "Two Variable Combined Timeseries":
        if selected_var1_ufk_for_combined_ts_sb and selected_var2_ufk_for_combined_ts_sb and selected_var1_ufk_for_combined_ts_sb != selected_var2_ufk_for_combined_ts_sb:
            if "Generic Two Variable Timeseries" in PLOT_SCRIPT_MAP and PLOT_SCRIPT_MAP["Generic Two Variable Timeseries"].exists():
                selected_plot_script_names_final = ["Generic Two Variable Timeseries"]
                # var1_dn = VARIABLE_DETAILS_MAP[selected_var1_ufk_for_combined_ts_sb]['display_name']
                # var2_dn = VARIABLE_DETAILS_MAP[selected_var2_ufk_for_combined_ts_sb]['display_name']
                # st.success(f"Plotting: {var1_dn} & {var2_dn}") # This can be verbose
            else: st.warning("Script for 'Generic Two Variable Timeseries' missing.")
    elif available_plot_scripts_sb:
        selected_plot_script_names_final = st.multiselect("Choose from available plot scripts:", available_plot_scripts_sb, default=[], key="plot_script_ms_common")
    elif st.session_state.view_all_plots_mode: st.warning("No plot scripts configured or found.")
    elif analysis_type_selected_sb: st.info("No plot scripts available for current guided selection.")

    if unavailable_plot_scripts_info_sb and not st.session_state.view_all_plots_mode and analysis_type_selected_sb != "Two Variable Combined Timeseries":
        with st.expander("Unavailable/missing scripts for selection", expanded=False):
            for plot_name_unavail in sorted(list(set(unavailable_plot_scripts_info_sb))): st.caption(f"- {plot_name_unavail}")

# --- Cyclone Dialog --- (Assumed correct from previous version)
if st.session_state.get('show_cyclone_dialog', False):
    st.dialog(title="Major Cyclones Affecting the Region (2019-2023)")
    st.header("Notable Cyclones") 
    if CYCLONES_DATA:
        cyclone_tab_names = [cyc['name'] for cyc in CYCLONES_DATA]
        if cyclone_tab_names:
            cyclone_display_tabs = st.tabs(cyclone_tab_names) 
            for i, tab_item in enumerate(cyclone_display_tabs):
                with tab_item:
                    cyc_data_item = CYCLONES_DATA[i]
                    st.write(f"**Region of Impact:** {cyc_data_item['region']}")
                    if cyc_data_item['image'].exists():
                        try: st.image(Image.open(cyc_data_item['image']), caption=cyc_data_item['name'], use_container_width=True) 
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
    if st.session_state.view_all_plots_mode: st.info("⬅️ In 'View All Plots Mode'. Select plot scripts from sidebar.")
    else: st.info("⬅️ Complete selections in sidebar to generate visualizations.")
else:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        
        def determine_vars_for_script_call(plot_script_name, analysis_type, 
                                           primary_var_ufk, sst_secondary_ufk, 
                                           ts_var1_ufk, ts_var2_ufk, 
                                           is_view_all, master_cfg, var_details_map, user_vars_filt):
            vars_to_pass_dict = {} # This will map file_key to details
            
            if is_view_all:
                if plot_script_name == "Correlation Matrix Plot":
                    for ufk_ in user_vars_filt: vars_to_pass_dict[var_details_map[ufk_]['file_key']] = master_cfg[var_details_map[ufk_]['file_key']]
                elif plot_script_name == "ECMWF Multipanel Overview":
                    for ufk_ in user_vars_filt: 
                        if var_details_map[ufk_]['data_file'] == ECMWF_FILE: vars_to_pass_dict[var_details_map[ufk_]['file_key']] = master_cfg[var_details_map[ufk_]['file_key']]
                elif plot_script_name == "Generic Two Variable Timeseries":
                    if len(user_vars_filt) >= 2:
                        vars_to_pass_dict[var_details_map[user_vars_filt[0]]['file_key']] = master_cfg[var_details_map[user_vars_filt[0]]['file_key']]
                        vars_to_pass_dict[var_details_map[user_vars_filt[1]]['file_key']] = master_cfg[var_details_map[user_vars_filt[1]]['file_key']]
                    else: st.warning("Generic Timeseries needs 2+ vars in View All mode."); return {}
                else: # Infer for other types
                    found = False
                    for at, mapping in ANALYSIS_PLOT_MAPPING.items():
                        if isinstance(mapping, dict):
                            for p_var, p_list in mapping.items():
                                if plot_script_name in p_list:
                                    if at == "Single Variable Analysis" and p_var in var_details_map:
                                        vars_to_pass_dict[var_details_map[p_var]['file_key']] = master_cfg[var_details_map[p_var]['file_key']]
                                    elif at == "SST Combined Analysis" and 'sst' in var_details_map and p_var in var_details_map:
                                        vars_to_pass_dict[var_details_map['sst']['file_key']] = master_cfg[var_details_map['sst']['file_key']]
                                        vars_to_pass_dict[var_details_map[p_var]['file_key']] = master_cfg[var_details_map[p_var]['file_key']]
                                    found = True; break
                        if found: break
                    if not found: st.caption(f"'{plot_script_name}' in View All: Passing all vars or script is self-contained.")
                        # Fallback: pass all if specific inference fails (script must handle)
                        # for ufk_all_fallback in user_vars_filt: vars_to_pass_dict[var_details_map[ufk_all_fallback]['file_key']] = master_cfg[var_details_map[ufk_all_fallback]['file_key']]


            elif analysis_type == "Single Variable Analysis" and primary_var_ufk:
                vars_to_pass_dict[var_details_map[primary_var_ufk]['file_key']] = master_cfg[var_details_map[primary_var_ufk]['file_key']]
            elif analysis_type == "SST Combined Analysis" and sst_secondary_ufk:
                vars_to_pass_dict[var_details_map['sst']['file_key']] = master_cfg[var_details_map['sst']['file_key']]
                vars_to_pass_dict[var_details_map[sst_secondary_ufk]['file_key']] = master_cfg[var_details_map[sst_secondary_ufk]['file_key']]
            elif analysis_type == "Two Variable Combined Timeseries" and ts_var1_ufk and ts_var2_ufk:
                vars_to_pass_dict[var_details_map[ts_var1_ufk]['file_key']] = master_cfg[var_details_map[ts_var1_ufk]['file_key']]
                vars_to_pass_dict[var_details_map[ts_var2_ufk]['file_key']] = master_cfg[var_details_map[ts_var2_ufk]['file_key']]
            elif analysis_type == "Correlation Matrix (All Variables)":
                for ufk_ in user_vars_filt: vars_to_pass_dict[var_details_map[ufk_]['file_key']] = master_cfg[var_details_map[ufk_]['file_key']]
            elif analysis_type == "ECMWF All Variables Overview":
                for ufk_ in user_vars_filt:
                    if var_details_map[ufk_]['data_file'] == ECMWF_FILE: vars_to_pass_dict[var_details_map[ufk_]['file_key']] = master_cfg[var_details_map[ufk_]['file_key']]
            
            return vars_to_pass_dict

        if len(selected_plot_script_names_final) == 1:
            plot_script_name_to_run_single = selected_plot_script_names_final[0]
            
            # Determine title context
            title_context_single = ""
            if st.session_state.view_all_plots_mode: title_context_single = "(from All Plots list)"
            elif analysis_type_selected_sb == "Single Variable Analysis" and selected_primary_var_ufk_sb: title_context_single = f"for {VARIABLE_DETAILS_MAP[selected_primary_var_ufk_sb]['display_name']}"
            elif analysis_type_selected_sb == "SST Combined Analysis" and selected_secondary_var_ufk_for_sst_sb: title_context_single = f"for SST & {VARIABLE_DETAILS_MAP[selected_secondary_var_ufk_for_sst_sb]['display_name']}"
            elif analysis_type_selected_sb == "Two Variable Combined Timeseries" and selected_var1_ufk_for_combined_ts_sb and selected_var2_ufk_for_combined_ts_sb:
                title_context_single = f"for {VARIABLE_DETAILS_MAP[selected_var1_ufk_for_combined_ts_sb]['display_name']} & {VARIABLE_DETAILS_MAP[selected_var2_ufk_for_combined_ts_sb]['display_name']}"
            elif analysis_type_selected_sb: title_context_single = f"({analysis_type_selected_sb})"

            st.subheader(f"Displaying: {plot_script_name_to_run_single} {title_context_single}")
            script_to_run_single = PLOT_SCRIPT_MAP.get(plot_script_name_to_run_single)

            if not script_to_run_single or not script_to_run_single.exists(): st.error(f"Script for '{plot_script_name_to_run_single}' not found/configured.")
            else:
                vars_for_script_call_single = determine_vars_for_script_call(
                    plot_script_name_to_run_single, analysis_type_selected_sb,
                    selected_primary_var_ufk_sb, selected_secondary_var_ufk_for_sst_sb,
                    selected_var1_ufk_for_combined_ts_sb, selected_var2_ufk_for_combined_ts_sb,
                    st.session_state.view_all_plots_mode, master_config, VARIABLE_DETAILS_MAP, USER_FRIENDLY_VARS_FILTERED
                )
                if not vars_for_script_call_single and plot_script_name_to_run_single not in ["Correlation Matrix Plot", "ECMWF Multipanel Overview"]: # Some general plots might not need vars if they know what to load
                    # For most scripts, especially specific single/dual var ones, this is an issue.
                     if not (st.session_state.view_all_plots_mode and plot_script_name_to_run_single not in ["Generic Two Variable Timeseries"]): # Allow generic to proceed if it can handle empty
                        st.warning(f"Could not determine required variables for '{plot_script_name_to_run_single}'. Plot generation might fail.")
                        # Do not proceed if critical vars are missing for non-general plots in guided mode.
                        if not st.session_state.view_all_plots_mode : return_from_func = True # Bit of a hack, should ideally use st.stop() or better flow
                
                output_img_name_single = f"{plot_script_name_to_run_single.replace(' ', '_').lower()}_plot.png"
                if plot_script_name_to_run_single == "Generic Two Variable Timeseries" and selected_var1_ufk_for_combined_ts_sb and selected_var2_ufk_for_combined_ts_sb:
                    output_img_name_single = f"timeseries_{selected_var1_ufk_for_combined_ts_sb}_vs_{selected_var2_ufk_for_combined_ts_sb}.png"
                elif plot_script_name_to_run_single == "Generic Two Variable Timeseries" and st.session_state.view_all_plots_mode and len(USER_FRIENDLY_VARS_FILTERED) >= 2:
                    # Default name for view all generic
                    output_img_name_single = f"timeseries_{USER_FRIENDLY_VARS_FILTERED[0]}_vs_{USER_FRIENDLY_VARS_FILTERED[1]}.png"


                output_img_path_single = tmp_path / output_img_name_single
                
                success_s, stdout_s, stderr_s, csv_path_s = run_plotting_script(script_to_run_single, vars_for_script_call_single, output_img_path_single, plot_script_name_to_run_single)
                if stdout_s: 
                    with st.expander(f"Output: {script_to_run_single.name}", expanded=False): st.code(stdout_s)
                if stderr_s: 
                    with st.expander(f"Errors/Warnings: {script_to_run_single.name}", expanded=not success_s): st.code(stderr_s, language='bash')

                if success_s and output_img_path_single.exists():
                    st.image(str(output_img_path_single), caption=f"Generated: {plot_script_name_to_run_single}", use_container_width=True)
                    with open(output_img_path_single, "rb") as f_img_s: st.download_button(f"Download Plot", f_img_s, output_img_name_single, "image/png", key=f"dl_s_{plot_script_name_to_run_single}")
                    if csv_path_s and csv_path_s.exists():
                        with open(csv_path_s, "rb") as f_csv_s: st.download_button(f"Download CSV", f_csv_s, csv_path_s.name, "text/csv", key=f"dl_csv_s_{plot_script_name_to_run_single}")
                    st.markdown("---"); st.subheader("Plot Explanation")
                    exp_key_s = plot_script_name_to_run_single
                    if plot_script_name_to_run_single == "Generic Two Variable Timeseries":
                        var1_for_key, var2_for_key = None, None
                        if analysis_type_selected_sb == "Two Variable Combined Timeseries": var1_for_key, var2_for_key = selected_var1_ufk_for_combined_ts_sb, selected_var2_ufk_for_combined_ts_sb
                        elif st.session_state.view_all_plots_mode and len(USER_FRIENDLY_VARS_FILTERED) >=2: var1_for_key,var2_for_key = USER_FRIENDLY_VARS_FILTERED[0], USER_FRIENDLY_VARS_FILTERED[1]
                        if var1_for_key and var2_for_key: exp_key_s = create_two_var_key(var1_for_key, var2_for_key)
                    st.markdown(PLOT_EXPLANATIONS.get(exp_key_s, DEFAULT_PLOT_EXPLANATION))
                elif success_s: st.error(f"Script ran but output not found: {output_img_path_single}.")
                # else: Error message from run_plotting_script already shown via stderr expander

        else: # Multiple plots in tabs
            tabs_objects = st.tabs([f"📊 {name}" for name in selected_plot_script_names_final])
            for i, tab_item_disp in enumerate(tabs_objects):
                with tab_item_disp:
                    current_plot_script_name_tab = selected_plot_script_names_final[i]
                    
                    title_context_tab = "" # Similar logic as single plot for context
                    # ... (determine title_context_tab based on analysis_type_selected_sb etc. for this current_plot_script_name_tab) ...
                    # For brevity, skipping detailed context determination per tab here, but it would be similar to single plot.
                    if st.session_state.view_all_plots_mode: title_context_tab = "(from All Plots list)"
                    elif analysis_type_selected_sb: title_context_tab = f"({analysis_type_selected_sb})"


                    st.markdown(f"#### {current_plot_script_name_tab} {title_context_tab}")
                    script_to_run_tab = PLOT_SCRIPT_MAP.get(current_plot_script_name_tab)

                    if not script_to_run_tab or not script_to_run_tab.exists(): st.error(f"Script for '{current_plot_script_name_tab}' not found."); continue
                    
                    vars_for_script_call_tab = determine_vars_for_script_call(
                        current_plot_script_name_tab, analysis_type_selected_sb,
                        selected_primary_var_ufk_sb, selected_secondary_var_ufk_for_sst_sb,
                        selected_var1_ufk_for_combined_ts_sb, selected_var2_ufk_for_combined_ts_sb,
                        st.session_state.view_all_plots_mode, master_config, VARIABLE_DETAILS_MAP, USER_FRIENDLY_VARS_FILTERED
                    )
                    if not vars_for_script_call_tab and current_plot_script_name_tab not in ["Correlation Matrix Plot", "ECMWF Multipanel Overview"]:
                         if not (st.session_state.view_all_plots_mode and current_plot_script_name_tab not in ["Generic Two Variable Timeseries"]):
                            st.warning(f"Vars for '{current_plot_script_name_tab}' in tab missing."); continue

                    output_img_name_tab_str = f"{current_plot_script_name_tab.replace(' ', '_').lower()}_plot.png"
                    if current_plot_script_name_tab == "Generic Two Variable Timeseries" and selected_var1_ufk_for_combined_ts_sb and selected_var2_ufk_for_combined_ts_sb:
                        output_img_name_tab_str = f"timeseries_{selected_var1_ufk_for_combined_ts_sb}_vs_{selected_var2_ufk_for_combined_ts_sb}.png"
                    elif current_plot_script_name_tab == "Generic Two Variable Timeseries" and st.session_state.view_all_plots_mode and len(USER_FRIENDLY_VARS_FILTERED) >= 2:
                        output_img_name_tab_str = f"timeseries_{USER_FRIENDLY_VARS_FILTERED[0]}_vs_{USER_FRIENDLY_VARS_FILTERED[1]}.png"
                    
                    output_img_path_tab_tmp = tmp_path / output_img_name_tab_str
                    
                    success_t, stdout_t, stderr_t, csv_path_t = run_plotting_script(script_to_run_tab, vars_for_script_call_tab, output_img_path_tab_tmp, current_plot_script_name_tab)
                    if stdout_t:
                        with st.expander(f"Output: {script_to_run_tab.name}", expanded=False): st.code(stdout_t)
                    if stderr_t:
                        with st.expander(f"Errors/Warnings: {script_to_run_tab.name}", expanded=not success_t): st.code(stderr_t, language='bash')

                    if success_t and output_img_path_tab_tmp.exists():
                        st.image(str(output_img_path_tab_tmp), caption=f"Generated: {current_plot_script_name_tab}", use_container_width=True)
                        with open(output_img_path_tab_tmp, "rb") as f_img_t: st.download_button(f"Download Plot", f_img_t, output_img_name_tab_str, "image/png", key=f"dl_t_{current_plot_script_name_tab}")
                        if csv_path_t and csv_path_t.exists():
                            with open(csv_path_t, "rb") as f_csv_t: st.download_button(f"Download CSV", f_csv_t, csv_path_t.name, "text/csv", key=f"dl_csv_t_{current_plot_script_name_tab}")
                        st.markdown("---"); st.subheader("Plot Explanation")
                        exp_key_t = current_plot_script_name_tab
                        if current_plot_script_name_tab == "Generic Two Variable Timeseries":
                            var1_key_tab, var2_key_tab = None, None
                            if analysis_type_selected_sb == "Two Variable Combined Timeseries": var1_key_tab, var2_key_tab = selected_var1_ufk_for_combined_ts_sb, selected_var2_ufk_for_combined_ts_sb
                            elif st.session_state.view_all_plots_mode and len(USER_FRIENDLY_VARS_FILTERED) >=2: var1_key_tab,var2_key_tab = USER_FRIENDLY_VARS_FILTERED[0], USER_FRIENDLY_VARS_FILTERED[1]
                            if var1_key_tab and var2_key_tab: exp_key_t = create_two_var_key(var1_key_tab, var2_key_tab)
                        st.markdown(PLOT_EXPLANATIONS.get(exp_key_t, DEFAULT_PLOT_EXPLANATION))
                    elif success_t: st.error(f"Script ran but output not found: {output_img_path_tab_tmp}.")

st.markdown("---")
st.markdown("""
<style>.footer-left {text-align: left; font-size: 16px; margin-top: 50px;} .footer-left img {height: 16px; vertical-align: middle; margin-left: 5px;}</style>
<div class="footer-left">Developed by: Riya Parmar <a href="https://www.linkedin.com/in/riya-parmar-44890a26b/" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/c/ca/LinkedIn_logo_initials.png" alt="LinkedIn"></a></div>
""", unsafe_allow_html=True)
