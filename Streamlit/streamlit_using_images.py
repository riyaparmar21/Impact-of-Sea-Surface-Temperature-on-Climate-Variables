#!/usr/bin/env python
# coding: utf-8

import streamlit as st
import xarray as xr # Not strictly needed if only displaying images, but often kept
import os
from pathlib import Path
from PIL import Image
import pandas as pd # Not strictly needed if only displaying images

# --- Page Configuration ---
st.set_page_config(page_title="SST Impact Analysis-images", layout="wide")

# --- File Definitions ---
BASE_DIR = Path(".")
NOAA_PRECIP_FILE = BASE_DIR / "precip_konkan.nc" # Kept for download feature
ECMWF_FILE = BASE_DIR / "ecmwf datas.nc"      # Kept for download feature
NOAA_SST_FILE = BASE_DIR / "noaa_sst_masked.nc"   # Kept for download feature
CYCLONE_IMAGES_DIR = BASE_DIR / "cyclone images"
PRE_GENERATED_PLOTS_DIR = BASE_DIR / "visualizations"

if not PRE_GENERATED_PLOTS_DIR.exists():
    PRE_GENERATED_PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    st.info(f"Created directory for pre-generated plots: '{PRE_GENERATED_PLOTS_DIR}'")

# --- Plot Image Mapping (Ensure these paths are correct and images exist) ---
PLOT_IMAGE_MAP = {
    # SST Combined Analysis Plots
    "SST-PRECIP Combined Analysis": PRE_GENERATED_PLOTS_DIR / "combinations/corr/precip and sst.png",
    "SST-2M TEMP Combined Analysis": PRE_GENERATED_PLOTS_DIR / "combinations/corr/2m temp and sst.png",
    "SST-MSL Combined Analysis": PRE_GENERATED_PLOTS_DIR / "combinations/corr/msl and sst.png",
    "SST-TCC Combined Analysis": PRE_GENERATED_PLOTS_DIR / "combinations/corr/tcc and sst.png",
    "SST-U10 Combined Analysis": PRE_GENERATED_PLOTS_DIR / "combinations/corr/u wind and sst.png",
    "SST-V10 Combined Analysis": PRE_GENERATED_PLOTS_DIR / "combinations/corr/v wind and sst.png",
    "SST-TCC Spatial Correlation": PRE_GENERATED_PLOTS_DIR / "combinations/corr/sst-tcc spatial corr w values.png",
    "Correlation Matrix Plot": {"image": PRE_GENERATED_PLOTS_DIR / "combinations/corr/corr_matrix.png"},
    "ECMWF Multipanel Overview": PRE_GENERATED_PLOTS_DIR / "ecmwf/ecmwf_all_variables_average_panel.png",
    "PRECIP Average Heatmap": PRE_GENERATED_PLOTS_DIR / "precip/noaa_precip_average_heatmap.png",
    "PRECIP Seasonal Heatmaps": PRE_GENERATED_PLOTS_DIR / "precip/noaa_precip_seasonal_averages.png",
    "PRECIP Timeseries": PRE_GENERATED_PLOTS_DIR / "precip/precip_timeseries.png",
    "PRECIP At Selected Grid Points": PRE_GENERATED_PLOTS_DIR / "precip/precip at selected grid points.png",
    "PRECIP Distribution Over Time": PRE_GENERATED_PLOTS_DIR / "precip/precip_distribution_time.png",
    "PRECIP Anomalies Timeseries": PRE_GENERATED_PLOTS_DIR / "precip/precip_anomalies_timeseries.png",
    "PRECIP 6-Point Moving Average": PRE_GENERATED_PLOTS_DIR / "precip/precip_moving_avg.png",
    "SST Average Heatmap": PRE_GENERATED_PLOTS_DIR / "sst noaa/noaa_sst_average_heatmap.png",
    "SST Anomalies Timeseries": PRE_GENERATED_PLOTS_DIR / "sst noaa/sst noaa anomalies (diff from climatological mean).png",
    "SST Area Average Timeseries": PRE_GENERATED_PLOTS_DIR / "sst noaa/sst noaa area avg.png",
    "SST Seasonal Heatmaps": PRE_GENERATED_PLOTS_DIR / "sst noaa/noaa_sst_seasonal_averages.png",
    "SST Seasonal Cycle Timeseries": PRE_GENERATED_PLOTS_DIR / "sst noaa/sst noaa seasonal cycle.png",
    "T2M Area Average Timeseries": PRE_GENERATED_PLOTS_DIR / "ecmwf/t2m/2m temp area avg.png",
    "T2M Monthly Mean Barplot": PRE_GENERATED_PLOTS_DIR / "ecmwf/t2m/2m temp mon mean bar plot.png",
    "T2M Seasonal Heatmaps": PRE_GENERATED_PLOTS_DIR / "ecmwf/t2m/seasonal avg heatmaps.png",
    "T2M Seasonal Timeseries": PRE_GENERATED_PLOTS_DIR / "ecmwf/2m temp seasonal cycle.png", # Original had this path typo, kept for consistency with original
    "MSL Average Heatmap": PRE_GENERATED_PLOTS_DIR / "ecmwf/msl/ecmwf_msl_average_heatmap.png",
    "MSL Area Average Timeseries": PRE_GENERATED_PLOTS_DIR / "ecmwf/msl/mean sea level pressure.png",
    "MSL Seasonal Heatmaps": PRE_GENERATED_PLOTS_DIR / "ecmwf/msl/msl seasonal avg heatmaps.png",
    "MSL Seasonal Timeseries": PRE_GENERATED_PLOTS_DIR / "ecmwf/msl/mean sea level pressure seasonal cycle.png",
    "TCC Average Heatmap": PRE_GENERATED_PLOTS_DIR / "ecmwf/tcc/ecmwf_tcc_average_heatmap.png",
    "TCC Area Average Timeseries": PRE_GENERATED_PLOTS_DIR / "ecmwf/tcc/total cloud cover (area avg).png",
    "TCC Seasonal Heatmaps": PRE_GENERATED_PLOTS_DIR / "ecmwf/tcc/tcc seasonal avg heatmaps.png",
    "TCC Seasonal Timeseries": PRE_GENERATED_PLOTS_DIR / "ecmwf/tcc/total cloud cover seasonal cycle.png",
    "TCC Seasonal Mean Barplot": PRE_GENERATED_PLOTS_DIR / "ecmwf/tcc/seasonal mean cloud cover.png",
    "U10 Average Heatmap": PRE_GENERATED_PLOTS_DIR / "ecmwf/u10/ecmwf_u10_average_heatmap.png",
    "U10 Area Average Timeseries": PRE_GENERATED_PLOTS_DIR / "ecmwf/u10/10m u wind (area avg).png",
    "U10 Seasonal Heatmaps": PRE_GENERATED_PLOTS_DIR / "ecmwf/u10/u10 seasonal avg heatmaps.png",
    "U10 Seasonal Timeseries": PRE_GENERATED_PLOTS_DIR / "ecmwf/u10/10m u wind seasonal cycle.png",
    "V10 Average Heatmap": PRE_GENERATED_PLOTS_DIR / "ecmwf/v10/ecmwf_v10_average_heatmap.png",
    "V10 Area Average Timeseries": PRE_GENERATED_PLOTS_DIR / "ecmwf/v10/10m v wind area avg.png",
    "V10 Seasonal Heatmaps": PRE_GENERATED_PLOTS_DIR / "ecmwf/v10/v10 seasonal avg heatmaps.png",
    "V10 Seasonal Timeseries": PRE_GENERATED_PLOTS_DIR / "ecmwf/v10/10m v wind seasonal cycle.png",
    "Wind Speed Area Average": PRE_GENERATED_PLOTS_DIR / "ecmwf/10m wind speed area avg.png",
    "Wind Vectors Average Heatmap": PRE_GENERATED_PLOTS_DIR / "ecmwf/ecmwf_wind_vectors_average.png",
}
PLOT_IMAGE_MAP = {k.strip(): v for k, v in PLOT_IMAGE_MAP.items()}

# --- PERMANENT PLOT EXPLANATIONS ---
DEFAULT_PLOT_EXPLANATION = "This pre-generated plot visualizes the selected data. (Further details about its generation and insights can be added here.)"
PLOT_EXPLANATIONS = {
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
Annotates SST vs. total cloud cover correlations with exact values (≈ 0.4–0.6 offshore). Warmer SSTs drive enhanced convection and cloud formation, a key feedback for radiation and precipitation processes that your model must capture.
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
Line plot of the seasonal cycle (monthly climatology) of 2 m temperature. Smooth curves highlight inter-month transitions, emphasizing how air temperature lags or leads SST—critical information for temporal feature engineering in your predictive models.
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
""",
}



# --- User-Friendly Variable Names ---
USER_FRIENDLY_VARS = sorted(['sst', 'precip', 'msl', 'tcc', 't2m', 'v10', 'u10'])
VARS_FOR_SST_COMBINED = sorted([v for v in USER_FRIENDLY_VARS if v != 'sst'])

# --- Define which plots belong to which analysis type ---
ANALYSIS_PLOT_MAPPING = {
    "Single Variable Analysis": {
        'precip': ["PRECIP Average Heatmap", "PRECIP Seasonal Heatmaps", "PRECIP Timeseries", "PRECIP At Selected Grid Points", "PRECIP Distribution Over Time", "PRECIP Anomalies Timeseries", "PRECIP 6-Point Moving Average"],
        'sst': ["SST Average Heatmap", "SST Anomalies Timeseries", "SST Area Average Timeseries", "SST Seasonal Heatmaps", "SST Seasonal Cycle Timeseries"],
        't2m': ["T2M Area Average Timeseries", "T2M Monthly Mean Barplot", "T2M Seasonal Heatmaps", "T2M Seasonal Timeseries"],
        'msl': ["MSL Average Heatmap", "MSL Area Average Timeseries", "MSL Seasonal Heatmaps", "MSL Seasonal Timeseries"],
        'tcc': ["TCC Average Heatmap", "TCC Area Average Timeseries", "TCC Seasonal Heatmaps", "TCC Seasonal Timeseries", "TCC Seasonal Mean Barplot"],
        'u10': ["U10 Average Heatmap", "U10 Area Average Timeseries", "U10 Seasonal Heatmaps", "U10 Seasonal Timeseries", "Wind Speed Area Average", "Wind Vectors Average Heatmap"],
        'v10': ["V10 Average Heatmap", "V10 Area Average Timeseries", "V10 Seasonal Heatmaps", "V10 Seasonal Timeseries", "Wind Speed Area Average", "Wind Vectors Average Heatmap"],
    },
    "SST Combined Analysis": {
        'precip': ["SST-PRECIP Combined Analysis"],
        't2m': ["SST-2M TEMP Combined Analysis"],
        'msl': ["SST-MSL Combined Analysis"],
        'tcc': ["SST-TCC Combined Analysis", "SST-TCC Spatial Correlation"],
        'u10': ["SST-U10 Combined Analysis"],
        'v10': ["SST-V10 Combined Analysis"],
    },
    "Correlation Matrix (All Variables)": ["Correlation Matrix Plot"],
    "ECMWF All Variables Overview": ["ECMWF Multipanel Overview"],
}
# Config Check
for analysis_type, var_plots in ANALYSIS_PLOT_MAPPING.items():
    if isinstance(var_plots, dict):
        for var, plot_list in var_plots.items():
            for plot_name in plot_list:
                if plot_name not in PLOT_IMAGE_MAP: st.error(f"Config Error (1): Plot '{plot_name}' in ANALYSIS_PLOT_MAPPING but not in PLOT_IMAGE_MAP.")
    elif isinstance(var_plots, list):
        for plot_name in var_plots:
            if plot_name not in PLOT_IMAGE_MAP: st.error(f"Config Error (2): Plot '{plot_name}' in ANALYSIS_PLOT_MAPPING but not in PLOT_IMAGE_MAP.")

# --- Cyclone Data (with permanent explanations) ---
CYCLONES_DATA = [
    {
        "name": "Cyclone Kyarr (2019)", 
        "region": "Western India, Oman, Yemen, Somalia", 
        "image": CYCLONE_IMAGES_DIR / "kyaar.jpeg",
        "explanation": """
Cyclone Kyarr was an extremely severe cyclonic storm that formed in the Arabian Sea in October 2019. It Cyclone Kyarr (October 2019) skirted the western Arabian Sea, raising SST anomalies (+0.8 °C) as it suppressed coastal upwelling. The associated deep low-pressure system (< 990 hPa) drew strong northerly v10 winds (exceeding 5 m/s) into the region. Although Kyarr stayed offshore, peripheral rainbands delivered light to moderate precipitation (50–80 mm/day) and increased total cloud cover by 15–25%, highlighting how distant storms still modulate Konkan’s climate variables.
        """
    },
    {
        "name": "Cyclone Nisarga (2020)", 
        "region": "Maharashtra (Alibag),Maharashtra,", 
        "image": CYCLONE_IMAGES_DIR / "nisarga.jpg",
        "explanation": """
Cyclone Nisarga (June 2020) made landfall near Mumbai, directly impacting the Konkan coast. Peak SSTs beneath the storm climbed by \~1 °C pre-landfall, fueling moist convection. Surface pressure plunged below 985 hPa, driving strong easterly u10 winds (> 6 m/s) and northerly v10 gusts. The system dumped over 300 mm of rain in 24 hours, creating extreme precipitation anomalies (+200 mm above climatology) and boosting cloud cover to near 100% during its passage.
        """
    },
    {
        "name": "Cyclone Tauktae (2021)", 
        "region": "Gujarat, Maharashtra, Karnataka, Goa, Kerala", 
        "image": CYCLONE_IMAGES_DIR / "tauktae.png",
        "explanation": """
Cyclone Tauktae (May 2021) tracked parallel to the Konkan shoreline, generating a pronounced SST cooling of \~0.7 °C due to vigorous mixing. The cyclone’s low-pressure eye (< 980 hPa) intensified regional pressure gradients, producing one of the season’s strongest winds—u10 and v10 peaked above 7 m/s. Heavy rainbands delivered 100–200 mm/day, creating sharp positive precipitation anomalies, while total cloud cover soared above 90%, underscoring the dramatic, multi-variable impacts.
        """
    },
    {
        "name": "Cyclone Biparjoy (2023)", 
        "region": "Gujarat (Did not make landfall as expected, skirted coast)", 
        "image": CYCLONE_IMAGES_DIR / "biparjoy.jpg",
        "explanation": """
Cyclone Biparjoy (June 2023) passed north of the Konkan coast, elevating SST by up to 1–1.5 °C in its wake due to ocean mixing and upwelling. The storm’s low-pressure core (< 996 hPa) deepened the regional pressure trough, enhancing onshore winds (u10 and v10 increased by \~2 m/s). Convection driven by the warm, moist air mass led to a short-lived but intense precipitation spike (\~150 mm/day) and 20–30% surge in cloud cover over coastal Maharashtra.
        """
    }
]
# --- Variable Details Mapping ---
VARIABLE_DETAILS_MAP = {
    'precip': {'file': NOAA_PRECIP_FILE, 'display_name': 'Precipitation'},
    'sst': {'file': NOAA_SST_FILE, 'display_name': 'Sea Surface Temperature'},
    'msl': {'file': ECMWF_FILE, 'display_name': 'Mean Sea Level pressure'},
    'tcc': {'file': ECMWF_FILE, 'display_name': 'Total Cloud Cover'},
    't2m': {'file': ECMWF_FILE, 'display_name': '2 meter temperature'},
    'v10': {'file': ECMWF_FILE, 'display_name': 'V wind'},
    'u10': {'file': ECMWF_FILE, 'display_name': 'U wind'},
}
DATASETS_TO_DOWNLOAD = {
    "NOAA Precipitation Data": NOAA_PRECIP_FILE,
    "NOAA SST Data": NOAA_SST_FILE,
    "ECMWF Combined Data (All Vars)": ECMWF_FILE,
}

# --- Main App ---
st.title("SST Impact Analysis on climate variables")
st.markdown("""
**Region:** Arabian Sea near Konkan Coast, Maharashtra (Approx. 16°N - 21°N, 69°E - 74°E)
""")
st.markdown("---")

# --- Sidebar ---
with st.sidebar:
    st.header("Tools & Info")
    if st.button("🌀 Cyclone Information (2019-2023)", key="cyclone_info_button_sidebar"):
        st.session_state.show_cyclone_dialog = not st.session_state.get('show_cyclone_dialog', False)

    with st.expander("📥 Download Datasets", expanded=False):
        for display_name, file_path in DATASETS_TO_DOWNLOAD.items():
            if file_path.exists():
                with open(file_path, "rb") as fp:
                    st.download_button(
                        label=f"Download {display_name} ({file_path.name})",
                        data=fp,
                        file_name=file_path.name,
                        mime="application/x-netcdf",
                        key=f"download_{display_name.replace(' ', '_')}"
                    )
            else:
                st.caption(f"{display_name} ({file_path.name}) not found at {file_path.resolve()}")
    
    st.markdown("---")
    st.header("Plotting Controls")

    if 'view_all_plots_mode' not in st.session_state:
        st.session_state.view_all_plots_mode = False 

    current_view_all_mode = st.session_state.view_all_plots_mode
    view_all_plots_toggle = st.checkbox(
        "🔓 View All Available Plot Images", 
        value=current_view_all_mode, 
        key="view_all_plots_toggle_cb"
    )

    if view_all_plots_toggle != current_view_all_mode:
        st.session_state.view_all_plots_mode = view_all_plots_toggle
        st.rerun() 

    if st.session_state.view_all_plots_mode:
        st.subheader("All Available Plots")
        st.caption("Select any plot(s) from the complete list of available images.")
        analysis_type = "View All Plots Mode" 
        selected_primary_var = None
        selected_secondary_var_for_sst = None
        available_plot_types_for_multiselect = sorted(list(PLOT_IMAGE_MAP.keys()))
        unavailable_plot_types_info = [] 
    else:
        analysis_type = st.radio(
            "Select Analysis Type:",
            options=["Single Variable Analysis", "SST Combined Analysis", "Correlation Matrix (All Variables)", "ECMWF All Variables Overview"],
            key="analysis_type_radio_guided",
            horizontal=True, 
        )
        selected_primary_var = None
        selected_secondary_var_for_sst = None

        if analysis_type == "Single Variable Analysis":
            selected_primary_var = st.selectbox(
                "Select Variable:",
                options=USER_FRIENDLY_VARS,
                format_func=lambda x: VARIABLE_DETAILS_MAP.get(x, {}).get('display_name', x.upper()),
                key="single_var_select_guided"
            )
        elif analysis_type == "SST Combined Analysis":
            selected_secondary_var_for_sst = st.selectbox(
                "Select Variable to Combine with SST:",
                options=VARS_FOR_SST_COMBINED,
                format_func=lambda x: VARIABLE_DETAILS_MAP.get(x, {}).get('display_name', x.upper()),
                key="sst_combined_var_select_guided"
            )
        elif analysis_type == "Correlation Matrix (All Variables)":
            st.caption(f"Shows correlation matrix for: {', '.join(USER_FRIENDLY_VARS)}.")
        elif analysis_type == "ECMWF All Variables Overview":
            st.caption("Shows a multipanel overview of key ECMWF variables.")

        available_plot_types_for_multiselect = []
        unavailable_plot_types_info = []
        if analysis_type == "Single Variable Analysis" and selected_primary_var:
            available_plot_types_for_multiselect = ANALYSIS_PLOT_MAPPING["Single Variable Analysis"].get(selected_primary_var, [])
            all_single_var_plots = [p for var_plots in ANALYSIS_PLOT_MAPPING["Single Variable Analysis"].values() for p in var_plots]
            unavailable_plot_types_info = [p for p in all_single_var_plots if p not in available_plot_types_for_multiselect and p in PLOT_IMAGE_MAP]
        elif analysis_type == "SST Combined Analysis" and selected_secondary_var_for_sst:
            available_plot_types_for_multiselect = ANALYSIS_PLOT_MAPPING["SST Combined Analysis"].get(selected_secondary_var_for_sst, [])
            all_sst_combined_plots = [p for var_plots in ANALYSIS_PLOT_MAPPING["SST Combined Analysis"].values() for p in var_plots]
            unavailable_plot_types_info = [p for p in all_sst_combined_plots if p not in available_plot_types_for_multiselect and p in PLOT_IMAGE_MAP]
        elif analysis_type == "Correlation Matrix (All Variables)":
            available_plot_types_for_multiselect = ANALYSIS_PLOT_MAPPING["Correlation Matrix (All Variables)"]
            unavailable_plot_types_info = [p for p in PLOT_IMAGE_MAP.keys() if p not in available_plot_types_for_multiselect]
        elif analysis_type == "ECMWF All Variables Overview":
            available_plot_types_for_multiselect = ANALYSIS_PLOT_MAPPING["ECMWF All Variables Overview"]
            unavailable_plot_types_info = [p for p in PLOT_IMAGE_MAP.keys() if p not in available_plot_types_for_multiselect]
        
        available_plot_types_for_multiselect = [p for p in available_plot_types_for_multiselect if p in PLOT_IMAGE_MAP]
        unavailable_plot_types_info = [p for p in unavailable_plot_types_info if p in PLOT_IMAGE_MAP and p not in available_plot_types_for_multiselect]

    st.subheader("Select Plot Type(s) to Display")
    selected_plot_types_list = []
    if available_plot_types_for_multiselect:
        selected_plot_types_list = st.multiselect(
            "Choose from available plot types:",
            options=available_plot_types_for_multiselect, 
            default=[],
            key="plot_type_multiselect_common"
        )
    elif st.session_state.view_all_plots_mode: 
        st.warning("No plot images found in `PLOT_IMAGE_MAP` configuration.")
    elif analysis_type: 
        st.info("No specific plot types configured or available for the current guided selection.")
    else: 
        st.info("Select an analysis mode and relevant options to see plots.")

    if unavailable_plot_types_info and not st.session_state.view_all_plots_mode: 
        with st.expander("Other plot types (unavailable for current guided selection)", expanded=False):
            for plot_name_exp in unavailable_plot_types_info: # Renamed plot_name to avoid conflict
                st.caption(f"- {plot_name_exp}")

# --- Cyclone Dialog (Corrected for older Streamlit versions) ---
if st.session_state.get('show_cyclone_dialog', False):
    st.dialog(title="Major Cyclones Affecting the Region (2019-2023)") 

    st.header("Notable Cyclones") 
    if CYCLONES_DATA:
        cyclone_tab_names = [cyc['name'] for cyc in CYCLONES_DATA]
        if cyclone_tab_names:
            cyclone_display_tabs = st.tabs(cyclone_tab_names) 
            for i, tab in enumerate(cyclone_display_tabs):
                with tab: 
                    cyc_data = CYCLONES_DATA[i]
                    st.write(f"**Region of Impact:** {cyc_data['region']}")
                    if cyc_data['image'].exists():
                        try:
                            pil_image = Image.open(cyc_data['image'])
                            st.image(pil_image, caption=cyc_data['name'], use_container_width=True) 
                        except Exception as e:
                            st.warning(f"Could not load image for {cyc_data['name']}: {e}") 
                    else:
                        st.warning(f"Image file not found: {cyc_data['image']}")
                    # Display permanent cyclone explanation
                    st.markdown("---") 
                    st.subheader("Details & Impact")
                    st.markdown(cyc_data.get("explanation", "Detailed explanation for this cyclone is pending."))
    else:
        st.write("No cyclone data configured.") 

    if st.button("Close Cyclone Info", key="close_cyclone_dialog_button_in_dialog_static"): # Unique key
        st.session_state.show_cyclone_dialog = False
        st.rerun()
# --- End Cyclone Dialog ---

st.header("Data Visualization")
explanation_context = "" # This variable is not strictly needed if explanations are static
# ... (explanation_context logic can be removed or simplified if not used by PLOT_EXPLANATIONS)

if not selected_plot_types_list:
    if st.session_state.view_all_plots_mode:
        st.info(f"⬅️ In 'View All Plots Mode'. Please select one or more plot types from the sidebar.")
    else:
        st.info(f"⬅️ Please complete selections in the sidebar (Analysis Type, Variables, Plot Types) to view visualizations.")
else:
    if len(selected_plot_types_list) == 1:
        plot_type_to_display = selected_plot_types_list[0]
        subheader_text = f"Displaying: {plot_type_to_display}"
        if not st.session_state.view_all_plots_mode: 
            subheader_text += f" (Analysis: {analysis_type})"
        st.subheader(subheader_text)
        
        plot_info_entry = PLOT_IMAGE_MAP.get(plot_type_to_display)
        image_path_to_display = None
        csv_path_for_download = None
        image_filename = "plot.png"

        if isinstance(plot_info_entry, dict): # For plots with associated CSV, like Correlation Matrix
            image_path_to_display = plot_info_entry.get("image")
            csv_path_for_download = plot_info_entry.get("csv") # Assuming you might add CSV paths here
        elif isinstance(plot_info_entry, Path):
            image_path_to_display = plot_info_entry
        
        if image_path_to_display:
            image_filename = image_path_to_display.name
            if image_path_to_display.exists():
                st.image(str(image_path_to_display), caption=f"Pre-generated: {plot_type_to_display}", use_container_width=True)
                with open(image_path_to_display, "rb") as file_img:
                    st.download_button(f"Download Plot ({image_filename})", file_img, image_filename, "image/png", key=f"dl_single_{plot_type_to_display.replace(' ', '_')}")
                if csv_path_for_download and Path(csv_path_for_download).exists(): # Ensure csv_path is Path object if string
                     with open(csv_path_for_download, "rb") as file_csv:
                        st.download_button(f"Download Data ({Path(csv_path_for_download).name})", file_csv, Path(csv_path_for_download).name, "text/csv", key=f"dl_csv_single_{plot_type_to_display.replace(' ', '_')}")
                
                # Display permanent plot explanation
                st.markdown("---")
                st.subheader("Plot Explanation")
                explanation_text = PLOT_EXPLANATIONS.get(plot_type_to_display, DEFAULT_PLOT_EXPLANATION)
                st.markdown(explanation_text)
            else:
                 st.error(f"Image file for '{plot_type_to_display}' not found at: {image_path_to_display}. Please ensure the file exists in '{PRE_GENERATED_PLOTS_DIR}'.")
        else:
            st.error(f"Configuration error for plot '{plot_type_to_display}'. Image path not defined in PLOT_IMAGE_MAP.")

    else: # Multiple plots selected, use tabs
        plot_tabs_titles = [f"📊 {pt}" for pt in selected_plot_types_list]
        if plot_tabs_titles:
            plot_tabs = st.tabs(plot_tabs_titles)
            for i, tab_content in enumerate(plot_tabs):
                with tab_content:
                    current_plot_type = selected_plot_types_list[i]
                    # tab_subheader_text = f"Displaying: {current_plot_type}" # Optional
                    # if not st.session_state.view_all_plots_mode:
                    #      tab_subheader_text += f" (Analysis: {analysis_type})"
                    # st.subheader(tab_subheader_text)

                    plot_info_entry = PLOT_IMAGE_MAP.get(current_plot_type)
                    image_path_to_display_tab = None
                    csv_path_for_download_tab = None
                    image_filename_tab = "plot.png"

                    if isinstance(plot_info_entry, dict):
                        image_path_to_display_tab = plot_info_entry.get("image")
                        csv_path_for_download_tab = plot_info_entry.get("csv")
                    elif isinstance(plot_info_entry, Path):
                        image_path_to_display_tab = plot_info_entry
                    
                    if image_path_to_display_tab:
                        image_filename_tab = image_path_to_display_tab.name
                        if image_path_to_display_tab.exists():
                            st.image(str(image_path_to_display_tab), caption=f"Pre-generated: {current_plot_type}", use_container_width=True)
                            with open(image_path_to_display_tab, "rb") as file_img:
                                st.download_button(label=f"Download Plot ({image_filename_tab})", data=file_img, file_name=image_filename_tab, mime="image/png", key=f"dl_img_tab_{current_plot_type.replace(' ', '_')}")
                            if csv_path_for_download_tab and Path(csv_path_for_download_tab).exists():
                                with open(csv_path_for_download_tab, "rb") as file_csv:
                                    st.download_button(label=f"Download Data ({Path(csv_path_for_download_tab).name})", data=file_csv, file_name=Path(csv_path_for_download_tab).name, mime="text/csv", key=f"dl_csv_tab_{current_plot_type.replace(' ', '_')}")
                            
                            # Display permanent plot explanation in tab
                            st.markdown("---")
                            st.subheader("Plot Explanation")
                            explanation_text_tab = PLOT_EXPLANATIONS.get(current_plot_type, DEFAULT_PLOT_EXPLANATION)
                            st.markdown(explanation_text_tab)
                        else:
                             st.error(f"Image file for '{current_plot_type}' not found at: {image_path_to_display_tab}. Please ensure the file exists in '{PRE_GENERATED_PLOTS_DIR}'.")
                    else:
                        st.error(f"Configuration error for plot '{current_plot_type}'. Image path not defined in PLOT_IMAGE_MAP.")
        else:
            st.info("No plot types selected to display in tabs.")

# --- Footer ---
st.markdown("---")
st.markdown("""
<style>.footer-left {text-align: left; font-size: 16px; margin-top: 50px;} .footer-left img {height: 16px; vertical-align: middle; margin-left: 5px;}</style>
<div class="footer-left">Developed by: Riya Parmar <a href="https://www.linkedin.com/in/riya-parmar-44890a26b/" target="_blank"><img src="https://upload.wikimedia.org/wikipedia/commons/c/ca/LinkedIn_logo_initials.png" alt="LinkedIn"></a></div>
""", unsafe_allow_html=True)