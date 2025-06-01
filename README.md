# Impact-of-Sea-Surface-Temperature-on-Climate-Variables
Impact of Sea Surface Temperature (SST) on climate variables in the Arabian Sea near the Konkan coast of Maharashtra, using global reanalysis and observational datasets. Features correlation analysis and an interactive Streamlit dashboard for visualization.


## Key Features & Achievements 

* Correlation analysis between SST and climate variables (precipitation, mean sea level pressure, wind components, temperature, cloud cover).
* Interactive dashboard built with Streamlit for exploratory data visualization.
* Comprehensive data processing pipeline using Python and xarray for NetCDF datasets.

## Data Sources & Acknowledgements

* **NOAA Precipitation**: National Centers for Environmental Information (NCEI) Climate Data Record (CDR) of global precipitation (2019–2025).
* **ECMWF Reanalysis (ERA5)**: Copernicus Climate Change Service (C3S) ERA5 monthly averaged fields for sea level pressure, 2 m temperature, wind components, and cloud cover.
* **NOAA SST**: NOAA Optimum Interpolation SST (OISST) dataset.
* **MOSDAC SST**: Meteorological and Oceanographic Satellite Data Archival Centre SST product (Note: contains missing values).

*Please ensure proper citation and attribution when using these datasets in your work.*

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/riyaparmar21/Impact-of-Sea-Surface-Temperature-on-Climate-Variables.git
cd Impact-of-Sea-Surface-Temperature-on-Climate-Variables
```

### 2. Create & Activate a Virtual Environment

Use Python 3.9 (recommended):

```bash
python3 -m venv venv
source venv/bin/activate      # on Linux/macOS
venv\Scripts\activate         # on Windows
```

### 3. Install Dependencies

```bash
pip install -r requirements.txt
```

## Usage

Launch the Streamlit dashboard:

```bash
streamlit run <streamlit_using_scripts>.py
```

### Note on Streamlit Apps

This repository includes two Streamlit applications:

* `streamlit_using_scripts.py`: Uses .py scripts to generate plots. The full-featured dashboard with all plots and functionalities. Recommended for comprehensive analysis.
* `streamlit_using_images.py`: Showcases the visualizations images stored locally. A simplified app showing fewer plots and limited functionalities. Not recommended for full usage.

## Dependencies

* numpy
* pandas
* xarray
* rioxarray
* matplotlib
* seaborn
* scikit-learn
* streamlit
* netCDF4
* h5netcdf
* plotly
* scipy
* cmocean
* os
* cartopy
  

## Directory Structure

```
├── datasets/
│   ├── ecmwf datas.nc      (used for final plotting and dashboard making)
│   ├── noaa_konkan_sst.nc
│   ├── noaa_sst_masked.nc  (used for final plotting and dashboard making)
│   ├── precip_konkan.nc    (used for final plotting and dashboard making)
│   └── precip.mon.mean.nc
│
├── notebooks/
│   ├── checking data redundancy and need for interpolation.ipynb
│   ├── Combined variables visualizations.ipynb
│   ├── comparing sst jan datasets.ipynb
│   ├── correlation matrix.ipynb
│   ├── cross-correlation.ipynb
│   ├── crs_reprojection.ipynb
│   ├── deleting variables.ipynb
│   ├── mosdac data preprocessing.ipynb
│   ├── NOAA SST nan values masking - for visualizations.ipynb
│   ├── preprocess lst, air temp and precip data.ipynb
│   ├── slicing and extracting noaa sst data.ipynb
│   ├── SST NOAA data extraction and cleaning.ipynb
│   └── Visualizations - individual datasets.ipynb
│
├── plotting_scripts_for_images/    (this scripts are used for streamlit_using_images.py file)
│   ├── 2m temp area avg time series.py
│   ├── 2m temp avg heatmap.py
│   ├── 2m temp monthly mean bar plot.py
│   ├── 2m temp seasonal avg heatmaps.py
│   ├── 2m temp seasonal time series.py
│   ├── 10m wind speed area avg.py
│   ├── correlation matrix plot.py
│   ├── msl area avg time series.py
│   ├── msl avg heatmap.py
│   ├── msl seasonal avg heatmaps.py
│   ├── msl seasonal time series.py
│   ├── multipanel ecmwf variables avg heatmaps.py
│   ├── precip seasonal avg heatmaps.py
│   ├── precipitation area averaged heatmap.py
│   ├── precipitation area avg time series.py
│   ├── precipitation at selected grid points.py
│   ├── precipitation seasonal time series.py
│   ├── seasonal mean cloud cover bar plot.py
│   ├── sst - msl combined analysis.py
│   ├── sst - t2m combined analysis.py
│   ├── sst - tcc combined analysis.py
│   ├── sst - tcc spatial correlation with values.py
│   ├── sst - u10 combined analysis.py
│   ├── sst - v10 combined analysis.py
│   ├── SST anomalies time series.py
│   ├── SST area avg time series.py
│   ├── SST avg heatmap.py
│   ├── SST seasonal avg heatmaps.py
│   ├── SST seasonal cycle time series.py
│   ├── sst-precip combined analysis.py
│   ├── total cloud cover area avg time series.py
│   ├── total cloud cover avg heatmap.py
│   ├── total cloud cover seasonal avg heatmaps.py
│   ├── total cloud cover seasonal time series.py
│   ├── u wind area avg time series.py
│   ├── u wind seasonal time series.py
│   └── u-10 wind avg heatmap.py
|   |-- u-10 wind seasonal avg heatmaps.py
|   |-- v wind area avg time series.py
|   |-- v wind seasonal time series.py
|   |-- v-10 wind avg heatmap.py
|   |-- v-10 wind seasonal avg heatmaps.py
|   |-- wind vectors avg heatmap.py
│
├── streamlit/
│   ├── cyclone images/
│   │   ├── biparjoy.jpg
│   │   ├── kyaar.jpeg
│   │   ├── nisarga.jpg
│   │   └── tauktae.png
│   │
│   ├── dashboard ss/
│   ├── dashboard ss1/
│   │
│   ├── plotting_scripts/      (this scripts are used for streamlit_using_scripts.py notebook)
│   ├── 2m temp area avg time series.py
│   ├── 2m temp avg heatmap.py
│   ├── 2m temp monthly mean bar plot.py
│   ├── 2m temp seasonal avg heatmaps.py
│   ├── 2m temp seasonal time series.py
│   ├── 10m wind speed area avg.py
│   ├── correlation matrix plot.py
│   ├── msl area avg time series.py
│   ├── msl avg heatmap.py
│   ├── msl seasonal avg heatmaps.py
│   ├── msl seasonal time series.py
│   ├── multipanel ecmwf variables avg heatmaps.py
│   ├── plot_combined_timeseries.py
│   ├── precip avg heatmap.py
│   ├── precip seasonal avg heatmaps.py
│   ├── precipitation area averaged heatmap.py
│   ├── precipitation area avg time series.py
│   ├── precipitation at selected grid points.py
│   ├── precipitation seasonal time series.py
│   ├── seasonal mean cloud cover bar plot.py
│   ├── sst - msl combined analysis.py
│   ├── sst - t2m combined analysis.py
│   ├── sst - tcc combined analysis.py
│   ├── sst - tcc spatial correlation with values.py
│   ├── sst - u10 combined analysis.py
│   ├── sst - v10 combined analysis.py
│   ├── SST anomalies time series.py
│   ├── SST area avg time series.py
│   ├── SST avg heatmap.py
│   ├── SST seasonal avg heatmaps.py
│   ├── SST seasonal cycle time series.py
│   ├── sst-precip combined analysis.py
│   ├── total cloud cover area avg time series.py
│   ├── total cloud cover avg heatmap.py
│   ├── total cloud cover seasonal avg heatmaps.py
│   ├── total cloud cover seasonal time series.py
│   ├── u wind area avg time series.py
│   ├── u wind seasonal time series.py
│   └── u-10 wind avg heatmap.py
|   |-- u-10 wind seasonal avg heatmaps.py
|   |-- v wind area avg time series.py
|   |-- v wind seasonal time series.py
|   |-- v-10 wind avg heatmap.py
|   |-- v-10 wind seasonal avg heatmaps.py
|   |-- wind vectors avg heatmap.py
│   │
│   ├── streamlit_using_images.py
│   └── streamlit_using_scripts.py           (use this one for all functionalities and plots)
│   │
│   ├── ecmwf datas.nc      (used for final plotting and dashboard making)
│   ├── noaa_konkan_sst.nc
│   ├── noaa_sst_masked.nc  (used for final plotting and dashboard making)
│   ├── precip_konkan.nc    (used for final plotting and dashboard making)
│   └── precip.mon.mean.nc
│
└── visualizations/
    ├── .DS_Store
    │
    ├── sst noaa/
    │   ├── .DS_Store
    │   ├── noaa_sst_seasonal_averages.png
    │   ├── sst noaa area avg.png
    │   ├── sst noaa anomalies (diff from climatological mean).png
    │   ├── noaa_sst_average_heatmap.png
    │   └── sst noaa seasonal cycle.png
    │
    ├── ecmwf/
    │   ├── .DS_Store
    │   ├── 10m wind speed area avg.png
    │   ├── avg wind speed vectors.png
    │   ├── ecmwf_wind_vectors_average.png
    │   ├── ecmwf_all_variables_average_panel.png
    │   └── tcc/  ← (nested folder for cloud cover)
    │       ├── seasonal mean cloud cover.png
    │       ├── total cloud cover seasonal cycle.png
    │       ├── ecmwf_tcc_average_heatmap.png
    │       ├── tcc seasonal avg heatmaps.png
    │       ├── total cloud cover (area avg).png
    │
    │   ├── msl/
    │   │   ├── msl seasonal avg heatmaps.png
    │   │   ├── mean sea level pressure.png
    │   │   ├── ecmwf_msl_average_heatmap.png
    │   │   └── mean sea level pressure seasonal cycle.png
    │   │
    │   ├── t2m/
    │   │   ├── 2m temp mon mean bar plot.png
    │   │   ├── 2m temp area avg.png
    │   │   ├── t2m seasonal avg heatmaps.png
    │   │   ├── 2m temp seasonal cycle.png
    │   │   └── ecmwf_t2m_average_heatmap.png
    │   │
    │   ├── v10/
    │   │   ├── v10 seasonal avg heatmaps.png
    │   │   ├── 10m v wind area avg.png
    │   │   ├── 10m v wind seasonal cycle.png
    │   │   └── ecmwf_v10_average_heatmap.png
    │   │
    │   └── u10/
    │       ├── 10m u wind seasonal cycle.png
    │       ├── ecmwf_u10_average_heatmap.png
    │       ├── u10 seasonal avg heatmaps.png
    │       └── 10m u wind (area avg).png
    │
    ├── combinations/
    │   ├── .DS_Store
    │   │
    │   └── 2 variables timeseries/
    │       ├── precip_vs_v10_timeseries.png
    │       ├── sst_vs_v10_timeseries.png
    │       ├── t2m_vs_u10_timeseries.png
    │       ├── msl_vs_precip_timeseries.png
    │       ├── msl_vs_sst_timeseries.png
    │       ├── msl_vs_u10_timeseries.png
    │       ├── sst_vs_tcc_timeseries.png
    │       ├── precip_vs_tcc_timeseries.png
    │       ├── msl_vs_t2m_timeseries.png
    │       ├── tcc_vs_u10_timeseries.png
    │       ├── t2m_vs_tcc_timeseries.png
    │       ├── msl_vs_v10_timeseries.png
    │       ├── tcc_vs_v10_timeseries.png
    │       ├── precip_vs_sst_timeseries.png
    │       ├── precip_vs_u10_timeseries.png
    │       ├── sst_vs_u10_timeseries.png
    │       └── sst_vs_t2m_timeseries.png

```

## Main Outputs


This project provides a comprehensive analysis of Sea Surface Temperature (SST) impacts on key climate variables over the Konkan coast region (approx. 16°N - 21°N, 69°E - 74°E), primarily utilizing NOAA and ECMWF datasets. The interactive Streamlit dashboard serves as the central platform for visualizing these analyses and exploring climate dynamics.

**Key outputs and insights you can expect from this dashboard include:**

### 1. Detailed Single Variable Climate Analysis:
Users can explore the characteristics of individual climate variables, including:
*   **Spatial Patterns:**
    *   **Average Heatmaps:** Visualize the time-averaged spatial distribution (e.g., `SST Average Heatmap`, `Precipitation Average Heatmap`, `MSL Average Heatmap`).
    *   **Seasonal Average Heatmaps:** Understand how spatial patterns change across meteorological seasons (Dec–Jan–Feb, Mar–Apr–May, Jun–Jul–Aug, Sep–Oct–Nov) for variables like SST, Precipitation, T2M, etc.
*   **Temporal Trends & Variability:**
    *   **Area-Averaged Time Series:** Track the evolution of spatially averaged variables over the analysis period (e.g., `SST Area Average Timeseries`, `Precipitation Timeseries`).
    *   **Seasonal Cycle Time Series:** Examine the climatological monthly progression for variables like SST and T2M, highlighting typical annual cycles.
    *   **Anomaly Time Series:** Identify deviations from the long-term mean (e.g., `SST Anomalies Timeseries`), crucial for understanding unusual climate events.
    *   **Specific Grid Point Analysis:** Investigate time series at selected grid points to understand local variations (e.g., `Precipitation At Selected Grid Points`).
*   **Statistical Distributions & Summaries:**
    *   **Monthly Mean Bar Plots:** Compare average monthly values (e.g., `T2M Monthly Mean Barplot`).
    *   **Seasonal Mean Bar Plots:** Compare average seasonal values (e.g., `TCC Seasonal Mean Barplot`).

### 2. SST-Driven Combined Variable Analysis:
The dashboard allows for in-depth exploration of how SST influences other critical climate parameters:
*   **SST vs. Other Variables (Precipitation, T2M, MSL, TCC, U10, V10):**
    *   **Combined Time Series Plots:** Simultaneously view the time series of SST and another selected variable on a twin-axis plot, resampled to a common monthly frequency for aligned comparison.
    *   **Scatter Plots & Correlation:** Accompanying scatter plots visualize the direct relationship, with a calculated Pearson correlation coefficient and p-value indicating the strength and significance of the linkage.
*   **SST-TCC Spatial Correlation:** A dedicated map showing the grid-point level temporal correlation between SST and Total Cloud Cover, revealing spatially varying relationships.

### 3. Inter-Variable Relationships (User-Defined Pairs):
A flexible analysis mode empowers users to explore relationships between *any two* selected climate variables:
*   **"Two Variable Combined Timeseries":** Select any pair of available variables (e.g., Precipitation vs. MSL, T2M vs. TCC) to generate:
    *   Aligned monthly mean time series on a twin-axis plot.
    *   A corresponding scatter plot with a trend line.
    *   The Pearson correlation coefficient and p-value.
    *   This allows for discovery of novel or less obvious interactions within the climate system of the region.

### 4. Holistic Climate System Overview:
*   **Correlation Matrix:** A comprehensive matrix displaying Pearson correlation coefficients between all primary climate variables, providing a quick overview of linear relationships across the dataset.
*   **ECMWF Multipanel Overview:** A single figure presenting time-averaged spatial maps of key ECMWF variables (MSL, TCC, T2M, Wind Speed & Vectors), offering a snapshot of the mean climate state.

### 5. Contextual Cyclone Information:
*   An interactive dialog provides details on major cyclones (Kyarr 2019, Nisarga 2020, Tauktae 2021, Biparjoy 2023) that impacted or were near the study region.
*   Includes:
    *   Region of impact.
    *   Satellite imagery.
    *   **Detailed textual explanations** on their characteristics, intensity, and observed or potential impacts on local climate variables (SST, precipitation, wind, etc.), providing context for observed anomalies in the data.

### 6. Data Accessibility:
*   Download buttons for the core NetCDF datasets (NOAA Precipitation, NOAA SST, ECMWF combined) used in the analyses.
*   Download capability for generated plot images and, where applicable (e.g., Correlation Matrix), the underlying numerical data in CSV format.

### 7. Customizable Visualizations:
*   The dashboard primarily uses dynamically generated plots based on user selections, ensuring visualizations are tailored to the specific analysis being performed.
*   **Permanent Plot Explanations:** Each plot type is accompanied by a detailed, pre-written explanation outlining what the plot shows, how to interpret it, and potential insights, making the dashboard a valuable learning and analysis tool.

---
This dashboard aims to facilitate a deeper understanding of the regional climate dynamics around the Konkan coast in the Arabian Sea, with a particular focus on the influence of Sea Surface Temperatures and the impact of significant cyclonic events.

## Project Status

**Complete**

## License

This project is licensed under the [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) License.

---

*Developed by Riya Parmar. This work originated from an internship project at the Space Applications Centre (SAC), ISRO.*

