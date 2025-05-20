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
# To be updated after cleaning and organizing files
```

## Main Outputs

```
# Expected outputs include:  
# - Correlation heatmaps between SST and climate variables  
# - Time-series plots of SST trends  
# - Interactive maps showing spatial overlays of climate variables  
# - Dashboard screenshots (included here for visual reference)
```

## Project Status

**Complete**

## License

This project is licensed under the [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) License.

---

*Developed by \[Your Name], Space Applications Centre, ISRO.*

---

Let me know if you want me to help fill in the directory structure or main outputs once you share your folder.
