# **Gabon Tree Exploration Project**

This folder contains the data and Jupyter Notebook scripts for analyzing TreeWatch data (Natkon dendrometer + sap flow sensor), combined with TOMST sensor data and phenological image data.


## FOLDER STRUCTURE 

Bush_climate_R_scripts/ -  R scripts by Bush et al. (2020)

scripts/                  -  Python decid package + Jupyter notebooks



## SUBFOLDERS DETAILS 

### Bush_climate_R_scripts 
Contains the R code and provided data in Bush_data/.

The R script Bush_climate_variables_writeout.R exports plotting data as CSVs for use in Python notebooks.

### data
contains all necessary data for running the scripts

### scripts 
Contains the decid package and the jupyter notebooks:

- 1_data_cleanup.ipynb 		— data cleaning, detrending, smoothing
- 2_diameter_comparison.ipynb 	— resampling, aggregation, seasonal comparison, phenology event detection
- 3_diurnal_changes.ipynb 	— explore diurnal variations during leaf flush periods
- 4_climate.ipynb 		— load Bush climate data and reproduce graphs
- 5_water_relations.ipynb 	— hydraulic measures (TWD, DGR, MDS)
- 6_phenocam_data_animation.ipynb — create and save MP4 animation of phenocam data
- 6_rain_and_flowering.ipynb	— identify rain events and bloom periods

## RUNNING THE PROJECT 
Recommended setup: Visual Studio Code (VSC) with a local Python environment and jupyter extension.
1. Clone the repo
2. Download the dataset from [link]
3. Put the data in `data/` folder
4. Create a new Python virtual environment in a local folder (e.g., .venv).
5. Install dependencies: 
	pip install -r requirements_local.txt
6. Open the Jupyter notebooks in the scripts/ folder and run them in order, or run `python run_all.py` to generate outputs



## CITATION: 
Bush ER, Jeffery K, Bunnefeld N, Tutin C, Musgrave R, Moussavou G, Mihindou V, Malhi Y, Lehmann D, Edzang Ndong J, Makaga L, Abernethy K. 2020. Rare ground data confirm significant warming and drying in western equatorial Africa. PeerJ 8:e8732 DOI 10.7717/peerj.8732
