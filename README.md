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

### outputs
contains subfolders needed for writing out csv files and figures for further use

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

### 🌐 Option 1: Google Colab (Easiest)
1. Open Colab 
Go to colab.research.google.com

2. Load the Setup Notebook
- Click File > Open notebook > GitHub
- Enter: runedeclercq/GabonTreeProject
- Select: scripts/0_SETUP_GOOGLE_COLAB.ipynb

3. Run the Code Blocks
- Step 1 and Step 3 are required
- Step 2 and Step 4 are optional
	- Step 4: Executes all notebooks in the folder, generating outputs and figures. These will appear in the Colab file browser under: GabonTreeProject > outputs
	- Step 2: Saves executed notebooks to your Google Drive under: GabonTreeProject > outputs

4. To Run Other Notebooks Individually
- After setup (steps 1 & 3), open other notebooks via: File > Open notebook > GitHub in another tabblad, it's essential that the 0_SETUP_GOOGLE_COLAB.ipynb notebook stays open and running to access the repository and data
- You can view code, outputs, and edit/run them freely.


Local setup: Visual Studio Code (VSC) with a local Python environment and jupyter extension.

2. Create a new Python virtual environment in a local folder (e.g., .venv).
3. Install dependencies: 
	pip install -r requirements_local.txt
4. Open the Jupyter notebooks in the scripts/ folder and run them in order, or run `python run_all.py` to generate outputs

### 💻 Option 2: Local Setup (VS Code + Python)
1. Clone the repository
git clone https://github.com/runedeclercq/GabonTreeProject.git
2. Create a Virtual Environment
Example using .venv folder:
python -m venv .venv
source .venv/bin/activate  # or .venv\Scripts\activate on Windows
3. Install Dependencies
pip install -r requirements_local.txt
4. Run the Project
- Open notebooks in scripts/ folder and run them in order OR
- Run everything at once: python scripts/run_all.py

## CITATION: 
Bush ER, Jeffery K, Bunnefeld N, Tutin C, Musgrave R, Moussavou G, Mihindou V, Malhi Y, Lehmann D, Edzang Ndong J, Makaga L, Abernethy K. 2020. Rare ground data confirm significant warming and drying in western equatorial Africa. PeerJ 8:e8732 DOI 10.7717/peerj.8732
