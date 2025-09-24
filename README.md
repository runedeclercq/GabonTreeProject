# **Gabon Tree Exploration Project**

This folder contains the data and Jupyter Notebook scripts for analyzing TreeWatch data (Natkon dendrometer + sap flow sensor), combined with TOMST sensor data and phenological image data.


## SUBFOLDERS 

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
- 6_rain_and_flowering.ipynb	— identify rain events and bloom periods


## RUNNING THE PROJECT 

### 🌐 Option 1: Google Colab (Easiest)
1. Open Colab: go to colab.research.google.com

2. Run all notebooks
- Click File > Open notebook > GitHub
- Enter: runedeclercq/GabonTreeProject
- Select: scripts/0_SETUP_GOOGLE_COLAB.ipynb
- Run the Code Blocks
	- Step 4: Executes all notebooks in the folder, generating outputs and figures in the colab runtime. These will appear in the Colab file browser (left window pane) under: GabonTreeProject > outputs
	- Step 2: Saves executed notebooks to your Google Drive under: GabonTreeProject > outputs

OR: Open the Notebooks individually to view code and outputs
- Open each notebook via File > Open notebook > GitHub > runedeclercq
- To check wether the repository is open and accessible, you should see the GabonTreeProject folder with subfolders in the colab folder (on the left navigation pane).
- You can view code and outputs.
(⚠️ RUNNING THE NOTEBOOKS INDIVIDUALLY IS NOT POSSIBLE YET, I'M WORKING ON IT)


⚠️ Colab sessions are isolated per tab. If your notebooks depend on outputs from previous ones, use %run to execute them sequentially in a single notebook, or run run_all.py (code block 4 in setup notebook).


### 💻 Option 2: Local Setup (VS Code + Python)
1. Clone the repository
git clone https://github.com/runedeclercq/GabonTreeProject.git
2. Create a Virtual Environment and activate it
3. Install Dependencies: use "pip install -r requirements_local.txt"
4. Run the Project
- Open notebooks in scripts/ folder and run them in order OR
- Run everything at once: python scripts/run_all.py

## CITATION: 
Bush ER, Jeffery K, Bunnefeld N, Tutin C, Musgrave R, Moussavou G, Mihindou V, Malhi Y, Lehmann D, Edzang Ndong J, Makaga L, Abernethy K. 2020. Rare ground data confirm significant warming and drying in western equatorial Africa. PeerJ 8:e8732 DOI 10.7717/peerj.8732
