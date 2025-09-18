#!/usr/bin/env python3
"""
GabonTreeProject Startup Script
- Sets up paths for decid_package
- Executes notebooks in order
- Handles large files gracefully
- Saves outputs to outputs/ folder
"""

import os
import sys
import urllib.request

# ---------------------------
# 1️⃣ Paths
# ---------------------------
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_PATH = os.path.join(REPO_ROOT, "scripts")
DATA_PATH = os.path.join(REPO_ROOT, "data")
OUTPUTS_PATH = os.path.join(REPO_ROOT, "outputs")

os.makedirs(OUTPUTS_PATH, exist_ok=True)
sys.path.append(SCRIPTS_PATH)  # make decid_package importable

# ---------------------------
# 2️⃣ Check for large files
# ---------------------------
phenocam_zip = os.path.join(DATA_PATH, "Phenocams.zip")
if not os.path.exists(phenocam_zip):
    print(f"⚠️  Phenocams.zip not found at {phenocam_zip}")
    print("Please download it manually or provide a direct download link.")
    # Optional automatic download (replace URL if available)
    # url = "YOUR_ONEDRIVE_OR_DROPBOX_LINK"
    # urllib.request.urlretrieve(url, phenocam_zip)

# ---------------------------
# 3️⃣ List of notebooks to execute
# ---------------------------
notebooks = [
    "scripts/1_data_cleanup.ipynb",
    "scripts/2_diameter_comparison.ipynb",
    "scripts/3_diurnal_changes.ipynb",
    "scripts/4_climate.ipynb",
    "scripts/5_water_relations.ipynb",
    "scripts/6_rain_and_flowering.ipynb",
]

# Optional animation notebook
phenocam_notebook = "scripts/6_phenocam_data_animation.ipynb"
if os.path.exists(phenocam_zip):
    notebooks.append(phenocam_notebook)
else:
    print(f"Skipping {phenocam_notebook} (requires Phenocams.zip)")

# ---------------------------
# 4️⃣ Execute notebooks sequentially
# ---------------------------
for nb in notebooks:
    nb_path = os.path.join(REPO_ROOT, nb)
    output_nb = os.path.join(OUTPUTS_PATH, f"executed_{os.path.basename(nb)}")
    print(f"🚀 Executing {nb_path} → {output_nb}")
    os.system(f'jupyter nbconvert --to notebook --execute "{nb_path}" --output "{output_nb}"')

print("✅ All notebooks executed. Outputs are in the outputs/ folder.")
