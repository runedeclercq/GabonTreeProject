#!/usr/bin/env python3
"""
GabonTreeProject Startup Script
- Sets up paths for decid_package
- Executes notebooks in order
- Handles large files gracefully
- Saves outputs to outputs/ folder or Drive
"""

import os
import sys
import urllib.request

# ---------------------------
# 1️⃣ Paths
# ---------------------------
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SCRIPTS_PATH = os.path.join(REPO_ROOT, "scripts")
DATA_PATH = os.path.join(REPO_ROOT, "data")

# Default outputs inside repo
OUTPUTS_PATH = os.path.join(REPO_ROOT, "outputs")

# If Google Drive is mounted, use it instead
drive_path = "/content/drive/MyDrive/GabonTreeProject/outputs"
if os.path.exists("/content/drive/MyDrive"):
    OUTPUTS_PATH = drive_path

os.makedirs(OUTPUTS_PATH, exist_ok=True)
sys.path.append(SCRIPTS_PATH)  # make decid_package importable

# print(REPO_ROOT)
# print(SCRIPTS_PATH)


# ---------------------------
# 2️⃣ Check for large files
# ---------------------------
# phenocam_zip = os.path.join(DATA_PATH, "Phenocams.zip")
# has_phenocam = True
# if not os.path.exists(phenocam_zip):
#     has_phenocam = False
#     print(f"⚠️  Phenocams.zip not found at {phenocam_zip}")
#     print("Skipping the animation notebook. You can download it manually or provide a link.")

# ---------------------------
# 3️⃣ List of notebooks to execute
# ---------------------------
notebooks = [
    "1_data_cleanup.ipynb",
    "2_diameter_comparison.ipynb",
    "3_diurnal_changes.ipynb",
    "4_climate.ipynb",
    "5_water_relations.ipynb",
    "6_rain_and_flowering.ipynb",
]

# Only include the animation notebook if Phenocams.zip exists
if has_phenocam:
    notebooks.append("6_phenocam_data_animation.ipynb")

# ---------------------------
# 4️⃣ Execute notebooks sequentially
# ---------------------------
for nb in notebooks:
    nb_path = os.path.join(SCRIPTS_PATH, nb)
    output_nb = os.path.join(OUTPUTS_PATH, f"executed_{os.path.basename(nb)}")
    print(f"/n🚀 Executing {nb_path} → {output_nb}/n")
    os.system(f'jupyter nbconvert --to notebook --execute "{nb_path}" --output "{output_nb}"')

print(f"/n/n✅ All notebooks executed. Outputs are in: {OUTPUTS_PATH}")
