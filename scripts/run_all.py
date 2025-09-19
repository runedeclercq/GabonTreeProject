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
from nbclient import NotebookClient
from nbformat import read, write
from pathlib import Path

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
if SCRIPTS_PATH not in sys.path:
    sys.path.insert(0, SCRIPTS_PATH)

print(REPO_ROOT)
print(SCRIPTS_PATH)


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
    "scripts/1_data_cleanup.ipynb",
    "scripts/2_diameter_comparison.ipynb",
    "scripts/3_diurnal_changes.ipynb",
    "scripts/4_climate.ipynb",
    "scripts/5_water_relations.ipynb",
    "scripts/6_rain_and_flowering.ipynb",
]

# Only include the animation notebook if Phenocams.zip exists
# if has_phenocam:
#     notebooks.append("6_phenocam_data_animation.ipynb")

# ---------------------------
# 4️⃣ Execute notebooks sequentially
# ---------------------------
print("Scripts path:", os.path.join(SCRIPTS_PATH, "decid_package"))
print("Contents:", os.listdir(os.path.join(SCRIPTS_PATH, "decid_package")))


os.environ["PIPELINE_MODE"] = "1"

import asyncio
if sys.platform.startswith("win"):
    asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

REPO_ROOT = Path(REPO_ROOT)
for nb in notebooks:
    nb_path = REPO_ROOT / nb
    print(f"\n🚀 Executing {nb_path}")

    with open(nb_path, "r", encoding="utf-8") as f:
        nb_node = read(f, as_version=4)

    client = NotebookClient(
        nb_node,
        timeout=600,
        kernel_name="python3",
        resources={"metadata": {"path": str(nb_path.parent)}},  # THIS LINE
    )

    client.execute()

    # output_path = REPO_ROOT / "outputs" / f"executed_{Path(nb).stem}.ipynb"
    # with open(output_path, "w", encoding="utf-8") as f:
    #     write(nb_node, f)

print("\n\n✅ All notebooks executed.")

