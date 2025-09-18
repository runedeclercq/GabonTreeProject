import os

# folder where executed notebooks will be saved
OUTPUTS_PATH = "outputs"
os.makedirs(OUTPUTS_PATH, exist_ok=True)

# list of notebooks
notebooks = [
    "scripts/1_data_cleanup.ipynb",
    "scripts/2_diameter_comparison.ipynb",
]

for nb in notebooks:
    output_nb = os.path.join(OUTPUTS_PATH, f"executed_{os.path.basename(nb)}")
    print(f"Executing {nb} ...")
    os.system(f"jupyter nbconvert --to notebook --execute {nb} --output {output_nb}")

