import os
from datetime import datetime

import matplotlib
matplotlib.rcParams['font.family'] = 'Segoe UI Emoji'

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
from PIL import Image
import numpy as np
from tqdm.notebook import tqdm
import zipfile
import io

def extract_timestamp(filename):
    base_name = os.path.splitext(filename)[0]
    parts = base_name.split('_')
    if len(parts) < 5:
        raise ValueError(f"Invalid filename format: {filename}")
    date_str = parts[-2]
    time_str = parts[-1]
    timestamp_str = f"{date_str} {time_str.replace('-', ':')}"
    return datetime.strptime(timestamp_str, "%Y-%m-%d %H:%M:%S")


def filter_images_by_time(zip_path, start_time, end_time):
    filtered = []
    with zipfile.ZipFile(zip_path, 'r') as z:
        for fname in z.namelist():
            if fname.endswith(('.jpg', '.png')) and 'ERROR' not in fname:
                try: 
                    ts = extract_timestamp(fname)
                    if start_time <= ts <= end_time:
                        filtered.append(fname)
                except Exception as e:
                    print(f"Skipping {fname}: {e}")
    return filtered


def animate_images(zip_path, image_files, interval=500, export_path=None):
    if not image_files:
        print("No images found for the selected time range.")
        return None

    import zipfile, io, numpy as np, matplotlib.pyplot as plt
    from PIL import Image
    import matplotlib.animation as animation
    from tqdm import tqdm

    z = zipfile.ZipFile(zip_path, 'r')  # keep zip open

    first_image = Image.open(io.BytesIO(z.read(image_files[0]))).convert("RGB")
    w, h = first_image.size

    fig, ax = plt.subplots(figsize=(w / 100, h / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    img_display = ax.imshow(np.array(first_image))

    pbar = tqdm(total=len(image_files), desc="Animating", unit="frame")

    def update(i, zip_obj=z):
        im = Image.open(io.BytesIO(zip_obj.read(image_files[i]))).convert("RGB")
        img_display.set_data(np.array(im))
        pbar.update(1)
        return (img_display,)

    ani = animation.FuncAnimation(fig, update, frames=len(image_files),
                                  interval=interval, blit=True, repeat=False)

    if export_path:
        writer = animation.FFMpegWriter(fps=1000 // interval, bitrate=5000)
        ani.save(export_path, writer=writer, dpi=300)
        print(f"🎞️ Animation saved to '{export_path}'")

    pbar.close()
    z.close()
    plt.close(fig)

    return ani


def parse_time_string(time_str):
    return datetime.strptime(time_str, "%Y-%m-%d %H:%M:%S")


def run_phenocam_animation(zip_path, start_str, end_str, interval=300,
                           export_path="outputs/phenocam_videos/phenocam.mp4"):

    zip_path = os.path.abspath(zip_path)
    export_path = os.path.abspath(export_path)
    os.makedirs(os.path.dirname(export_path), exist_ok=True)

    start_time = parse_time_string(start_str)
    end_time = parse_time_string(end_str)

    filtered_images = filter_images_by_time(zip_path, start_time, end_time)
    print(f"Found {len(filtered_images)} images between {start_time} and {end_time}")

    animate_images(zip_path, filtered_images, interval, export_path=export_path)
