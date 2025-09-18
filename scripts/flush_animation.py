import pandas as pd
import os
import matplotlib
import imageio_ffmpeg

# Set ffmpeg path explicitly
matplotlib.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()
matplotlib.use('Agg')  # Use non-interactive backend for file output

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter

# Optional: print the ffmpeg path to confirm
print("FFmpeg path:", imageio_ffmpeg.get_ffmpeg_exe())


animation.FFMpegWriter = animation.writers['ffmpeg']



def load_flush_dfs(folder='../outputs/flush_data'):
    # Get the absolute path to the folder relative to this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    full_path = os.path.join(script_dir, folder)

    if not os.path.exists(full_path):
        raise FileNotFoundError(f"Folder not found: {full_path}")

    dfs = {}
    for filename in os.listdir(full_path):
        if filename.endswith('.csv'):
            sensor_name = filename.replace('.csv', '')
            dfs[sensor_name] = pd.read_csv(os.path.join(full_path, filename))
    return dfs

script_dir = os.path.dirname(os.path.abspath(__file__))
output_path_default = os.path.join(script_dir,'..', 'outputs', 'sensor_animation.mp4')
output_path_default = os.path.normpath(output_path_default)

def animate_diurnal_bins(df, bin_size_days=3, sensor_name='Sensor', output_path = output_path_default):
    df = df.copy()
    df = df.dropna(subset=['DateTime', 'detrended_daily_mean'])
    df['DateTime'] = pd.to_datetime(df['DateTime'])
    df['hour'] = df['DateTime'].dt.hour
    df['bin'] = ((df['DateTime'].dt.floor('D') - df['DateTime'].dt.floor('D').min()) // pd.Timedelta(days=bin_size_days)).astype(int)

    grouped = df.groupby(['bin', 'hour'])['detrended_daily_mean'].mean().reset_index()
    bin_ids = sorted(grouped['bin'].unique())

    fig, ax = plt.subplots(figsize=(10, 6))
    line, = ax.plot([], [], lw=2)
    ax.set_xlim(0, 23)
    ax.set_ylim(grouped['detrended_daily_mean'].min(), grouped['detrended_daily_mean'].max())
    ax.set_xlabel('Hour of Day')
    ax.set_ylabel('Avg Detrended Value')
    title = ax.set_title('')

    def update(frame):
        bin_id = bin_ids[frame]
        subset = grouped[grouped['bin'] == bin_id]
        line.set_data(subset['hour'], subset['detrended_daily_mean'])
        title.set_text(f'{sensor_name} - Bin {bin_id}')
        return line, title

    anim = animation.FuncAnimation(fig, update, frames=len(bin_ids), interval=800, blit=True)
    writer = FFMpegWriter(fps=2, metadata=dict(artist=sensor_name), bitrate=1800)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    anim.save(output_path, writer=writer)
    print(f"Animation saved to {output_path}")


# usage
sensor_name = 'sapflow'
flush_dfs = load_flush_dfs()
script_dir = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.join(script_dir,'...', 'outputs', f'{sensor_name}_animation.mp4')
output_path = os.path.normpath(output_path)

animate_diurnal_bins(flush_dfs[sensor_name], bin_size_days=3, sensor_name=sensor_name, output_path=output_path)

