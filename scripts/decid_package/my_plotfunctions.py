import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import calendar
import matplotlib.cm as cm
import numpy as np
import matplotlib.patches as patches
import string

def get_color_maps():
    sensor_colors = {
        'TOMST': 'tab:orange',
        'Natkon': 'tab:blue',
        'Sap Flow': 'tab:green'
    }

    phen_colors = {
        'steady':      "#6c7375",   # darker gray-blue
        'pre-flush':   "#58b2d3",   # deeper sky blue
        'flush':       "#1c80be",   # stronger blue
        'post-flush':  "#043b7a"    # bold navy
    }


    bloom_colors = {
        'pre-bloom':     '#f6e3b4',
        'major_bloom':   '#f4b942',
        'inter_bloom':   '#e8a735',
        'minor_bloom':   '#d98c2b',
        'post-bloom':    '#c97a2b'
    }

#     bloom_colors = {
#     'pre-bloom':   '#f9d8a6',  # very light orange
#     'major_bloom': '#f4a261',  # vivid orange
#     'inter-bloom': '#e76f51',  # medium orange-red
#     'minor-bloom': '#d65a3b',  # deeper, muted orange
#     'post-bloom':  '#9c472b'   # dark, desaturated brown-orange
# }


    season_colors = {
        'wet': 'deepskyblue',
        'dry': 'gold'
    }

    return {
        'sensor': sensor_colors,
        'phenology': phen_colors,
        'bloom': bloom_colors,
        'season': season_colors
    }


def plot_sap_flow(df, datetime_col='DateTime', original_col='Sap_flow',
                  trend_col='Trend', detrended_col='Sap_detrended_shifted', title_suffix=''):
    # Ensure DateTime is a column
    if datetime_col in df.index.names:
        df = df.reset_index()

    plt.figure(figsize=(14,5))
    
    # Plot original sap flow
    plt.plot(df[datetime_col], df[original_col], label='Original', color='tab:grey', alpha=0.7)

    
    # Plot trend
    if trend_col in df.columns:
        plt.plot(df[datetime_col], df[trend_col], label=f'{trend_col}', color='tab:blue', linewidth=2, linestyle = '--')
    
    # Plot detrended sap flow
    if detrended_col in df.columns:
        plt.plot(df[datetime_col], df[detrended_col], label='Detrended', color='forestgreen', alpha=0.8)
    
    plt.xlabel('DateTime')
    plt.ylabel('Sap Flow (l/h)')
    # plt.title(f'Sap Flow: Original, Trend, and Detrended {title_suffix}')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('../outputs/figures/sap_flow_detrending.png', dpi=600, bbox_inches='tight')
    plt.show()


def plot_diurnal_profiles_by_category(sensor_dfs: dict, value_cols: dict, category: str, ax_title, color_map: dict = None, 
                                      save = False, filename = 'diurnal_profiles_by_category.png', uncertainty = False):
    """
    Plot diurnal variation (mean daily cycle) of sensor data grouped by a given category.

    Parameters:
        sensor_dfs: dict of DataFrames with DateTime index
        value_cols: dict of column names for detrended values
        category: str, one of 'season', 'period', 'phenology'
    """

    dfs = []
    for key, df in sensor_dfs.items():
        tmp = df.copy()
        tmp['Hour'] = tmp['DateTime'].dt.hour + tmp['DateTime'].dt.minute / 60
        tmp['Sensor'] = key.capitalize()
        tmp['Value'] = tmp[value_cols[key]]
        tmp['Category'] = tmp[category]
        dfs.append(tmp[['Hour', 'Sensor', 'Value', 'Category']])

    combined = pd.concat(dfs)
    diurnal_stats = combined.groupby(['Sensor', 'Category', 'Hour'])['Value'].agg(['mean', 'std']).reset_index()
    diurnal_stats.columns = ['Sensor', 'Category', 'Hour', 'mean', 'std']

    categories_sorted = sorted(diurnal_stats['Category'].dropna().unique())
    if color_map is None:
        colors = cm.viridis(np.linspace(0, 1, len(categories_sorted)))
    else:
        colors = color_map
    if isinstance(colors, dict):
        category_color_dict = colors
    else:
        category_color_dict = {cat: colors[i] for i, cat in enumerate(categories_sorted)}


    plt.figure(figsize=(6 * len(sensor_dfs), 6))
    sns.set_style("whitegrid")

    for i, sensor in enumerate(diurnal_stats['Sensor'].unique(), 1):
        plt.subplot(1, len(sensor_dfs), i)
        sensor_df = diurnal_stats[diurnal_stats['Sensor'] == sensor]
        for cat in categories_sorted:
            if cat == 'unknown':
                continue
            else:
                dfc = sensor_df[sensor_df['Category'] == cat]
                plt.plot(dfc['Hour'], dfc['mean'], color=category_color_dict[cat], label=cat, alpha = 1, zorder = 2)
                if uncertainty == True:
                    plt.fill_between(dfc['Hour'], dfc['mean'] - dfc['std'], dfc['mean'] + dfc['std'],
                                    color=category_color_dict[cat], alpha=0.1, edgecolor = 'none', 
                                    zorder = 1)

        # plt.title(f'Diurnal Profile by {category.capitalize()}: {sensor}')
        plt.xlabel('hour of the day')
        plt.xticks(np.arange(0, 25, 6))
        plt.title(f'{sensor}')
        plt.xlim(0, 24)
        plt.ylabel(ax_title[sensor])
        if i == 1:
            plt.legend(loc = 'lower left', frameon = False)
        plt.grid(True, linestyle='--', alpha=0.5)
        if len(sensor_dfs) > 1:
            plt.text(0.05, 0.95, f'({string.ascii_lowercase[i-1]})', 
                transform = plt.gca().transAxes, 
                ha='left', va='top', fontsize = 16, 
                label = None)
        
        plt.tight_layout()
    if save: 
        plt.savefig(f'../outputs/figures/{filename}', dpi=600, bbox_inches='tight')
    plt.show()


def subplot_diurnal_profiles_by_category(sensor: str, 
                                         category: str,
                                         ax: plt.Axes,
                                         ax_title: dict, 
                                         sensor_dfs= dict,
                                         value_col: str = 'detrended_daily_mean',
                                         color_map: dict = None, 
                                        #  save: bool = False, 
                                        #  filename : str = 'diurnal_profiles_by_category.png', 
                                         uncertainty : bool = False, 
                                         legend: bool = True):
    """
    Plot diurnal variation (mean daily cycle) of sensor data grouped by a given category.

    Parameters:
        sensor_dfs: dict of DataFrames with DateTime index
        value_cols: dict of column names for detrended values
        category: str, one of 'season', 'period', 'phenology'
    """

    df = sensor_dfs[sensor].copy()
    df['Hour'] = df['DateTime'].dt.hour + df['DateTime'].dt.minute / 60
    df['Sensor'] = sensor.capitalize()
    df['Value'] = df[value_col]
    df['Category'] = df[category]

    diurnal_stats = df.groupby(['Category', 'Hour'])['Value'].agg(['mean', 'std']).reset_index()
    categories_sorted = sorted(diurnal_stats['Category'].dropna().unique())

    if color_map is None:
        colors = cm.viridis(np.linspace(0, 1, len(categories_sorted)))
        category_color_dict = {cat: colors[i] for i, cat in enumerate(categories_sorted)}
    elif isinstance(color_map, dict):
        category_color_dict = color_map
    else:
        category_color_dict = {cat: color_map[i] for i, cat in enumerate(categories_sorted)}


    # Plot on passed axis
    for cat in categories_sorted:
        if cat == 'unknown':
            continue
        dfc = diurnal_stats[diurnal_stats['Category'] == cat]
        ax.plot(dfc['Hour'], dfc['mean'], color=category_color_dict[cat], label=cat, alpha=1, zorder=2)
        if uncertainty:
            ax.fill_between(dfc['Hour'],
                            dfc['mean'] - dfc['std'],
                            dfc['mean'] + dfc['std'],
                            color=category_color_dict[cat],
                            alpha=0.1,
                            edgecolor='none',
                            zorder=1)

        # plt.title(f'Diurnal Profile by {category.capitalize()}: {sensor}')
    ax.set_ylabel(ax_title[sensor])
    # ax.set_title(sensor)
    ax.set_xlim(0, 24)
    ax.grid(True, linestyle='--', alpha=0.5)
    if legend == True:
        ax.legend(loc='lower left', frameon = False)
        ax.set_xticks([])
    else:
        ax.set_xlabel('Hour of the day')
        ax.set_xticks(np.arange(0, 25, 6))




def plot_diurnal_profiles(sensor_dfs: dict, value_cols: dict, ax_titles: dict, 
                          save = False, filename = 'diurnal_profiles.png'):
    """
    Plot diurnal variation (mean daily cycle) of sensor data per month.

    Parameters:
        sensor_dfs: dict of DataFrames with DateTime index, e.g. {'natkon': df1, 'tomst': df2}
        value_cols: dict of column names for detrended values, e.g. {'natkon': 'detrended_daily_mean'}
    """
    dfs = []
    for key, df in sensor_dfs.items():
        tmp = df.copy()
        tmp['Hour'] = tmp['DateTime'].dt.hour + tmp['DateTime'].dt.minute / 60
        tmp['Month'] = tmp['DateTime'].dt.month
        tmp['Sensor'] = key.capitalize()
        tmp['Value'] = tmp[value_cols[key]]
        dfs.append(tmp[['Hour', 'Month', 'Sensor', 'Value']])

    combined = pd.concat(dfs)
    diurnal = combined.groupby(['Sensor', 'Month', 'Hour'])['Value'].mean().reset_index()

    months_sorted = sorted(diurnal['Month'].unique())
    colors = cm.viridis(np.linspace(0, 1, len(months_sorted)))
    month_color_dict = {month: colors[i] for i, month in enumerate(months_sorted)}
    month_labels = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']

    plt.figure(figsize=(6 * len(sensor_dfs), 6))
    sns.set_style("whitegrid")

    for i, sensor in enumerate(diurnal['Sensor'].unique(), 1):
        plt.subplot(1, len(sensor_dfs), i)
        sensor_df = diurnal[diurnal['Sensor'] == sensor]
        for month in sorted(sensor_df['Month'].unique()):
            dfm = sensor_df[sensor_df['Month'] == month]
            plt.plot(dfm['Hour'], dfm['Value'], color=month_color_dict[month],
                     label=month_labels[month - 1])
        # plt.title(f'Diurnal Profile: {sensor}')
        if len(sensor_dfs) > 1:
            plt.text(0.05, 0.95, f'({string.ascii_lowercase[i-1]})', 
                    transform = plt.gca().transAxes, 
                    ha='left', va='top', fontsize = 16, 
                    label = None)
        plt.xlabel('hour of the day')
        plt.xlim(0, 24)
        plt.xticks([0, 6, 12, 18, 24])
        plt.title(f'{sensor}')
        plt.ylabel(ax_titles[sensor])
        plt.grid(True, linestyle='--', alpha=0.5)
        if i == len(sensor_dfs):
            plt.legend(title = 'month',loc = 'lower right', ncol = 2)
        plt.tight_layout()
    if save:
        plt.savefig(f'../outputs/figures/{filename}', dpi=600, bbox_inches='tight')
    plt.show()


def plot_diurnal_by_season(sensor_dfs: dict, value_cols: dict, ax_titles: dict, 
                           save = False, filename = 'diurnal_seasons.png'):
    """
    Plot diurnal variation (mean daily cycle) of sensor data aggregated by tropical seasons.

    Parameters:
        sensor_dfs: dict of DataFrames with DateTime index
        value_cols: dict of column names for detrended values
    """
    season_dict = {
        10: 'ON', 11: 'ON',
        12: 'DJF', 1: 'DJF', 2: 'DJF',
        3: 'MAM', 4: 'MAM', 5: 'MAM',
        6: 'JJAS', 7: 'JJAS', 8: 'JJAS', 9: 'JJAS'
    }

    dfs = []
    for key, df in sensor_dfs.items():
        tmp = df.copy()
        tmp['Hour'] = tmp['DateTime'].dt.hour + tmp['DateTime'].dt.minute / 60
        tmp['Month'] = tmp['DateTime'].dt.month
        tmp['Season'] = tmp['Month'].map(season_dict)
        tmp['Sensor'] = key.capitalize()
        tmp['Value'] = tmp[value_cols[key]]
        dfs.append(tmp[['Hour', 'Season', 'Sensor', 'Value']])

    combined = pd.concat(dfs)
    diurnal = combined.groupby(['Sensor', 'Season', 'Hour'])['Value'].mean().reset_index()

    season_colors = {
        'JJAS': '#fdae61',  # Dry (June–Sept) – warm orange
        'DJF':  '#f46d43',  # Dry (Dec–Feb) – deeper orange
        'ON':   '#74add1',  # Wet (Oct–Nov) – soft blue
        'MAM':  '#4575b4'   # Wet (Mar–May) – deeper blue
    }


    plt.figure(figsize=(6 * len(sensor_dfs), 6))
    sns.set_style("whitegrid")

    for i, sensor in enumerate(diurnal['Sensor'].unique(), 1):
        ax = plt.subplot(1, len(sensor_dfs), i)
        sensor_df = diurnal[diurnal['Sensor'] == sensor]
        for season in ['ON', 'DJF', 'MAM', 'JJAS']:
            df_season = sensor_df[sensor_df['Season'] == season]
            if sensor == 'Sapflow':
                ax.plot(df_season['Hour'], df_season['Value'], color=season_colors[season],
                    label=season, linewidth=2.5)
            else: 
                ax.plot(df_season['Hour'], df_season['Value']*10, color=season_colors[season],
                    label=season, linewidth=2.5)
        # ax.set_title(f'Diurnal Profile: {sensor}')
        ax.set_xlabel('hour of the day')
        ax.set_ylabel(ax_titles[sensor])
        ax.set_xlim(0, 24)
        ax.grid(True, linestyle='--', alpha=0.5)
        if i == len(sensor_dfs):
            ax.legend(title='Season', loc = 'lower right')
        ax.set_xticks([0, 6, 12, 18, 24])
        ax.set_title(f'{sensor}')
        ax.set_xlim(0, 24)

        if len(sensor_dfs) > 1:
            ax.text(0.05, 0.95, f'({string.ascii_lowercase[i-1]})',
                    transform=ax.transAxes,
                    fontsize=16, va='top', ha='left')
    
    plt.tight_layout()
    if save is True:
        plt.savefig(f'../outputs/figures/{filename}', dpi=600, bbox_inches='tight')
    plt.show()



def main_plot(df_all_norm, sensor_colors, season_colors, agg_to_plot='20min', ):
    """
    Plots normalized stem diameter and sap flow data with seasonal shading.

    Parameters:
        df_all_norm (pd.DataFrame): DataFrame containing normalized sensor data.
        agg_to_plot (str): Aggregation type ('20min', 'roll10', or custom).
    """
    
    # --- Define seasons ---
    seasons = {
        'JJAS': {'start': -1, 'end': 120, 'type':'dry'},  # June–Sept
        'ON':   {'start': 121, 'end': 182, 'type':'wet'}, # Oct–Nov
        'DJF':  {'start': 183, 'end': 271, 'type':'dry'}, # Dec–Feb
        'MAM':  {'start': 272, 'end': 364, 'type':'wet'}  # Mar–May
    }
   
    
    # --- Column mapping ---
    if agg_to_plot == '20min':
        col_map = {s: 'smoothed_signal' for s in ['TOMST', 'Natkon', 'Sap Flow']}
    elif agg_to_plot == 'roll10d':
        col_map = {s: 'smoothed_signal_roll10d' for s in ['TOMST', 'Natkon', 'Sap Flow']}
    else:
        col_map = {s: f'smoothed_signal_mean_{agg_to_plot}' for s in ['TOMST', 'Natkon', 'Sap Flow']}
    
    # --- Create figure ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14,10), sharex=True)
    
    # --- Top subplot: dendrometer sensors ---
    for sensor in ['TOMST', 'Natkon']:
        col = col_map[sensor] + '_pct_period'
        subset = df_all_norm[df_all_norm['Sensor'] == sensor]
        
        for period in sorted(subset['Period'].unique()):
            period_subset = subset[subset['Period'] == period]
            
            if sensor == 'TOMST':
                flat_mask = period_subset['flat_flag'] == 1
                normal_mask = ~flat_mask
                
                # Normal data
                ax1.plot(
                    period_subset.loc[normal_mask, 'Day_of_Period'],
                    period_subset.loc[normal_mask, col],
                    color='tab:orange',
                    alpha=0.9 if period == 0 else 0.6,
                    label=f'TOMST Period {period+1} (normal)',
                    zorder=2 if period == 0 else 5
                )
                
                # Flat segments
                flat_points = period_subset.loc[flat_mask, ['Day_of_Period', col]]
                for _, seg in flat_points.groupby((flat_points['Day_of_Period'].diff() > 1).cumsum()):
                    ax1.plot(seg['Day_of_Period'], seg[col], color='red', zorder=3)
            
            else:
                ax1.plot(
                    period_subset['Day_of_Period'],
                    period_subset[col],
                    label=f'{sensor} Period {period+1}',
                    alpha=0.9 if period == 0 else 0.6,
                    color=sensor_colors[sensor],
                    zorder=2 if period == 0 else 5
                )
    
    # Season shading/labels
    for season, info in seasons.items():
        ax1.axvspan(info['start'], info['end'], color=season_colors[info['type']], alpha=0.05, zorder=0)
        ax1.axvline(x=info['start'], color='gray', linestyle='--', alpha=0.7)
        ax1.text(info['start']+1, ax1.get_ylim()[1]*0.945, season,
                 verticalalignment='top', fontsize=9, color='gray')
    
    ax1.set_xlim(0, 365)
    ax1.set_ylim(0, 1)
    ax1.set_ylabel('Percentage of 365-day Growth_cumulative Growth (%)')
    ax1.set_title(f'Normalized (Growth_cumulative) Stem Diameter [{agg_to_plot} aggregation]')
    ax1.legend()
    
    # --- Bottom subplot: sap flow ---
    col = col_map['Sap Flow']
    subset = df_all_norm[df_all_norm['Sensor'] == 'Sap Flow']
    
    for period in sorted(subset['Period'].unique()):
        period_subset = subset[subset['Period'] == period]
        ax2.plot(
            period_subset['Day_of_Period'],
            period_subset[col],
            label=f'Sap Flow Period {period+1}',
            alpha=0.7 if period == 0 else 0.9,
            color=sensor_colors['Sap Flow'] if period == 0 else 'forestgreen',
            zorder=2 if period == 0 else 4
        )
    
    for season, info in seasons.items():
        ax2.axvspan(info['start'], info['end'], color=season_colors[info['type']], alpha=0.05, zorder=0)
        ax2.axvline(x=info['start'], color='gray', linestyle='--', alpha=0.7)
        ax2.text(info['start']+1, ax2.get_ylim()[1]*0.95, season, fontsize=9, color='gray')
    
    ax2.set_xlabel('Day of Period (0–364)')
    ax2.set_ylabel('Detrended Sap Flow (raw units)')
    ax2.set_title(f'Sap Flow [{agg_to_plot} aggregation]')
    ax2.legend()
    
    plt.tight_layout()
    plt.show()



def add_season_shading(seasons, season_colors, ax, season_text, ypos=0.99, ):
    """
    Add seasonal shading to a plot, clipped to current x-axis limits.
    """
    
    xlim = ax.get_xlim()  # current visible x range
    
    for season, info in seasons.items():
        # Clip start and end to current xlim
        start = max(info['start'], xlim[0])
        end   = min(info['end'], xlim[1])
        if start >= end:
            continue  # season outside visible area
        
        ax.axvspan(start, end, color=season_colors[info['type']], alpha=0.05, zorder=0)
        ax.axvline(start, color='gray', linestyle='--', alpha=0.7)
        if season_text: 
            ax.text(
                start + 1, ypos,
                season,
                transform=ax.get_xaxis_transform(),  # x in data, y in axes
                va='top', ha='left',
                fontsize=9, color='gray'
            )

import matplotlib.dates as mdates

def add_season_shading_dates(season_dates, ax, season_text=True, ypos=0.99, xlim = None):
    season_colors = {'wet': 'deepskyblue', 'dry': 'gold'}
    if xlim is None:
        xlim = ax.get_xlim()

    for season, info in season_dates.items():
        # Convert timestamps to matplotlib float format
        start = mdates.date2num(info['start'])
        end   = mdates.date2num(info['end'])

        # Clip to current x-axis limits
        start_clipped = max(start, xlim[0])
        end_clipped   = min(end, xlim[1])
        if start_clipped >= end_clipped:
            continue

        ax.axvspan(start_clipped, end_clipped, color=season_colors[info['type']], alpha=0.05, zorder=0)
        ax.axvline(start_clipped, color='gray', linestyle='--', alpha=0.7)

        if season_text:
            ax.text(
                start_clipped + 2, ypos,
                season,
                transform=ax.get_xaxis_transform(),
                va='top', ha='left',
                fontsize=9, color='gray'
            )


def plot_sensor(df, ax, sensor, ycol, ylabel, title, 
                sensor_colors, season_colors,
                seasons = {
                    'JJAS': {'start': 0, 'end': 121, 'type':'dry'},  # June–Sept
                    'ON':   {'start': 121, 'end': 183, 'type':'wet'}, # Oct–Nov
                    'DJF':  {'start': 183, 'end': 274, 'type':'dry'}, # Dec–Feb
                    'MAM':  {'start': 274, 'end': 364, 'type':'wet'}  # Mar–May
                },
                xlim_dates = None, season_shading = True, 
                season_text = True):
    subset = df[df['Sensor'] == sensor].dropna(axis = 1, how = 'all')

    for period in sorted(subset['Period'].unique()):
        period_subset = subset[subset['Period'] == period]
        alpha = 0.9 if period == 0 else 0.6
        color = sensor_colors[sensor] if period == 0 else ("forestgreen" if sensor == 'Sap Flow' else sensor_colors[sensor])

        if sensor == 'TOMST':
            flat_mask = period_subset['flat_flag'] == 1
            normal = period_subset.loc[~flat_mask]
            ax.plot(
                normal['Day_of_Period'], normal[ycol],
                color=color, alpha=alpha, label=f'{sensor} Period {period+1}', zorder=2 if period == 0 else 5
            )

            flat_points = period_subset.loc[flat_mask, ['Day_of_Period', ycol]].sort_values('Day_of_Period')
            if not flat_points.empty:
                # Group contiguous Day_of_Period stretches in the flat points themselves
                groups = (flat_points['Day_of_Period'].diff() > 1).cumsum()
                for _, seg in flat_points.groupby(groups):
                    ax.plot(seg['Day_of_Period'], seg[ycol], color='red', zorder=1)

        else:
            ax.plot(period_subset['Day_of_Period'], period_subset[ycol],
                    color=color, alpha=alpha, label=f'{sensor} Period {period+1}', zorder=2 if period == 0 else 5)

    if season_shading:
        add_season_shading(seasons, season_colors, ax, season_text = season_text)

    ax.set_ylabel(ylabel)
    # ax.set_title(title)
    ax.legend(loc = 'lower right',frameon=False)

    # Apply zoom if xlim_dates is given
    if xlim_dates is not None:
        start_date, end_date = xlim_dates
        mask = (subset['DateTime'] >= start_date) & (subset['DateTime'] <= end_date)
        if mask.any():
            min_day = subset.loc[mask, 'Day_of_Period'].min()
            max_day = subset.loc[mask, 'Day_of_Period'].max()
            ax.set_xlim(min_day, max_day)
    



def add_event_rects(axes, events, color="red", alpha=0.2):
    """
    Add rectangles to highlight events on multiple subplots.
    
    axes   : list of matplotlib axes
    events : list of (xstart, xend) tuples
    """
    for ax in axes:
        ylim = ax.get_ylim()
        for xstart, xend in events:
            rect = patches.Rectangle(
                (xstart, ylim[0]),
                xend - xstart,
                ylim[1] - ylim[0],
                linewidth=1.5,
                edgecolor=color,
                facecolor=color,
                alpha=alpha,
                zorder=1
            )
            ax.add_patch(rect)




def plot_event_timeline(ax, timeline_df, base_colors={0: [0.2, 0.6, 0.8], 1: [0.8, 0.4, 0.2]},
                        height=0.25, y_gap=0.05):
    """
    Plot event timeline with overlapping events stacked and colored by period.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to plot on.
    timeline_df : pd.DataFrame
        DataFrame containing events with columns:
        'Boundary' ('start'/'end'), 'Period', 'Day_of_Period', 'Event_Type'.
    base_colors : dict
        RGB color for each period (keys are period numbers).
    height : float
        Height of each rectangle.
    y_gap : float
        Vertical gap between stacked overlapping events.
    """
    
    y_base_period = {0: 0.5, 1: 0.1}  # starting vertical positions per period

    # extract start/end rows
    starts = timeline_df[timeline_df["Boundary"] == "start"].reset_index(drop=True)
    ends   = timeline_df[timeline_df["Boundary"] == "end"].reset_index(drop=True)

    # store events per period
    period_events = {}
    for period in [0, 1]:
        events_p = starts[starts.Period==period].reset_index(drop=True)
        ends_p   = ends[ends.Period==period].reset_index(drop=True)

        n_events = len(events_p)
        row_indices = np.zeros(n_events, dtype=int)
        last_end = []

        for i in range(n_events):
            for r, e_end in enumerate(last_end):
                if events_p.loc[i,'Day_of_Period'] > e_end:
                    row_indices[i] = r
                    last_end[r] = ends_p.loc[i,'Day_of_Period']
                    break
            else:
                row_indices[i] = len(last_end)
                last_end.append(ends_p.loc[i,'Day_of_Period'])

        period_events[period] = (events_p, ends_p, row_indices)

    # plotting rectangles
    for period, (starts_p, ends_p, rows) in period_events.items():
        base_color = np.array(base_colors[period])
        max_row = max(rows)+1 if len(rows) else 1

        for start_row, end_row, r in zip(starts_p.itertuples(), ends_p.itertuples(), rows):
            shade_factor = 0.4 + 0.6*(1 - r/(max_row))  # higher rows lighter
            color = np.clip(base_color * shade_factor, 0, 1)

            y = y_base_period[period] + r*(height + y_gap)

            rect = plt.Rectangle(
                (start_row.Day_of_Period, y),
                end_row.Day_of_Period - start_row.Day_of_Period,
                height,
                color=color,
                alpha=0.7,
                hatch='/////'
            )
            ax.add_patch(rect)

            # label event type
            ax.text(
                (start_row.Day_of_Period + end_row.Day_of_Period)/2,
                y + height/2,
                start_row.Event_Type,
                ha='center', va='center',
                fontsize=8, fontweight='bold', color='black'
            )

    # axis formatting
    ax.set_ylim(0, 1)
    ax.set_xlabel('Day of Period (0–364)')
    ax.set_ylabel('Events')

    # create legend for periods
    legend_patches = [patches.Patch(color=base_colors[0], label='Period 1'),
                      patches.Patch(color=base_colors[1], label='Period 2')]
    ax.legend(handles=legend_patches, loc='upper right')
