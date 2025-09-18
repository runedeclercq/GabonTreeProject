from matplotlib.ticker import MultipleLocator
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import calendar
import matplotlib.cm as cm
import numpy as np
from importlib import reload

import pandas as pd
import numpy as np

def clean_sensor_data(dfs, names, short_gap_limit):
    """
    Full time series and interpolates time series sensor data with gap handling.

    Parameters:
        dfs (list of pd.DataFrame): List of sensor dataframes.
        names (list of str): Sensor names corresponding to each dataframe.
        resolution (list of int): Sampling resolution in minutes for each sensor.
        large_gap_threshold (int): Threshold in hours to truncate long gaps.
        short_gap_limit (dict): Max consecutive missing rows to interpolate per sensor.

    Returns:
        list of pd.DataFrame: Cleaned and interpolated sensor dataframes.
    """
    cleaned_dfs = []

    for df_orig, name in zip(dfs, names):
        df_proc = df_orig.sort_values('DateTime').reset_index(drop=True)

        # Identify time gaps
        # df_proc['time_diff'] = df_proc['DateTime'].diff().shift(-1)
        # df_proc['large_gap'] = (df_proc['time_diff'] > pd.Timedelta(hours=large_gap_threshold)).fillna(False)


        # Create full time index
        res = (df_proc.set_index('DateTime')
               .index.to_series().diff()
               .dt.total_seconds() / 60).dropna().round().astype(int).mode()[0]

        full_time_index = pd.date_range(start=df_proc['DateTime'].min(),
                                        end=df_proc['DateTime'].max(),
                                        freq=f'{res}min')

        df_proc = (df_proc.groupby('DateTime').mean().reset_index()
                   .set_index('DateTime')
                   .reindex(full_time_index)
                   .rename_axis('DateTime')
                   .reset_index())

        # Drop time_diff column
        # df_proc = df_proc.drop(columns=['time_diff'], errors='ignore')


        valid_cols = []
        for col in df_proc.columns:
            if pd.api.types.is_numeric_dtype(df_proc[col]):
                valid_cols.append(col)

        df_proc = df_proc.set_index('DateTime')
        df_proc[valid_cols] = df_proc[valid_cols].interpolate(method='time', limit=short_gap_limit[name])
        df_proc = df_proc.reset_index()

        cleaned_dfs.append(df_proc)

    return cleaned_dfs


def detect_large_gaps(df, val_col, time_column='DateTime', large_gap_threshold = 12):
    """
    Detects periods of large gaps in a DataFrame and returns a summary DataFrame.

    Parameters:
    - df (pd.DataFrame): Input DataFrame containing a boolean column indicating gaps.
    - gap_column (str): Name of the column indicating large gaps (default 'large_gap').
    - time_column (str): Name of the datetime column (default 'DateTime').

    Returns:
    - pd.DataFrame: Summary of gap periods with start, end, duration, and notes.
    """
    df['missing'] = df[val_col].isna()
    res = (df.set_index(time_column).index.to_series().diff().dt.total_seconds() / 60).dropna().round().astype(int).mode()[0]

    # Recalculate large gaps based on consecutive missing rows
    df['gap_block'] = (df['missing'] != df['missing'].shift()).cumsum()
    gap_sizes = df.groupby('gap_block')['missing'].sum()
    large_gap_blocks = gap_sizes[gap_sizes * res >= large_gap_threshold * 60].index

    df['large_gap'] = df['gap_block'].isin(large_gap_blocks).astype(bool)
    # Ensure gap column is boolean
    df['large_gap'] = df['large_gap'].astype(bool)

    # Identify start and end of gaps
    df['gap_start'] = df['large_gap'] & ~df['large_gap'].shift(1).fillna(False).infer_objects(copy=False)
    df['gap_end'] = ~df['large_gap'] & df['large_gap'].shift(1).fillna(False).infer_objects(copy=False)

    # Extract timestamps
    gap_starts = df.loc[df['gap_start'], time_column].reset_index(drop=True)
    gap_ends = df.loc[df['gap_end'], time_column].reset_index(drop=True)

    # Handle edge case: gap starts but never ends
    if len(gap_starts) > len(gap_ends):
        gap_ends = pd.concat([gap_ends, pd.Series([df[time_column].iloc[-1]])], ignore_index=True)

    # Create summary DataFrame
    gap_ranges = list(zip(gap_starts, gap_ends))
    gap_summary = pd.DataFrame(gap_ranges, columns=['Gap Start', 'Gap End'])

    # Add duration and notes
    gap_summary['Duration'] = gap_summary['Gap End'] - gap_summary['Gap Start']
    gap_summary['Notes'] = ''

    # Clean up temporary columns
    df.drop(columns=['gap_start', 'gap_end'], inplace=True)

    return gap_summary


def get_equatorial_season(date):
    m = date.month
    if m in [10, 11]:
        return 'ON'      # Rainy season
    elif m in [12, 1, 2]:
        return 'DJF'     # Minor dry season
    elif m in [3, 4, 5]:
        return 'MAM'     # Rainy season
    elif m in [6, 7, 8, 9]:
        return 'JJAS'    # Major dry season

def compute_mds_from_features(df: pd.DataFrame, value_col='smoothed_signal',
                              morning=(4,8), afternoon=(12,18)):
    df = df.copy()
    df['date'] = df.index.date
    df['hour'] = df.index.hour + df.index.minute/60
    mds_list = []

    for d, g in df.groupby('date'):
        D_max = g[(g['hour'] >= morning[0]) & (g['hour'] < morning[1])][value_col].max()
        D_min = g[(g['hour'] >= afternoon[0]) & (g['hour'] < afternoon[1])][value_col].min()
        mds_list.append({'date': d, 'MDS': D_max - D_min})

    return pd.DataFrame(mds_list).set_index('date')


def compute_dgr(df, value_col='smoothed_signal', morning=(4,8)):
    """
    Compute Daily Growth Rate (DGR) based on early-morning maximums.
    Returns a DataFrame with daily DGR.
    """
    df = df.copy()
    df['date'] = df.index.date
    df['hour'] = df.index.hour + df.index.minute / 60
    
    morning_max = df[(df['hour'] >= morning[0]) & (df['hour'] < morning[1])] \
                    .groupby('date')[value_col].max()
    
    dgr = morning_max.diff()  # day-to-day difference
    dgr.name = 'DGR'
    return dgr

def compute_metrics(df, value_col='detrended_daily_mean', morning=(0,8), afternoon=(12,18)):
    """
    Compute MDS, DGR, TWD, TWDnorm, and MDSnorm from dendrometer/sap flow data.
    
    Parameters
    ----------
    df : pd.DataFrame
        Time-indexed dataframe with sensor measurements.
    value_col : str
        Column with sensor values.
    morning : tuple
        Hours to define morning max/min (for MDS and TWD).
    afternoon : tuple
        Hours to define afternoon min (for MDS).
    
    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - MDS
        - DGR
        - TWD (pre-dawn)
        - TWDnorm
        - MDSnorm
    """
    df = df.copy()
    df['date'] = df.index.date
    df['hour'] = df.index.hour + df.index.minute / 60
    
    # --- Compute MDS ---
    mds_list = []
    for d, g in df.groupby('date'):
        D_max = g[(g['hour'] >= morning[0]) & (g['hour'] < morning[1])][value_col].max()
        D_min = g[(g['hour'] >= afternoon[0]) & (g['hour'] < afternoon[1])][value_col].min()
        mds_list.append({'date': d, 'MDS': D_max - D_min})
    mds_df = pd.DataFrame(mds_list).set_index('date')
    
    # --- Compute DGR ---
    morning_max = df[(df['hour'] >= morning[0]) & (df['hour'] < morning[1])] \
                    .groupby('date')[value_col].max()
    dgr = morning_max.diff()
    dgr.name = 'DGR'
    
    # --- Compute TWD ---
    twd = (df[value_col].cummax() - df[value_col]).clip(lower=0)
    twd.name = 'TWD'

    df['TWD'] = twd

    # --- Compute pre-dawn TWD ---
    pre_dawn = df[(df['hour'] >= morning[0]) & (df['hour'] < morning[1])]
    pre_dawn_twd = pre_dawn.groupby('date')['TWD'].min()
    pre_dawn_twd.name = 'pre_dawn_TWD'

    # --- Compute normalized metrics ---
    MDS_max = mds_df['MDS'].quantile(0.99)  # 99th percentile over multi-year
    twd_norm = pre_dawn_twd / MDS_max
    twd_norm.name = 'TWDnorm'
    
    mds_norm = mds_df['MDS'] / MDS_max
    mds_norm.name = 'MDSnorm'
    
    filtered_dgr = dgr.where((twd_norm < 0.6) & (mds_norm < 0.5))
    filtered_dgr.name = 'filtered_DGR'

    # --- Compute daily amplitude ---
    daily_max = df.groupby('date')[value_col].max()
    daily_min = df.groupby('date')[value_col].min()
    daily_amp = (daily_max - daily_min).rename('daily_amplitude')

    # --- Combine all metrics ---
    daily_metrics = pd.concat([
        mds_df, 
        dgr, 
        pre_dawn_twd, 
        twd_norm, 
        mds_norm, 
        filtered_dgr, 
        daily_amp
    ], axis=1)
    
    time_metrics = df[['TWD']].copy()


    daily_metrics = daily_metrics.reset_index()
    daily_metrics['date'] = pd.to_datetime(daily_metrics['date']).dt.date

    return daily_metrics, time_metrics


def smooth_function(df: pd.DataFrame, value_col: str, 
                         short_window = '1h', 
                         baseline_window = '3D', 
                         long_trend_window = '15D',
                         amplitude_window = '24h') -> pd.DataFrame:
    """
    Extract features with separate short-term and long-term detrended signals.
    
    - smoothed_signal: rolling mean of 'short_window' -> 1h
    - detrended_daily_mean: for daily changes (smoothed window minus daily mean)
    - long_term_detrended: captures multi-day changes, preserves flushes
    - rolling median: smooths short-term noise
    - rolling amplitude: approximate diurnal variation
    - day-to-day diff of baseline
    - daily summaries (min/max/amp)

    """
    out = df.copy()

    out.set_index('DateTime', inplace=True)
    
    raw = out[value_col]
    # compute sensor resolution as timedifference between 2 following timestamps (in minutes)
    res = (out.index.to_series().diff().dt.total_seconds() / 60).dropna().round().astype(int).mode()[0]
    
    # windows to number of points
    short_points = int(pd.Timedelta(short_window).total_seconds() / 60 / res)
    baseline_points = int(pd.Timedelta(baseline_window).total_seconds() / 60 / res)
    long_trend_points = int(pd.Timedelta(long_trend_window).total_seconds() / 60 / res)
    amplitude_points = int(pd.Timedelta(amplitude_window).total_seconds() / 60 / res)
   
    # rolling median to smooth high-frequency noise
    out["smoothed_signal"] = raw.rolling(short_points, center=True, min_periods=1).median()
    
    v = out["smoothed_signal"] # use smoothed signal
    
    # # daily mean subtraction for diurnal signal
    daily_mean_smoothed = v.resample('D').mean()
    out['daily_mean'] = daily_mean_smoothed

    out["detrended_daily_mean"] = v - daily_mean_smoothed.reindex(out.index, method='ffill')

    ## --- rolling mean for trend extraction
    # 1 short trend
    out["short_trend"] = v.rolling(baseline_points, center=True, min_periods=baseline_points//2).mean()
    out["detrended_short"] = v - out["short_trend"]
    
    # 2 long-term trend (preserves multi-day flushes)
    out["long_trend"] = v.rolling(long_trend_points, center=True, min_periods=long_trend_points//2).mean()
    out["detrended_long"] = v - out["long_trend"]

    # Day-to-day diff of short trend
    daily_baseline = out["short_trend"].resample('D').mean()
    out["short_trend_dday_diff"] = daily_baseline.diff().reindex(out.index, method='ffill')
    
    # daily difference (daily growth)
    # out["daily_mean_dday_diff"] = out["daily_mean"].diff().reindex(out.index, method='ffill')   

    # Daily min/max/amp (summary per day)
    # daily_change = out["detrended_daily_mean"]
    # out["daily_min"] = daily_change.min().reindex(out.index, method='ffill')
    # out["daily_max"] = daily_change.max().reindex(out.index, method='ffill')
    # out["daily_amp"] = out["daily_max"] - out["daily_min"]
    
    return out

def extract_waveform_metrics(df: pd.DataFrame, signal_col: str = "detrended_daily_mean") -> pd.DataFrame:

    """
    Extract daily waveform metrics from a time-series signal.

    Parameters:
        df : pd.DataFrame with datetime index
        signal_col : str, name of the column containing the raw signal

    Returns:
        pd.DataFrame with DateTime index and columns:
            - Rmax: daily maximum
            - Rmin: daily minimum
            - Tmax: time of day when Rmax occurs
            - Tmin: time of day when Rmin occurs
            - DeltaR: amplitude (Rmax - Rmin)
            - R_mean: mean of detrended signal per day
    """
    rmax_list, rmin_list, tmax_list, tmin_list, amp_list, rmean_list = [], [], [], [], [], []
    dates = []

    # Group by calendar day
    for day, group in df.groupby(df.index.date):
        signal = group[signal_col].dropna()
        if signal.empty or signal.is_monotonic_increasing or signal.is_monotonic_decreasing:
            continue  # Skip days without clear waveform

        rmax = signal.max()
        rmin = signal.min()
        tmax = signal.idxmax().time()
        tmin = signal.idxmin().time()
        amp = rmax - rmin
        rmean = signal.mean()

        rmax_list.append(rmax)
        rmin_list.append(rmin)
        tmax_list.append(tmax)
        tmin_list.append(tmin)
        amp_list.append(amp)
        rmean_list.append(rmean)
        dates.append(pd.Timestamp(day))

    waveform_df = pd.DataFrame({
        'DateTime': dates,
        'Rmax': rmax_list,
        'Rmin': rmin_list,
        'Tmax': tmax_list,
        'Tmin': tmin_list,
        'DeltaR': amp_list,
        'R_mean': rmean_list
    }).set_index('DateTime')

    return waveform_df



def summarize_multi_levels(df, value_cols, agg_levels=['D','W','M'], max_missing_pct=0.2, rolling_window = None):

    df = df.set_index('DateTime').sort_index()

    # Add rolling columns if requested
    if rolling_window is not None:
        res = (df.index.to_series().diff().dt.total_seconds() / 60).dropna().round().astype(int).mode()[0]
        rolling_days = int(pd.Timedelta(rolling_window).total_seconds() / 60 / res)

        for col in value_cols:
            df[f"{col}_roll{rolling_window}"] = df[col].rolling(
                rolling_days, center=True, min_periods=1).mean()

    result = df.copy()

    full_index = df.index  # preserve all 20-min timestamps
            
    for level in agg_levels:
        period = level
        if level=='M':
            period='MS'

        # Aggregate
        agg = df[value_cols].resample(period, label='left').agg(['mean','min','max'])

        # Flatten MultiIndex
        agg.columns = ['_'.join(col).strip() for col in agg.columns.values]

        # Mask missing data
        for col in value_cols:
            counts = df[col].resample(period).count()
            total = df[col].resample(period).size()
            mask = counts >= total*(1-max_missing_pct)
            for stat in ['mean','min','max']:
                col_name = f"{col}_{stat}"
                if col_name in agg.columns:
                    agg[col_name] = agg[col_name].where(mask, np.nan)

        # Forward-fill to 20-min timestamps
        agg_flat = agg.reindex(full_index, method='ffill')
        
        # Add suffix for aggregation level
        agg_flat = agg_flat.add_suffix(f'_{level}')

        # Merge
        result = result.merge(agg_flat, left_index=True, right_index=True)

    result = result.reset_index()
    return result



def add_period_normalization(df, rescale = False, value_cols = ['Growth_cumulative'], start_time=None, period_days=365):
    """
    Normalize multiple value columns per-period and rescale later periods so that
    the end-of-period percent matches the corresponding day in period 0.

    - df : DataFrame with 'DateTime' column (datetime64)
    - value_cols : list of column names (e.g. ['Growth_cumulative_mean_D','Growth_mean_D'])
    - start_time : optional pd.Timestamp (defaults to df['DateTime'].min())
    - period_days : number of days per period (default 365)

    Returns a copy of df with added columns for each value_col:
      - <col>_pct                 : normalized 0-100% per Period
      - <col>_pct_rescaled        : same as _pct but rescaled for Periods > 0
      plus 'Days_since_start','Period','Day_of_Period'.
    """
    df = df.copy()
    if start_time is None:
        start_time = df['DateTime'].min()

    # period indexing
    df['Days_since_start'] = (df['DateTime'] - start_time).dt.total_seconds() / (24 * 3600)
    df['Period'] = (df['Days_since_start'] // period_days).astype(int)
    df['Day_of_Period'] = df['Days_since_start'] % period_days

    # ensure data sorted by DateTime for asof / reindex operations
    df = df.sort_values('DateTime').reset_index(drop=True)

    if rescale:
        for col in value_cols:
            # compute per-period min/max
            # Step 1: global scaling 0–100%
            overall_min = df[col].min()
            overall_max = df[col].max()
            overall_range = overall_max - overall_min

            if overall_range != 0:
                df[f'{col}_pct_global'] = (df[col] - overall_min) / overall_range
            else:
                df[f'{col}_pct_global'] = 0.0

            # Step 2: per-period relative growth for plotting
            df[f'{col}_pct_period'] = np.nan  # initialize

            for p in df['Period'].unique():
                mask = df['Period'] == p
                period_min = df.loc[mask, f'{col}_pct_global'].min()
                period_max = df.loc[mask, f'{col}_pct_global'].max()
                period_range = period_max - period_min

                # scale within period from 0 → period's contribution to total growth
                df.loc[mask, f'{col}_pct_period'] = ((df.loc[mask, f'{col}_pct_global'] - period_min) 
                                                    / period_range * period_range if period_range != 0 else 0)

    return df
