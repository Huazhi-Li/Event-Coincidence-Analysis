# -*- coding: utf-8 -*-
"""
Function to convert discrete event dates into a binary time series, and function to conver time series to event series.

@author: Huazhi Li (huazhi.li@vu.nl)
"""
import numpy as np
import pandas as pd
from datetime import timedelta

def ts(date_range, event_dates, locations, pad=182, freq="weekly"):
    """
    date_range:     the start and end dates of the time period.
    event_dates:    start dates of events.
    buffer:         buffering days added to the start and end of the time series, defaulft set to 182 days (half a year).
    freq:           frequency of the time series, now only weekly data is supported.
    """
    
    if freq != "weekly":
        raise NotImplementedError("Only weekly frequency supported") # this can be further developed by considering data with different frequencies 
    
    start = pd.to_datetime(date_range[0]).tz_localize(None).normalize() - timedelta(days=pad)
    end   = pd.to_datetime(date_range[1]).tz_localize(None).normalize() + timedelta(days=pad)

    dates = pd.date_range(start=start, end=end, freq="7D")

    # build dataframe for events
    events_df = pd.DataFrame({
        "date": pd.Series(event_dates),
        "location": locations
    })
    
    events_df["date"] = (
    pd.to_datetime(events_df["date"], errors="coerce")
      .dt.tz_convert(None)
      .dt.normalize()
    )

    # assign week bins
    bins = list(dates) + [dates[-1] + pd.Timedelta(days=7)]
    events_df["week_bin"] = pd.cut(
        events_df["date"],
        bins=bins,
        right=False,
        include_lowest=True
    )
    
    # binary series
    # binary_series = (
    #     events_df["week_bin"]
    #     .value_counts(sort=False)
    #     .reindex(pd.IntervalIndex(events_df["week_bin"].cat.categories))
    #     .fillna(0)
    #     .astype(int)
    # )
    all_bins = pd.IntervalIndex.from_breaks(bins, closed="left")

    counts = (
        events_df["week_bin"]
        .value_counts(sort=False)
        .reindex(all_bins)
        .fillna(0)
        .astype(int)
    )
    
    binary_series = (counts > 0).astype(int)
    
    # aggregated locations per week
    weekly_locations = (
        events_df
        .groupby("week_bin",observed=False)["location"]
        .apply(list)
        .reindex(pd.IntervalIndex(events_df["week_bin"].cat.categories))
        .apply(lambda x: x if isinstance(x, list) else [])
        .tolist()
    )
    
    return binary_series, weekly_locations

def ts2es(time_sries):
    es = np.where(time_sries == 1)[0]
    span = (0, len(time_sries)-1)
    return es, span
