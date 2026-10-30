"""
pulse_peak_screening.py - Multi-criteria screening of continuous-flow pulse peaks
(多标准连续流脉冲峰筛选).

Extracted verbatim from the SC-MetaboFlow source code:
    SCMeTA-English-R/SCMeTA/sc-software/sc-software/processcell/processData.py
    (functions clean_data, segmented_normalization,
     find_peaks_and_valleys, find_valley_before_peak_rt; lines 20-123)

Screening criteria (multi-criteria):
1. Signal cleaning: negative chromatogram intensities are set to 0
   (clean_data).
2. Segmented normalization: the chromatogram is divided into rolling
   segments of `segment_length` minutes and each segment is divided by
   its own maximum, so pulses of different heights are comparable
   (segmented_normalization).
3. Local-extremum criterion: scipy.signal.argrelextrema with an absolute
   intensity threshold (`threshold`, default 0.1) and a neighbourhood
   window (`order`, default 1). A point is a pulse PEAK when it exceeds
   the threshold and all its neighbours; a point is a VALLEY (baseline
   position) when it is below the threshold, positive, and lower than
   all its neighbours (find_peaks_and_valleys).
4. Pairing criterion: every peak is paired with its nearest preceding
   valley, which supplies the local baseline for that pulse
   (find_valley_before_peak_rt).

Dependencies: numpy, scipy
"""

import numpy as np
from scipy.signal import argrelextrema


def clean_data(data) -> np.ndarray:
    """Chinese clean_data_with_zero: copy data and set negative values to 0."""
    copy_data = np.copy(data)
    copy_data[copy_data < 0] = 0
    return copy_data


def segmented_normalization(
    rt: np.ndarray,
    intensity: np.ndarray,
    segment_length: float = 5.0,
) -> np.ndarray:
    """Chinese segmented_normalization.

    Walks the retention-time array from the start; each segment spans
    [rt[start_idx], rt[start_idx] + segment_length); the segment is divided
    by its maximum (if > 0). Matches the Chinese implementation exactly.
    """
    rt = np.asarray(rt, dtype=np.float64)
    intensity = np.asarray(intensity, dtype=np.float64)
    normalized = np.zeros_like(intensity, dtype=float)

    if len(rt) == 0:
        return normalized

    start_idx = 0
    n = len(rt)
    while start_idx < n:
        start_time = rt[start_idx]
        end_time = start_time + segment_length
        end_idx = np.searchsorted(rt, end_time, side='left')
        if end_idx == start_idx:
            end_idx = start_idx + 1
        segment_vals = intensity[start_idx:end_idx]
        max_val = np.max(segment_vals)
        if max_val > 0:
            normalized[start_idx:end_idx] = segment_vals / max_val
        start_idx = end_idx

    return normalized


def find_peaks_and_valleys(data_x, data_y: np.ndarray, threshold: float, order: int) -> tuple:
    """Chinese find_peaks_and_valleys (scipy argrelextrema version).

    Peaks:  a > threshold AND a > all neighbors within `order`.
    Valleys: a < threshold AND a > 0 AND a < all neighbors within `order`.

    Args:
        data_x: horizontal axis (unused; kept for signature parity with the
            Chinese version).
        data_y: normalized intensity array (0..1 per segment).
        threshold: absolute threshold (Chinese default 0.1).
        order: argrelextrema window order (Chinese default 1).

    Returns:
        Tuple of (peak_index_list, valley_index_list).
    """
    data_y = np.asarray(data_y, dtype=np.float64)

    def greater(a, b):
        return (a > threshold) & (a > b)

    def less(a, b):
        return (a < threshold) & (a > 0) & (a < b)

    peak_idxs = argrelextrema(data_y, greater, order=order)[0]
    valley_idxs = argrelextrema(data_y, less, order=order)[0]

    return peak_idxs.tolist(), valley_idxs.tolist()


def find_valley_before_peak_rt(rt: np.ndarray, peaks: list, valleys: list) -> list:
    """Chinese find_valley_before_peak_rt.

    For each peak index, pick the nearest valley index that precedes it
    (np.searchsorted). If the peak precedes all valleys, the FIRST valley
    is used; if the peak follows all valleys, the LAST valley is used.

    Args:
        rt: retention time array.
        peaks: peak indices.
        valleys: valley indices.

    Returns:
        List of tuples: (peak_idx, rt_peak, chosen_valley_idx, rt_valley).
    """
    result = []
    if not peaks or not valleys:
        return result

    valleys_sorted = sorted(valleys)
    for peak_idx in peaks:
        pos = np.searchsorted(valleys_sorted, peak_idx)
        if pos == 0:
            chosen_valley = valleys_sorted[0]
        elif pos >= len(valleys_sorted):
            chosen_valley = valleys_sorted[-1]
        else:
            chosen_valley = valleys_sorted[pos - 1]
        rt_peak = rt[peak_idx]
        rt_valley = rt[chosen_valley]
        result.append((peak_idx, rt_peak, chosen_valley, rt_valley))
    return result
