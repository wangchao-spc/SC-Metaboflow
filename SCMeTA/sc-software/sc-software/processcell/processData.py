"""
processData.py - Signal processing functions for chromatogram data analysis.

Fully aligned with the Chinese version (SC-MetaboFlow) processData module.
Logic ported from the Chinese build's bytecode (verified via disassembly,
2026-09-29): PeakProcessor.process_single + processData functions.

Provides functions for:
- Data cleaning (negative values -> 0, Chinese clean_data_with_zero)
- Segmented normalization (rolling searchsorted segments, Chinese version)
- Peak and valley detection (scipy argrelextrema, Chinese version)
- Peak-valley matching (Chinese find_valley_before_peak_rt)
- Advanced C13 isotope exclusion (Chinese advanced_exclude_isotopes)
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


def estimate_carbon_count_fast(mz: float) -> int:
    """Chinese estimate_carbon_count_fast."""
    if mz < 100:
        carbon_count = max(1, int(mz / 15))
    elif mz < 300:
        carbon_count = max(1, int(mz / 18))
    else:
        carbon_count = max(1, int(mz / 22))
    return min(carbon_count, 30)


def mark_c13_isotopes_ahead(
    base_index: int,
    centroids: list,
    skip_indices: set,
    mass_tolerance_ppm: float,
    intensity_tolerance: float,
    min_isotope_probability: float,
):
    """Chinese mark_c13_isotopes_ahead.

    For a kept centroid at base_index, pre-mark C13 isotope candidates at
    higher m/z (up to 20 positions ahead) so they are skipped later.
    """
    mz_base, intensity_base = centroids[base_index]
    search_end = min(len(centroids), base_index + 20)
    for k in range(base_index + 1, search_end):
        if k in skip_indices:
            continue
        mz_candidate, intensity_candidate = centroids[k]
        mass_diff = mz_candidate - mz_base
        if mass_diff < 0.95 or mass_diff > 1.05:
            if mass_diff > 2.0:
                return  # centroids are m/z-sorted: no more +1 candidates
            continue
        ppm_error = abs(mass_diff - 1.00335) / 1.00335 * 1000000.0
        if ppm_error > mass_tolerance_ppm:
            continue
        estimated_carbons = estimate_carbon_count_fast(mz_base)
        theoretical_ratio = min(estimated_carbons * 0.011, 0.1)
        actual_ratio = intensity_candidate / intensity_base if intensity_base > 0 else 0
        if theoretical_ratio > 0:
            ratio_error = abs(actual_ratio - theoretical_ratio) / theoretical_ratio
        else:
            ratio_error = 1.0
        if ratio_error > intensity_tolerance:
            continue
        confidence = (1.0 - ppm_error / mass_tolerance_ppm) * (1.0 - ratio_error)
        if confidence < min_isotope_probability:
            continue
        skip_indices.add(k)


def advanced_exclude_isotopes(
    centroids: list,
    mass_tolerance_ppm: float = 5.0,
    intensity_tolerance: float = 0.3,
    min_isotope_probability: float = 0.6,
) -> list:
    """Chinese advanced_exclude_isotopes.

    centroids must be sorted by m/z (PeakProcessor sorts by lambda x: x[0]).
    A centroid is dropped when a high-confidence C13 mother (m/z - ~1.00335,
    up to 20 positions back) is found; kept centroids pre-mark their own
    C13 candidates ahead. Returns the filtered centroid list.
    """
    if len(centroids) < 2:
        return centroids

    C13_MASS_DIFF = 1.00335
    C13_ABUNDANCE = 0.011

    filtered_centroids = []
    skip_indices = set()

    for i, (mz_base, intensity_base) in enumerate(centroids):
        if i in skip_indices:
            continue
        is_c13_isotope = False
        search_start = max(0, i - 20)
        for j in range(search_start, i):
            if j in skip_indices:
                continue
            mz_mother, intensity_mother = centroids[j]
            actual_mass_diff = mz_base - mz_mother
            if actual_mass_diff < 0.95 or actual_mass_diff > 1.05:
                continue
            ppm_error = abs(actual_mass_diff - C13_MASS_DIFF) / C13_MASS_DIFF * 1000000.0
            if ppm_error > mass_tolerance_ppm:
                continue
            estimated_carbons = estimate_carbon_count_fast(mz_mother)
            theoretical_ratio = estimated_carbons * C13_ABUNDANCE
            if theoretical_ratio > 0.5:
                theoretical_ratio = 0.1
            actual_ratio = intensity_base / intensity_mother if intensity_mother > 0 else 0
            if theoretical_ratio > 0:
                ratio_error = abs(actual_ratio - theoretical_ratio) / theoretical_ratio
            else:
                ratio_error = 1.0
            if ratio_error > intensity_tolerance:
                continue
            mass_confidence = 1.0 - ppm_error / mass_tolerance_ppm
            intensity_confidence = 1.0 - ratio_error
            ratio_reasonableness = 1.0
            if actual_ratio > theoretical_ratio * 3 or actual_ratio < theoretical_ratio * 0.3:
                ratio_reasonableness = 0.5
            overall_confidence = (
                mass_confidence * 0.5 + intensity_confidence * 0.4
                + ratio_reasonableness * 0.1
            )
            if overall_confidence < min_isotope_probability:
                continue
            is_c13_isotope = True
            break
        if is_c13_isotope:
            continue
        filtered_centroids.append((mz_base, intensity_base))
        mark_c13_isotopes_ahead(
            i, centroids, skip_indices,
            mass_tolerance_ppm, intensity_tolerance, min_isotope_probability,
        )
    return filtered_centroids


def is_isotope(peak_mz: float, next_peak_mz: float, threshold: float = 1.00335) -> bool:
    """Chinese is_isotope: |(next_peak_mz - peak_mz) - threshold| < 0.001."""
    return abs(next_peak_mz - peak_mz - threshold) < 0.001


def exclude_isotopes(centroids: list) -> list:
    """Chinese exclude_isotopes.

    centroids: list of (mz, intensity) tuples, sorted by m/z ascending.
    Adjacent isotope pairs: keep the first of each pair, skip the second.

    Returns:
        Filtered centroid list.
    """
    non_iso = []
    i = 0
    n = len(centroids)
    while i < n:
        mz_i, inten_i = centroids[i]
        if i < n - 1:
            mz_next, inten_next = centroids[i + 1]
            if is_isotope(mz_i, mz_next):
                non_iso.append((mz_i, inten_i))
                i += 2
                continue
        non_iso.append((mz_i, inten_i))
        i += 1
    return non_iso
