"""
isotope_filtering.py - Weighted C13 isotope filtering algorithm
(加权13C同位素过滤算法).

Extracted verbatim from the SC-MetaboFlow source code:
    SCMeTA-English-R/SCMeTA/sc-software/sc-software/processcell/processData.py
    (functions estimate_carbon_count_fast, mark_c13_isotopes_ahead,
     advanced_exclude_isotopes, is_isotope, exclude_isotopes; lines 126-277)

Algorithm
---------
For each centroid (sorted by m/z), the filter searches up to 20 positions
back for a candidate C13 mother ion at m/z - ~1.00335 (C13 mass shift).
A centroid is dropped as a C13 isotopic false peak when a
high-confidence mother is found. The confidence is a WEIGHTED score:

    overall_confidence = mass_confidence * 0.5
                       + intensity_confidence * 0.4
                       + ratio_reasonableness * 0.1

where mass_confidence = 1 - ppm_error / mass_tolerance_ppm,
intensity_confidence = 1 - |actual - theoretical| / theoretical ratio
(theoretical ratio estimated from the carbon count of the mother,
C13 natural abundance 1.1%), and ratio_reasonableness penalizes ratios
that deviate more than 3x / 0.3x from the theoretical value.
Centroids above min_isotope_probability (default 0.6) are excluded;
kept centroids pre-mark their own C13 candidates ahead
(mark_c13_isotopes_ahead) so they are skipped later.

A simpler adjacent-pair variant (exclude_isotopes / is_isotope), used in
the Chinese original when only neighbouring peaks are checked, is also
included.

Dependencies: numpy only
"""

import numpy as np


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
