"""
run_example.py - End-to-end example of the SC-MetaboFlow core algorithms.

Runs the three algorithm modules on the small example files in
example_data/ (no Thermo RAW files or GUI needed):

  Step 1  pulse_peak_screening.py : multi-criteria continuous-flow pulse
                                    peak screening (clean -> segmented
                                    normalization -> peak/valley detection
                                    -> peak-valley pairing)
  Step 2  baseline_subtraction.py : paired peak-valley local baseline
                                    subtraction on the first peak detected
  Step 3  isotope_filtering.py    : weighted C13 isotope filtering on a
                                    centroid list

Outputs are written to example_data/output/.

Requirements: numpy, pandas, scipy (see dependencies/requirements_v2.txt)
Usage:        python run_example.py
"""

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)  # make the three algorithm modules importable

from pulse_peak_screening import (  # noqa: E402
    clean_data, segmented_normalization,
    find_peaks_and_valleys, find_valley_before_peak_rt,
)
from baseline_subtraction import paired_peak_valley_baseline_subtraction  # noqa: E402
from isotope_filtering import advanced_exclude_isotopes  # noqa: E402

DATA = os.path.join(HERE, "example_data")
OUT = os.path.join(DATA, "output")
os.makedirs(OUT, exist_ok=True)


def banner(title):
    line = "=" * 64
    print(f"\n{line}\n{title}\n{line}")


# ---------------------------------------------------------------
# Step 1: pulse peak screening
# ---------------------------------------------------------------
banner("Step 1  Pulse peak screening (pulse_peak_screening.py)")

chrom = pd.read_csv(os.path.join(DATA, "example_chromatogram.csv"))
rt = chrom["rt_min"].to_numpy()
y_raw = chrom["intensity"].to_numpy()

y_clean = clean_data(y_raw)                          # clip negative values to 0
y_norm = segmented_normalization(rt, y_clean, segment_length=1.0)  # per-minute segments
peaks, valleys = find_peaks_and_valleys(rt, y_norm, threshold=0.1, order=1)
pairs = find_valley_before_peak_rt(rt, peaks, valleys)

n_clipped = int((y_raw < 0).sum())
print(f"Input            : {len(rt)} scans over {rt[0]:.2f}-{rt[-1]:.2f} min "
      f"(6 cells, one pulse per minute)")
print(f"clean_data       : {n_clipped} negative values clipped to 0")
print(f"peaks detected   : {len(peaks)}  ->  rt = {np.round(rt[peaks], 3).tolist()}")
print(f"valleys detected : {len(valleys)}")
print(f"peak-valley pairs: {len(pairs)}  ->  {pairs}")

flags = np.array([""] * len(rt), dtype=object)
for i in peaks:
    flags[i] = "peak"
for i in valleys:
    flags[i] = "valley"
out_chrom = pd.DataFrame({
    "rt_min": rt, "raw": y_raw, "cleaned": y_clean,
    "normalized": np.round(y_norm, 6), "mark": flags,
})
out_chrom.to_csv(os.path.join(OUT, "example_chromatogram_processed.csv"), index=False)
print(f"Saved: example_data/output/example_chromatogram_processed.csv")

# ---------------------------------------------------------------
# Step 2: paired peak-valley baseline subtraction
# ---------------------------------------------------------------
banner("Step 2  Paired peak-valley baseline subtraction (baseline_subtraction.py)")

if len(pairs) == 0:
    print("No peak-valley pair found - cannot demonstrate baseline subtraction.")
else:
    # Build MS1 spectra on a shared m/z grid for the first peak (apex)
    # and its preceding valley, as the real pipeline does.
    mz_grid = np.round(np.arange(99.5, 103.51, 0.05), 2)
    baseline = 300.0 + 50.0 * np.sin(mz_grid * 4.0)
    analyte = (800.0 * (np.abs(mz_grid - 100.0) < 0.026)
               + 500.0 * (np.abs(mz_grid - 100.6) < 0.026)
               + 300.0 * (np.abs(mz_grid - 101.2) < 0.026)
               + 200.0 * (np.abs(mz_grid - 102.8) < 0.026))
    i_apex = baseline + analyte
    i_valley = baseline.copy()

    mz_sub, i_sub = paired_peak_valley_baseline_subtraction(
        mz_grid, i_apex, mz_grid, i_valley)

    print(f"First peak-valley pair used: peak #{pairs[0][0]} at rt "
          f"{pairs[0][1]:.2f} <- valley #{pairs[0][2]} at rt {pairs[0][3]:.2f}")
    print(f"Total intensity - apex   : {i_apex.sum():10.1f}")
    print(f"Total intensity - valley : {i_valley.sum():10.1f}")
    print(f"Total intensity - after  : {i_sub.sum():10.1f}   "
          f"(background removed, analyte peaks kept)")
    kept = mz_grid[i_sub > 0.1 * np.max(i_sub)]
    print(f"Analyte peaks retained at m/z: {kept.tolist()}")
    pd.DataFrame({"mz": mz_sub, "apex": i_apex, "valley": i_valley,
                  "subtracted": i_sub}).to_csv(
        os.path.join(OUT, "example_baseline_subtraction.csv"), index=False)
    print(f"Saved: example_data/output/example_baseline_subtraction.csv")

# ---------------------------------------------------------------
# Step 3: weighted C13 isotope filtering
# ---------------------------------------------------------------
banner("Step 3  Weighted C13 isotope filtering (isotope_filtering.py)")

cent_df = pd.read_csv(os.path.join(DATA, "example_centroids.csv"))
centroids = list(zip(cent_df["mz"], cent_df["intensity"]))
filtered = advanced_exclude_isotopes(centroids, mass_tolerance_ppm=5.0)

removed = [c for c in centroids if c not in filtered]
print(f"Input centroids    : {[(round(m, 4), int(i)) for m, i in centroids]}")
print(f"Removed as C13     : {[(round(m, 4), int(i)) for m, i in removed]}"
      f"  (mother ion m/z 300.0, C13 child at +1.00335)")
print(f"Kept centroids     : {[(round(m, 4), int(i)) for m, i in filtered]}")
pd.DataFrame(filtered, columns=["mz", "intensity"]).to_csv(
    os.path.join(OUT, "example_centroids_filtered.csv"), index=False)
print(f"Saved: example_data/output/example_centroids_filtered.csv")

banner("Example completed - see example_data/output/ for all result files")
