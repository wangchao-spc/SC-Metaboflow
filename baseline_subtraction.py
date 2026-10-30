"""
baseline_subtraction.py - Paired peak-valley local baseline subtraction
(配对峰谷局部基线扣除算法).

Extracted from the SC-MetaboFlow source code:
- English app: SCMeTA-English-R/sc_metaboflow_app_v2.py, PreprocessingWorker
  background-subtraction block (lines 878-894)
- Chinese original: SCMeTA/sc-software/sc-software/processcell/main7.py
  (lines 73-127)

Algorithm
---------
For each chromatographic pulse peak, the MS1 spectrum at its preceding
valley is taken as the local baseline. Because the valley spectrum and the
peak spectrum have different numbers of data points, the valley spectrum
is linearly interpolated onto the peak spectrum's m/z grid
(scipy.interpolate.interp1d, kind='linear'; outside the valley's m/z
range the baseline is taken as 0). The baseline-subtracted spectrum is:

    intensity_final = max(peak_intensity - interpolated_baseline, 0)

i.e. negative residuals are clipped to zero; m/z positions are unchanged.
(The Chinese original uses kind='slinear' and clips the m/z grid to the
valley's range instead - same algorithm, see main7.py lines 85-121.)
In the full pipeline, if interpolation raises an exception the raw peak
spectrum is used as fallback (sc_metaboflow_app_v2.py lines 891-892).

Dependencies: numpy, scipy
"""

import numpy as np
from scipy.interpolate import interp1d


def paired_peak_valley_baseline_subtraction(
    mz_peak, intensity_peak, mz_valley, intensity_valley
):
    """Subtract the local baseline (valley spectrum) from a peak spectrum.

    Args:
        mz_peak: m/z array of the peak-apex MS1 spectrum.
        intensity_peak: intensity array of the peak-apex MS1 spectrum.
        mz_valley: m/z array of the preceding valley's MS1 spectrum.
        intensity_valley: intensity array of the valley's MS1 spectrum.

    Returns:
        (mz_final, intensity_final): m/z grid is unchanged (= mz_peak);
        intensity_final = max(peak - interpolated_baseline, 0).
    """
    mz_p = np.asarray(mz_peak, dtype=np.float64)
    i_p = np.asarray(intensity_peak, dtype=np.float64)
    mz_v = np.asarray(mz_valley, dtype=np.float64)
    i_v = np.asarray(intensity_valley, dtype=np.float64)

    if len(mz_p) == 0 or len(mz_v) == 0:
        return mz_p, i_p

    # Linear interpolation of the valley (baseline) spectrum onto the
    # peak spectrum's m/z grid (sc_metaboflow_app_v2.py lines 883-886)
    f_background = interp1d(
        mz_v, i_v, kind='linear', fill_value=0, bounds_error=False
    )
    background_at_peak_positions = f_background(mz_p)

    # net spectrum = peak - baseline, negatives clipped to 0
    # (sc_metaboflow_app_v2.py lines 887-888)
    intensity_final = np.maximum(i_p - background_at_peak_positions, 0)
    return mz_p, intensity_final
