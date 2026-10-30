# SC-MetaboFlow — Repository Documentation

## Software Overview

SC-MetaboFlow is an automated data-processing software for continuous-flow pulse data from single-cell mass spectrometry, integrating core algorithms including multi-criteria pulse peak screening, paired peak-valley local baseline subtraction, and weighted ¹³C isotope filtering. The software links preprocessing, data cleaning, group splitting, and multivariate statistical visualization into a single graphical workflow, and provides quantitative performance benchmarks across different algorithm configurations (Basic / Background subtract / Isotope filter / Full optimized).

---

## How to Install Dependencies

### Python Environment (Required)

**Python version requirement: ≥ 3.12** (the development and testing environment was Python 3.14.6).

Recommended installation using a virtual environment (Windows PowerShell / CMD):

```bat
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r dependencies\requirements.txt
```

Alternatively, run the one-click setup script (creates the virtual environment and installs all dependencies automatically):

```bat
dependencies\install.bat
```

> `requirements.txt` is the fully pinned version (pip freeze; the exact versions used to validate the software, including PyQt6, pandas, numpy, scipy, scikit-learn, matplotlib, fisher_py, pycombat, umap-learn, etc.).
> If the pinned versions fail to install on your platform, use the loose versions instead: `pip install -r dependencies\requirements_v2.txt`.

### R Environment (Optional; only required for the R Analysis tab)

1. Install R from https://cran.r-project.org/ and add it to the system PATH;
2. The R package list is provided in `dependencies\requirements_r.txt`. On first use of the R Analysis tab, the software installs the following packages automatically:
   `Seurat`, `monocle`, `EigenMS`, `ggplot2`, `dplyr`, `Matrix`.

### Regarding environment.yml / packages.R

This project does not use Conda; dependencies are managed entirely with pip + venv (`requirements.txt`), so no `environment.yml` is needed. The R package list is already included in `requirements_r.txt`, so a separate `packages.R` is not needed.

---

## How to Run Example Data

### 1. Launch the GUI (main program)

```bat
python sc_metaboflow_main.py
```

Runtime requirements: the `SCMeTA` compute backend is bundled in this repository (`SCMeTA\`, code-only copy ~1.4 MB; the original data folders are excluded) and is located automatically at startup. Reading Thermo RAW data (fisher_py) on Windows requires .NET Framework 4.8.

### 2. Command-Line Usage of the Core Algorithms (Examples)

The three algorithm modules can be used independently of the GUI (numpy / scipy only):

**① Multi-criteria continuous-flow pulse peak screening** (`pulse_peak_screening.py`):

```bat
python -c "import numpy as np; from pulse_peak_screening import clean_data, segmented_normalization, find_peaks_and_valleys, find_valley_before_peak_rt; rt=np.linspace(0,10,101); y=0.01+0.003*np.sin(rt); y[40]=1.0; y[39]=y[41]=0.85; y[38]=y[42]=0.7; y=clean_data(y); yn=segmented_normalization(rt,y); peaks,valleys=find_peaks_and_valleys(rt,yn,threshold=0.1,order=1); print('peaks:',peaks,'valleys:',valleys,'pairs:',find_valley_before_peak_rt(rt,peaks,valleys))"
```

**② Paired peak-valley local baseline subtraction** (`baseline_subtraction.py`):

```bat
python -c "import numpy as np; from baseline_subtraction import paired_peak_valley_baseline_subtraction; mz_p=np.array([100.,101.,102.]); i_p=np.array([500.,800.,600.]); mz_v=np.array([100.5,101.5]); i_v=np.array([100.,200.]); print(paired_peak_valley_baseline_subtraction(mz_p,i_p,mz_v,i_v))"
```

**③ Weighted ¹³C isotope filtering** (`isotope_filtering.py`):

```bat
python -c "from isotope_filtering import advanced_exclude_isotopes; centroids=[(300.0,10000.0),(301.00335,1430.0),(305.0,800.0)]; print(advanced_exclude_isotopes(centroids))"
```

### 3. Run the Included Example Data End-to-End

A small example dataset is bundled in `example_data\` — a 6-cell continuous-flow chromatogram (`example_chromatogram.csv`) and a centroid list (`example_centroids.csv`), no Thermo RAW files needed. One command runs all three algorithm modules on it in the same order as the software's preprocessing pipeline:

```bat
python run_example.py
```

Expected outcome:

1. **Pulse peak screening** — the chromatogram is cleaned and segment-normalized, the 6 cell pulses are detected at rt = 1.0–6.0 min, and each pulse is paired with its preceding baseline valley;
2. **Baseline subtraction** — the valley spectrum is subtracted from the first pulse's MS1 spectrum (total intensity 25774.4 → 1800.0; the four analyte peaks at m/z 100.0 / 100.6 / 101.2 / 102.8 are retained);
3. **Isotope filtering** — the ¹³C isotopic peak at m/z 301.00335 (child of the m/z 300.0 mother ion) is removed from the centroid list.

All intermediate results are written to `example_data\output\`.

---

## File Inventory

| File                               | Description                                                                                      |
| ---------------------------------- | ------------------------------------------------------------------------------------------------ |
| `sc_metaboflow_main.py`            | Main program that links the entire workflow (GUI + preprocessing pipeline)                       |
| `baseline_subtraction.py`          | Paired peak-valley local baseline subtraction algorithm                                          |
| `pulse_peak_screening.py`          | Multi-criteria continuous-flow pulse peak screening                                              |
| `isotope_filtering.py`             | Weighted ¹³C isotope filtering algorithm                                                         |
| `SCMeTA\`                          | Bundled compute backend (code-only copy of the SCMeTA package, located automatically at startup) |
| `r_bridge.py`                      | R bridge module for the optional R Analysis tab                                                  |
| `run_example.py`                   | End-to-end example runner for the three algorithm modules                                        |
| `example_data\`                    | Small example dataset (chromatogram + centroid list) and its `output\` folder                    |
| `dependencies\requirements.txt`    | Pinned Python dependencies (primary)                                                             |
| `dependencies\requirements_v2.txt` | Loose Python dependencies (fallback)                                                             |
| `dependencies\requirements_r.txt`  | Python + R dependencies and R installation notes                                                 |
| `dependencies\install.bat`         | One-click setup script (Windows)                                                                 |
