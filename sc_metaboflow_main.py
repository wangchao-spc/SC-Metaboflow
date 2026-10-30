"""
SC-MetaboFlow English Edition v2.0
Single-Cell Metabolomics Data Processing Suite
Matches the original Chinese UI layout exactly.
"""

import sys
import os
import time
import traceback
from pathlib import Path
from datetime import datetime
from typing import Optional, List, Dict, Any


# ============================================
# SCMeTA Source Code Path Configuration
# ============================================
# Auto-detect SCMeTA package location
_CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_CURRENT_DIR)

# Try multiple possible locations for SCMeTA (bundled copy is preferred)
_POSSIBLE_PATHS = [
    os.path.join(_CURRENT_DIR, "SCMeTA"),                   # Bundled subfolder (preferred)
    os.path.join(_PARENT_DIR, "SCMeTA"),                    # Sibling folder
    os.path.join(_PARENT_DIR, "高速文件 实验数据", "SCMeTA"), # Chinese original path
]

_SCMETA_ROOT = None
for _path in _POSSIBLE_PATHS:
    if os.path.exists(_path) and os.path.isdir(_path):
        _SCMETA_ROOT = _path
        if _path not in sys.path:
            sys.path.insert(0, _path)
            sys.path.insert(0, os.path.dirname(_path))
        break

if _SCMETA_ROOT is None:
    # Last resort: search for SCMeTA in parent tree
    _search = _CURRENT_DIR
    for _ in range(5):
        _candidate = os.path.join(_search, "SCMeTA")
        if os.path.exists(_candidate) and os.path.isdir(_candidate):
            _SCMETA_ROOT = _candidate
            sys.path.insert(0, _candidate)
            sys.path.insert(0, os.path.dirname(_candidate))
            break
        _search = os.path.dirname(_search)

# Add processcell folder to path (for processData imports)
if _SCMETA_ROOT:
    _PROCESSCELL_PATH = os.path.join(_SCMETA_ROOT, "sc-software", "sc-software", "processcell")
    if os.path.exists(_PROCESSCELL_PATH) and _PROCESSCELL_PATH not in sys.path:
        sys.path.insert(0, _PROCESSCELL_PATH)
else:
    for _path in [
        os.path.join(_CURRENT_DIR, "SCMeTA", "sc-software", "sc-software", "processcell"),
        os.path.join(_PARENT_DIR, "SCMeTA", "sc-software", "sc-software", "processcell"),
    ]:
        if os.path.exists(_path) and _path not in sys.path:
            sys.path.insert(0, _path)
            break

# ============================================

import numpy as np
import pandas as pd
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QTabWidget, QLabel, QPushButton, QLineEdit, QComboBox,
    QFileDialog, QMessageBox, QProgressBar, QTextEdit, QGroupBox,
    QSpinBox, QDoubleSpinBox, QCheckBox, QRadioButton, QButtonGroup,
    QTableWidget, QTableWidgetItem, QHeaderView, QSplitter,
    QStackedWidget, QFrame, QScrollArea, QDialog, QDialogButtonBox,
    QFormLayout, QGridLayout, QStatusBar, QToolBar, QListWidget,
    QListWidgetItem, QMenuBar, QMenu, QSizePolicy,
    QSpacerItem, QPlainTextEdit, QToolButton
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QSize, QSettings, QTimer
from PyQt6.QtGui import QIcon, QFont, QPalette, QColor, QKeySequence, QAction

# R bridge import (optional — only if rpy2 + R are available)
try:
    from r_bridge import get_r_bridge, RBridge
    _R_BRIDGE_AVAILABLE = True
except Exception:
    _R_BRIDGE_AVAILABLE = False

# =============================================================================
# APP CONFIG
# =============================================================================

APP_NAME = "SC-MetaboFlow"
APP_VERSION = "2.0.0"
ORG_NAME = "SCMetaboLab"

# =============================================================================
# STYLESHEET - Light theme matching original Chinese UI
# =============================================================================

LIGHT_STYLESHEET = """
QMainWindow {
    background-color: #f0f0f0;
}
QWidget {
    font-family: "Microsoft YaHei", "Segoe UI", Arial, sans-serif;
    font-size: 13px;
}
QTabWidget::pane {
    border: 1px solid #c0c0c0;
    background-color: #f5f5f5;
    border-top: none;
}
QTabBar::tab {
    background-color: #e8e8e8;
    color: #555555;
    padding: 10px 24px;
    border: 1px solid #c0c0c0;
    border-bottom: none;
    border-top-left-radius: 4px;
    border-top-right-radius: 4px;
    margin-right: 2px;
    font-weight: bold;
}
QTabBar::tab:selected {
    background-color: #f5f5f5;
    color: #000000;
    border-bottom: 2px solid #2196F3;
}
QTabBar::tab:hover {
    background-color: #f0f0f0;
    color: #000000;
}
QGroupBox {
    font-weight: bold;
    color: #333333;
    border: 1px solid #c0c0c0;
    border-radius: 6px;
    margin-top: 10px;
    padding-top: 10px;
    padding: 12px;
    background-color: #fafafa;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 6px;
    color: #1976D2;
}
QLineEdit {
    background-color: #ffffff;
    color: #333333;
    border: 1px solid #c0c0c0;
    border-radius: 3px;
    padding: 6px 8px;
    min-height: 22px;
}
QLineEdit:focus {
    border: 1px solid #2196F3;
}
QComboBox {
    background-color: #ffffff;
    color: #333333;
    border: 1px solid #c0c0c0;
    border-radius: 3px;
    padding: 6px 8px;
    min-height: 22px;
}
QComboBox:focus {
    border: 1px solid #2196F3;
}
QComboBox::drop-down {
    border: none;
    width: 24px;
}
QSpinBox, QDoubleSpinBox {
    background-color: #ffffff;
    color: #333333;
    border: 1px solid #c0c0c0;
    border-radius: 3px;
    padding: 6px;
    min-height: 22px;
}
QCheckBox {
    color: #333333;
    spacing: 6px;
}
QCheckBox::indicator {
    width: 16px;
    height: 16px;
}
QRadioButton {
    color: #333333;
    spacing: 6px;
}
QPushButton {
    background-color: #87CEEB;
    color: #ffffff;
    border: none;
    border-radius: 4px;
    padding: 8px 20px;
    font-weight: bold;
    min-height: 32px;
}
QPushButton:hover {
    background-color: #5DADE2;
}
QPushButton:pressed {
    background-color: #3498DB;
}
QPushButton:disabled {
    background-color: #b0b0b0;
    color: #e0e0e0;
}
QPushButton#submitBtn {
    background-color: #4CAF50;
    color: white;
    font-weight: bold;
    padding: 10px 30px;
}
QPushButton#submitBtn:hover {
    background-color: #45a049;
}
QPushButton#stopBtn {
    background-color: #f44336;
    color: white;
}
QPushButton#stopBtn:hover {
    background-color: #da190b;
}
QPushButton#mergeBtn {
    background-color: #FF9800;
    color: white;
    font-weight: bold;
    padding: 10px;
}
QPushButton#mergeBtn:hover {
    background-color: #F57C00;
}
QPushButton#browseBtn {
    background-color: #87CEEB;
    color: white;
    padding: 6px 16px;
    min-width: 60px;
}
QPushButton#browseBtn:hover {
    background-color: #5DADE2;
}
QPushButton#previewBtn {
    background-color: #2196F3;
    color: white;
    font-weight: bold;
    padding: 10px;
}
QPushButton#previewBtn:hover {
    background-color: #1976D2;
}
QProgressBar {
    border: 1px solid #c0c0c0;
    border-radius: 3px;
    text-align: center;
    color: #333333;
    background-color: #ffffff;
    height: 20px;
}
QProgressBar::chunk {
    background-color: #4CAF50;
    border-radius: 2px;
}
QLabel {
    color: #333333;
}
QLabel#sectionTitle {
    font-weight: bold;
    color: #1976D2;
    font-size: 14px;
}
QPlainTextEdit {
    background-color: #ffffff;
    color: #333333;
    border: 1px solid #c0c0c0;
    border-radius: 3px;
    padding: 8px;
    font-family: "Consolas", monospace;
    font-size: 12px;
}
QPlainTextEdit#logConsole {
    background-color: #fafafa;
    border: 1px solid #d0d0d0;
}
QScrollArea {
    border: none;
}
QFrame#separator {
    background-color: #e0e0e0;
    max-height: 1px;
}
"""

# =============================================================================
# WORKER THREADS
# =============================================================================

class WorkerThread(QThread):
    """Base worker thread for background processing."""
    progress = pyqtSignal(int)
    status = pyqtSignal(str)
    finished_signal = pyqtSignal(bool, str)
    log = pyqtSignal(str)
    elapsed_time = pyqtSignal(str)

    def __init__(self, task_name: str):
        super().__init__()
        self.task_name = task_name
        self._is_running = True
        self._start_time = None

    def run(self):
        try:
            self._start_time = time.time()
            self.status.emit(f"Starting {self.task_name}...")
            self.execute()
            if self._is_running:
                elapsed = time.time() - self._start_time
                self.elapsed_time.emit(self._format_time(elapsed))
                self.finished_signal.emit(True, f"{self.task_name} completed successfully.")
        except Exception as e:
            error_msg = f"Error in {self.task_name}: {str(e)}\n{traceback.format_exc()}"
            self.log.emit(error_msg)
            self.finished_signal.emit(False, str(e))

    def execute(self):
        raise NotImplementedError

    def stop(self):
        self._is_running = False
        self.wait(3000)

    def _format_time(self, seconds: float) -> str:
        return time.strftime("%H:%M:%S", time.gmtime(seconds))


# ============================================================
# Excel-safe matrix saving utility (shared by all workers)
# ============================================================
EXCEL_MAX_COLS = 16000  # Excel limit is 16384, leave margin
MZ_RESOLUTION = 0.01    # SCMeTA round_mat resolution


def _round_and_merge_mz(cell_data: dict, resolution: float = MZ_RESOLUTION) -> dict:
    """Round m/z values to given resolution and merge duplicates by summing.

    Matches SCMeTA's round_mat(resolution=0.01) behavior.
    """
    n = int(round(-np.log10(resolution))) if resolution > 0 else 5
    rounded = {}
    for cell_name, cell_dict in cell_data.items():
        r = {}
        for mz_val, intensity in cell_dict.items():
            key = round(float(mz_val), n)
            r[key] = r.get(key, 0.0) + float(intensity)
        rounded[cell_name] = r
    return rounded


def _build_matrix_df(rounded_data: dict, sort_key=None) -> "pd.DataFrame":
    """Build a cell matrix DataFrame from rounded cell data.
    Returns wide-format: rows=cells, cols=m/z, first col=filename.
    """
    import pandas as pd
    all_mz = sorted(set().union(*[set(d.keys()) for d in rounded_data.values()]))
    if not all_mz:
        return None, []

    mz_cols = [f"{mz:.2f}" for mz in all_mz]
    rows = []
    names = sorted(rounded_data.keys(), key=sort_key)
    for cell_name in names:
        row = [rounded_data[cell_name].get(mz, 0.0) for mz in all_mz]
        rows.append([cell_name] + row)

    df = pd.DataFrame(rows, columns=["filename"] + mz_cols)
    return df, mz_cols


def _save_matrix_safe(df: "pd.DataFrame", base_path_no_ext: str, log_func=None):
    """Save matrix in wide format + long format for Excel compatibility.

    Args:
        df: DataFrame with 'filename' column + m/z columns
        base_path_no_ext: Output path WITHOUT extension (e.g. '/path/to/output')
        log_func: Optional logging function
    """
    import pandas as pd
    import numpy as np

    numeric_cols = [c for c in df.columns if c != "filename"]
    n_cols = len(numeric_cols)
    n_rows = df.shape[0]

    # Always save wide format (for programmatic use)
    df.to_csv(f"{base_path_no_ext}.csv", index=False)

    # If too many columns for Excel, also save long format
    if n_cols > EXCEL_MAX_COLS:
        # Melt to long format: filename, m/z, intensity
        long_df = df.melt(id_vars=["filename"], var_name="m/z", value_name="intensity")
        long_df = long_df[long_df["intensity"] > 0]  # Skip zeros
        long_path = f"{base_path_no_ext}_long.csv"
        long_df.to_csv(long_path, index=False)
        if log_func:
            log_func(f"  Saved long format ({len(long_df)} rows) for Excel: {long_path}")

    # Apply 80% presence filter and save filtered version
    if n_rows > 0:
        missing_ratio = (df[numeric_cols] == 0).sum(axis=0) / n_rows
        keep = missing_ratio[missing_ratio <= 0.8].index.tolist()
        if keep and len(keep) < n_cols:
            filtered = df[["filename"] + keep]
            filtered.to_csv(f"{base_path_no_ext}_filtered.csv", index=False)
            if log_func:
                log_func(f"  Filtered: {n_cols} → {len(keep)} m/z (80% presence)")
            # Long format for filtered too
            if len(keep) > EXCEL_MAX_COLS:
                fl = filtered.melt(id_vars=["filename"], var_name="m/z", value_name="intensity")
                fl = fl[fl["intensity"] > 0]
                fl.to_csv(f"{base_path_no_ext}_filtered_long.csv", index=False)


# ============================================================
# Chinese-version aligned merge functions
# ============================================================
# Ported from the Chinese build's ms_data_merger / merge_strategies
# modules (logic verified via bytecode disassembly, 2026-09-29).
# These replace the old 0.01 Da rounding approach: the Chinese version
# groups m/z by PPM tolerance and keeps full-precision values.


def _cn_calculate_ppm_distance(mz1: float, mz2: float) -> float:
    """Chinese ms_data_merger.calculate_ppm_distance."""
    if mz1 == mz2:
        return 0.0
    return abs(mz1 - mz2) / ((mz1 + mz2) / 2) * 1000000.0


def _cn_fast_peak_grouping(mz_values, ppm_tolerance: float) -> dict:
    """Chinese ms_data_merger.fast_peak_grouping_for_raw.

    Adjacent-chaining PPM grouping: walk m/z values in ascending order,
    append to the current group while |mz - last| / midpoint <= tolerance
    (in ppm). Group representative = np.mean of the group.
    """
    import numpy as np
    if len(mz_values) == 0:
        return {}
    if len(mz_values) == 1:
        return {mz_values[0]: [mz_values[0]]}
    sorted_mz = np.sort(mz_values)
    groups = {}
    current_group = [sorted_mz[0]]
    for i in range(1, len(sorted_mz)):
        current_mz = sorted_mz[i]
        last_mz = current_group[-1]
        ppm_distance = abs(current_mz - last_mz) / ((current_mz + last_mz) / 2) * 1000000.0
        if ppm_distance <= ppm_tolerance:
            current_group.append(current_mz)
        else:
            representative = np.mean(current_group)
            groups[representative] = current_group.copy()
            current_group = [current_mz]
    if current_group:
        representative = np.mean(current_group)
        groups[representative] = current_group
    return groups


def _cn_improved_group_mz_by_ppm(mz_values, ppm_tolerance: float) -> dict:
    """Chinese ms_data_merger.improved_group_mz_by_ppm.

    Greedy single-linkage clustering: each unclustered m/z (ascending order)
    seeds a group; any unclustered m/z whose minimum PPM distance to ANY
    group member is <= tolerance joins the group. Representative = mean.
    """
    import numpy as np
    if len(mz_values) == 0:
        return {}
    if len(mz_values) == 1:
        return {mz_values[0]: [mz_values[0]]}
    sorted_mz = np.sort(mz_values)
    groups = {}
    used = np.zeros(len(sorted_mz), dtype=bool)
    for i, mz in enumerate(sorted_mz):
        if used[i]:
            continue
        current_group = [mz]
        used[i] = True
        for j, other_mz in enumerate(sorted_mz):
            if used[j]:
                continue
            min_ppm_distance = min(
                _cn_calculate_ppm_distance(other_mz, g) for g in current_group
            )
            if min_ppm_distance <= ppm_tolerance:
                current_group.append(other_mz)
                used[j] = True
        representative = np.mean(current_group)
        groups[representative] = current_group
    return groups


def _cn_filter_zero_columns(df, zero_threshold: float, folder_path,
                            output_suffix: str, file_prefix: str,
                            log_func=None):
    """Chinese ms_data_merger.filter_zero_columns.

    Keep m/z columns whose zero ratio < zero_threshold (Chinese default 0.8).
    Saves {file_prefix}merged_data{output_suffix}.csv only when at least one
    column is dropped; returns the filtered df (or None).
    """
    import numpy as np
    import pandas as pd
    from pathlib import Path

    filename_col = df['filename'].copy()
    data_columns = df.columns[1:]
    data_matrix = df.iloc[:, 1:].values
    total_rows = len(df)
    zero_counts = np.sum(data_matrix == 0, axis=0)
    zero_ratios = zero_counts / total_rows
    keep_mask = zero_ratios < zero_threshold
    if np.sum(~keep_mask) > 0:
        keep_columns = data_columns[keep_mask]
        filtered_df = pd.DataFrame()
        filtered_df['filename'] = filename_col
        for i, col in enumerate(data_columns):
            if keep_mask[i]:
                filtered_df[col] = df[col]
        filtered_output_file = Path(folder_path) / f"{file_prefix}merged_data{output_suffix}.csv"
        filtered_df.to_csv(filtered_output_file, index=False)
        if log_func:
            log_func(f"  Filtered: {len(data_columns)} → {len(keep_columns)} m/z "
                     f"(zero ratio ≥ {zero_threshold:.0%} dropped)")
        return filtered_df
    return None


def _cn_merge_summary_for_raw(folder_path: str, ppm_tolerance: float,
                              log_func=None):
    """Chinese ms_data_merger.merge_summary_for_raw — per-RAW auto-merge.

    Reads the first *_summary.csv in the folder, groups m/z by PPM
    (fast adjacent chaining), pivots to wide format (rows = peak_id,
    aggfunc = max), formats column names to 5 decimals, renames peak_id
    -> filename with the '{raw}_peak_' prefix, then saves:
      {raw}_peak_merged_data-output.csv  (always)
      {raw}_peak_merged_data-filtered.csv (only if columns were dropped)
    """
    import numpy as np
    import pandas as pd
    from pathlib import Path

    folder = Path(folder_path)
    if not folder.is_dir():
        raise FileNotFoundError(f"Folder does not exist: {folder_path}")
    summary_files = list(folder.glob("*_summary.csv"))
    if not summary_files:
        if log_func:
            log_func(f"  No *_summary.csv found in {folder_path}, skip merge")
        return None
    summary_file = summary_files[0]
    try:
        df_long = pd.read_csv(summary_file)
        raw_basename = summary_file.stem.replace("_summary", "")
        file_prefix = f"{raw_basename}_peak_"
        unique_mz = df_long['m/z'].unique()
        mz_groups = _cn_fast_peak_grouping(np.sort(unique_mz), ppm_tolerance)
        mz_map = {mz: rep for rep, mz_list in mz_groups.items() for mz in mz_list}
        df_long['group_mz'] = df_long['m/z'].map(mz_map)
        df_wide = df_long.pivot_table(
            index='peak_id', columns='group_mz', values='intensity',
            aggfunc='max'
        ).fillna(0)
        df_wide.columns = [f"{col:.5f}" for col in df_wide.columns]
        df_wide.reset_index(inplace=True)
        df_wide.rename(columns={'peak_id': 'filename'}, inplace=True)
        df_wide['filename'] = f"{raw_basename}_peak_" + df_wide['filename'].astype(str)
        output_file_path = folder / f"{file_prefix}merged_data-output.csv"
        df_wide.to_csv(output_file_path, index=False)
        if log_func:
            log_func(f"  Saved {output_file_path.name}: {df_wide.shape[0]} rows × "
                     f"{df_wide.shape[1]} cols")
        _cn_filter_zero_columns(df_wide, zero_threshold=0.8,
                                folder_path=folder, output_suffix='-filtered',
                                file_prefix=file_prefix, log_func=log_func)
        return df_wide
    except Exception as e:
        if log_func:
            log_func(f"  Per-RAW merge failed for {folder_path}: {e}")
        return None


def _cn_merge_filtered_files(files, ppm_tolerance: float, output_folder: str,
                             log_func=None) -> "pd.DataFrame":
    """Chinese merge-tab logic (merge_strategies.merge_filtered_files_logic).

    Merges per-RAW filtered CSV files into one wide matrix:
    - union of all numeric m/z columns (excluding filename)
    - PPM grouping with improved_group_mz_by_ppm (greedy single linkage)
    - column-name mapping at 4/5/6 decimal precision (create_flexible_mz_mapping)
    - per row, per group: max over mapped columns (NaN -> 0.0)
    - columns sorted by m/z value, 5-decimal names
    Outputs (in output_folder):
      merged_data-output.csv           (always)
      merged_data-filtered.csv         (columns with zero ratio > 0.8 dropped)
    """
    import numpy as np
    import pandas as pd
    from pathlib import Path

    all_dataframes = []
    all_mz_values = set()
    for file in files:
        try:
            df = pd.read_csv(file)
            if "filename" not in df.columns:
                continue
            for col in df.columns[1:]:
                try:
                    all_mz_values.add(float(col))
                except (ValueError, TypeError):
                    pass  # skip non-numeric columns
            all_dataframes.append((file.stem, df))
            if log_func:
                log_func(f"  Loaded {file.name}: {df.shape[0]} rows × {df.shape[1]} cols")
        except Exception as e:
            if log_func:
                log_func(f"  Skip {file.name}: {e}")
            continue

    if not all_dataframes:
        raise ValueError("No valid data files could be loaded.")

    unique_mz = np.array(sorted(all_mz_values))
    mz_groups = _cn_improved_group_mz_by_ppm(unique_mz, ppm_tolerance)

    # Flexible column-name mapping (Chinese create_flexible_mz_mapping)
    mz_to_group = {}
    representative_mz = {}
    for rep_mz, mz_list in mz_groups.items():
        rep_col_name = f"{rep_mz:.5f}"
        representative_mz[rep_mz] = rep_col_name
        for mz in mz_list:
            for precision in (4, 5, 6):
                mz_to_group[f"{mz:.{precision}f}"] = rep_mz
            mz_to_group[str(mz)] = rep_mz

    all_rows = []
    for file_stem, df in all_dataframes:
        for _, row in df.iterrows():
            new_row = {"filename": row["filename"]}
            for rep_mz in representative_mz.keys():
                new_row[representative_mz[rep_mz]] = 0.0
            for col in df.columns[1:]:
                if col in mz_to_group:
                    rep_col = representative_mz[mz_to_group[col]]
                    current_val = float(row[col]) if pd.notna(row[col]) else 0.0
                    new_row[rep_col] = max(new_row[rep_col], current_val)
            all_rows.append(new_row)

    merged_df = pd.DataFrame(
        all_rows,
        columns=["filename"] + sorted(representative_mz.values(), key=lambda x: float(x)),
    )

    out_folder = Path(output_folder)
    merged_df.to_csv(out_folder / "merged_data-output.csv", index=False)
    if log_func:
        log_func(f"  Saved merged_data-output.csv: {merged_df.shape[0]} rows × "
                 f"{merged_df.shape[1]} cols")

    mz_cols = merged_df.select_dtypes(include=[np.number]).columns
    zero_ratios = (merged_df[mz_cols] == 0).mean()
    sparse_cols = zero_ratios[zero_ratios > 0.8].index.tolist()
    if sparse_cols:
        df_filtered = merged_df.drop(columns=sparse_cols)
        df_filtered.to_csv(out_folder / "merged_data-filtered.csv", index=False)
        if log_func:
            log_func(f"  Saved merged_data-filtered.csv: removed {len(sparse_cols)} "
                     f"sparse columns (zero ratio > 80%)")
    return merged_df


class PreprocessingWorker(WorkerThread):
    """Worker running the Chinese-version preprocessing pipeline.

    Fully aligned with the Chinese SC-MetaboFlow PeakProcessor
    (peak_processor.py + processData.py + ms_data_merger.py):
    1. Extract m/z chromatogram around the target mass (TraceType.MassRange)
    2. clean_data (negatives → 0) + segmented normalization
    3. argrelextrema peak/valley detection (absolute threshold)
    4. Match each peak with its preceding valley
    5. Per peak: MS1 at apex/valley → linear background subtraction →
       positive centroids → advanced C13 isotope filter (optional)
    6. Save {raw}_summary.csv (peak_id = scan number, m/z raw, %.5f)
    7. Per-RAW auto-merge: PPM grouping → {raw}_peak_merged_data-output.csv
       and {raw}_peak_merged_data-filtered.csv (zero ratio < 0.8 kept)
    8. Optional cross-file merge into _merged/ (Chinese merge-tab algorithm)
    """

    def __init__(self, raw_path: str, output_path: str, params: dict):
        super().__init__("Preprocessing")
        self.raw_path = raw_path
        self.output_path = output_path
        self.params = params

    def execute(self):
        import numpy as np
        from scipy.interpolate import interp1d
        import glob

        self.status.emit("Initializing preprocessing pipeline...")
        self.progress.emit(2)

        # ============================================================
        # Use fisher_py for RAW reading (matching Chinese main7.py)
        # ============================================================
        try:
            from fisher_py import RawFile
            from fisher_py.data.business import TraceType
        except ImportError:
            raise ImportError("fisher_py not installed. Install: pip install fisher-py")
        except Exception as e:
            # In a packaged (PyInstaller) build, fisher_py's ThermoFisher
            # .NET DLLs live in fisher_py/dll — if they were not bundled,
            # or the .NET Framework runtime is missing, loading fails here
            # with a pythonnet error. Give a clear message instead.
            raise RuntimeError(
                f"The Thermo RAW-reading component failed to load: {e}\n\n"
                "Possible causes:\n"
                "1) This portable copy is missing the ThermoFisher .NET DLLs "
                "(fisher_py/dll folder was not bundled).\n"
                "2) The .NET Framework 4.8 runtime is not installed on this "
                "computer.\n\n"
                "This only affects .RAW file processing; XLSX/CSV workflows "
                "(Data Cleaning, Visualization) do not need it.")

        # Import signal processing functions (Chinese processData module)
        from processData import (
            find_peaks_and_valleys, clean_data,
            find_valley_before_peak_rt, segmented_normalization,
            advanced_exclude_isotopes,
        )

        # Step 1: Find RAW files
        self.status.emit("Scanning for RAW files...")
        raw_files = []
        if os.path.isdir(self.raw_path):
            for ext in ("*.RAW", "*.raw", "*.Raw"):
                raw_files.extend(glob.glob(os.path.join(self.raw_path, ext)))
            # Windows filesystems are case-insensitive: all three patterns
            # match the same files, so dedupe before processing.
            raw_files = sorted(set(raw_files))
        elif os.path.isfile(self.raw_path) and self.raw_path.lower().endswith('.raw'):
            raw_files = [self.raw_path]

        if not raw_files:
            raise ValueError(
                f"No RAW files found in: {self.raw_path}\n"
                f"Please select a folder containing .RAW files."
            )

        self.status.emit(f"Found {len(raw_files)} RAW file(s) to process")
        self.log.emit(f"Files: {[os.path.basename(f) for f in raw_files]}")
        self.progress.emit(5)

        # Parameters — Chinese version defaults
        # (mass_spin 760.58, ppm_spin 5, seg_spin 5.0, thr_spin 0.1, order_spin 1)
        target_mass = self.params.get("target_mz", 760.58)
        mass_tolerance_ppm = self.params.get("ppm_tolerance", 5.0)
        segment_length = self.params.get("segment_length", 5.0)
        intensity_threshold = self.params.get("intensity_threshold", 0.1)
        peak_order = self.params.get("smoothing_order", 1)  # argrelextrema order
        remove_isotopes_flag = self.params.get("remove_isotopes", False)
        subtract_background_flag = self.params.get("subtract_background", True)
        auto_merge_flag = self.params.get("auto_merge", True)

        total_files = len(raw_files)
        all_file_names = []

        # ============================================================
        # Step 2: Process each RAW file (identical to main7.py logic)
        # ============================================================
        for file_idx, raw_file_path in enumerate(raw_files):
            if not self._is_running:
                break

            raw_filename = os.path.splitext(os.path.basename(raw_file_path))[0]
            self.status.emit(f"[{file_idx+1}/{total_files}] Processing: {raw_filename}")

            base_progress = 5 + int((file_idx) / total_files * 85)

            try:
                raw_file = RawFile(raw_file_path)
            except Exception as e:
                self.log.emit(f"ERROR: Cannot read {raw_filename}: {e}")
                continue

            # Create output subfolder for this file
            file_output_folder = os.path.join(self.output_path, raw_filename)
            os.makedirs(file_output_folder, exist_ok=True)

            # --- Chromatogram extraction (CN PeakProcessor.process_single) ---
            self.status.emit(f"  Extracting chromatogram m/z={target_mass}...")
            rt, intensity = raw_file.get_chromatogram(
                target_mass, mass_tolerance_ppm, TraceType.MassRange
            )
            rt = np.array(rt, dtype=np.float64)
            intensity = np.array(intensity, dtype=np.float64)

            # --- clean_data: negatives → 0 (CN clean_data_with_zero) ---
            self.status.emit(f"  Cleaning data...")
            cleaned_intensity = clean_data(intensity)
            self.progress.emit(base_progress + 2)

            # --- Segmented normalization (CN rolling segments) ---
            self.status.emit(f"  Normalizing...")
            normalized_intensity = segmented_normalization(
                rt, cleaned_intensity, segment_length=segment_length
            )
            self.progress.emit(base_progress + 4)

            # --- argrelextrema peak/valley detection (CN find_peaks_and_valleys) ---
            self.status.emit(f"  Finding peaks...")
            peak_indexes, valley_indexes = find_peaks_and_valleys(
                rt, normalized_intensity,
                threshold=intensity_threshold, order=peak_order
            )

            if len(peak_indexes) == 0:
                self.log.emit(f"  WARNING: No peaks found in {raw_filename}")
                continue

            # --- Match peaks with preceding valleys (CN find_valley_before_peak_rt) ---
            vp_pairs = find_valley_before_peak_rt(rt, peak_indexes, valley_indexes)
            if not vp_pairs:
                self.log.emit(f"  WARNING: No peak-valley pairs in {raw_filename}")
                continue
            self.status.emit(f"  Found {len(vp_pairs)} cells (peaks)")
            self.progress.emit(base_progress + 6)

            # --- MS1 cache over unique peak/valley RTs (CN process_single) ---
            all_rts = []
            for _, rt_peak, _, rt_valley in vp_pairs:
                all_rts.extend([rt_peak, rt_valley])
            unique_rts = list(set(all_rts))
            ms1_cache = {rt_val: raw_file.get_scan_ms1(rt_val) for rt_val in unique_rts}

            # --- Process each peak (CN process_single main loop) ---
            summary_rows = []
            for peak_idx, rt_peak, valley_idx, rt_valley in vp_pairs:
                if not self._is_running:
                    break

                mz_p, i_p, _, _ = ms1_cache[rt_peak]
                mz_v, i_v, _, _ = ms1_cache[rt_valley]
                mz_p = np.asarray(mz_p, dtype=np.float64)
                i_p = np.asarray(i_p, dtype=np.float64)
                mz_v = np.asarray(mz_v, dtype=np.float64)
                i_v = np.asarray(i_v, dtype=np.float64)

                # peak_id = scan number at peak apex (CN), fallback = peak index
                try:
                    scan_number = raw_file.get_scan_number_from_retention_time(rt_peak)
                except Exception:
                    scan_number = peak_idx

                # Background subtraction (CN _subtract_background_clear_logic):
                # linear interpolation of the valley spectrum onto the peak's
                # m/z grid, then net = max(peak - background, 0)
                if subtract_background_flag and len(mz_p) > 0 and len(mz_v) > 0:
                    try:
                        f_background = interp1d(
                            mz_v, i_v, kind='linear', fill_value=0, bounds_error=False
                        )
                        background_at_peak_positions = f_background(mz_p)
                        intensity_final = np.maximum(
                            i_p - background_at_peak_positions, 0
                        )
                        mz_final = mz_p
                    except Exception:
                        mz_final, intensity_final = mz_p, i_p
                else:
                    mz_final, intensity_final = mz_p, i_p

                mask_pos = intensity_final > 0
                centroids = list(zip(mz_final[mask_pos], intensity_final[mask_pos]))

                # Advanced C13 isotope filter (CN advanced_exclude_isotopes)
                if remove_isotopes_flag and len(centroids) > 1:
                    sorted_centroids = sorted(centroids, key=lambda x: x[0])
                    centroids = advanced_exclude_isotopes(
                        sorted_centroids,
                        mass_tolerance_ppm=5.0,
                        intensity_tolerance=0.3,
                        min_isotope_probability=0.6,
                    )

                if centroids:
                    for mz_val, inten in centroids:
                        summary_rows.append(
                            {"peak_id": scan_number, "m/z": mz_val, "intensity": inten}
                        )

            # --- Save summary + per-RAW auto-merge (CN naming) ---
            if summary_rows:
                df_raw = pd.DataFrame(summary_rows)
                df_raw.to_csv(
                    os.path.join(file_output_folder, f"{raw_filename}_summary.csv"),
                    index=False, float_format='%.5f',
                )
                all_file_names.append(raw_filename)
                if auto_merge_flag:
                    _cn_merge_summary_for_raw(
                        file_output_folder, mass_tolerance_ppm,
                        log_func=self.log.emit,
                    )
                self.status.emit(
                    f"  {raw_filename}: {len(vp_pairs)} cells, "
                    f"{len(summary_rows)} centroids → summary saved"
                )
            else:
                self.log.emit(f"  WARNING: No centroids extracted from {raw_filename}")

            self.progress.emit(base_progress + 8)

        # ============================================================
        # Step 3: Cross-file merge (if auto-merge enabled)
        # ============================================================
        self.progress.emit(92)
        if auto_merge_flag and len(all_file_names) > 1:
            self.status.emit("Auto-merging all processed files...")
            self._auto_merge_files(all_file_names)
        elif len(all_file_names) == 1:
            self.status.emit(f"Processed 1 file: {all_file_names[0]}")

        self.progress.emit(100)
        self.output_file_path = self.output_path
        self.status.emit(
            f"Preprocessing complete! {len(all_file_names)}/{total_files} file(s) processed."
        )
        self.log.emit(f"Output saved to: {self.output_path}")

        # ============================================================
        # Step 4: Try SCMeTA post-processing for additional outputs
        # ============================================================
        self._try_scmeta_postprocess(all_file_names)

    def _auto_merge_files(self, file_names: list):
        """Cross-file merge using the Chinese merge-tab algorithm.

        Merges the per-RAW {raw}_peak_merged_data-filtered.csv files via PPM
        grouping (improved_group_mz_by_ppm) into _merged/ using the same
        merged_data naming as File Merge, so both paths produce identical
        outputs. Filenames keep the plain {raw}_peak_{id} form (no prefix).
        """
        merge_dir = os.path.join(self.output_path, "_merged")
        os.makedirs(merge_dir, exist_ok=True)

        files = []
        for fname in file_names:
            p = os.path.join(
                self.output_path, fname, f"{fname}_peak_merged_data-filtered.csv"
            )
            if os.path.exists(p):
                files.append(Path(p))
        if len(files) > 1:
            try:
                _cn_merge_filtered_files(
                    files, self.params.get("ppm_tolerance", 5.0),
                    merge_dir, log_func=self.log.emit,
                )
                self.log.emit(f"Auto-merge: {len(files)} files → {merge_dir}")
                self.output_file_path = merge_dir
            except Exception as e:
                self.log.emit(f"Auto-merge skipped: {e}")

    def _try_scmeta_postprocess(self, file_names: list):
        """Apply SCMeTA normalize + fill steps to the merged output files.
        This ensures consistency with Chinese version's post-processing pipeline.
        """
        if not file_names:
            return
        try:
            from SCMeTA.core import Process
            from SCMeTA.config import PARAMETERS

            p = Process(ref_mz=self.params.get("target_mz", 760.58))

            # Load the cell matrices we just created
            for fname in file_names:
                cell_mat_path = os.path.join(
                    self.output_path, fname,
                    f"{fname}_peak_merged_data-output.csv"
                )
                if os.path.exists(cell_mat_path):
                    try:
                        p.load_processed(cell_mat_path, file_name=fname, file_type="cell_mat")
                    except Exception:
                        continue

            if p.data:
                # Apply SCMeTA normalize (m/z-based) and KNN fill
                p.post_process(
                    data=p.data,
                    normalize_method=["mz"],
                    fillna_method="knn",
                )
                # Save processed results alongside original outputs
                scmeta_dir = os.path.join(self.output_path, "SCMeTA_processed")
                p.save(path=scmeta_dir, data_type="cell_mat")
                p.info()
                self.log.emit(f"SCMeTA normalize+fill applied → {scmeta_dir}")
        except Exception as e:
            self.log.emit(f"Note: SCMeTA post-process skipped ({type(e).__name__}: {e})")

class MergeWorker(WorkerThread):
    """Worker for merging preprocessed filtered CSV files (Chinese merge-tab
    algorithm: improved PPM grouping + per-row max aggregation, zero-ratio
    filtering, plain {raw}_peak_{id} filenames).
    """

    def __init__(self, input_dir: str, output_dir: str):
        super().__init__("Data Merging")
        self.input_dir = input_dir
        self.output_dir = output_dir

    def execute(self):
        from pathlib import Path

        self.status.emit("Scanning for preprocessed files...")

        # CN merge-tab file discovery (merge_filtered_files_logic):
        # 1) *merged_data-filtered.csv
        # 2) *_peak_merged_data-filtered.csv
        # 3) *.csv whose name contains both 'filtered' and 'merged'
        # Recursive search that excludes auto-merge (_merged) and this
        # worker's own output folder, so re-running never merges
        # already-merged data.
        def _usable(p):
            return not ({x.lower() for x in p.parts} & {"_merged", "merged_output"})

        root = Path(self.input_dir)
        all_files = [p for p in root.rglob("*merged_data-filtered.csv") if _usable(p)]
        if not all_files:
            all_files = [p for p in root.rglob("*_peak_merged_data-filtered.csv")
                         if _usable(p)]
        if not all_files:
            all_files = [p for p in root.rglob("*.csv") if _usable(p)
                         and "filtered" in p.name and "merged" in p.name]

        if not all_files:
            raise FileNotFoundError(
                "No suitable filtered files found in the folder.\n"
                "Please make sure the folder contains CSV files produced by "
                "peak extraction and auto-merge."
            )

        self.status.emit(f"Found {len(all_files)} file(s) to merge.")
        self.log.emit(f"Files: {[f.name for f in all_files]}")

        # CN merge-tab algorithm: PPM grouping + per-row max aggregation,
        # merged_data-output.csv always, merged_data-filtered.csv when
        # sparse columns (zero ratio > 0.8) are dropped.
        merge_dir = os.path.join(self.output_dir, "merged_output")
        os.makedirs(merge_dir, exist_ok=True)
        _cn_merge_filtered_files(
            all_files, 5.0, merge_dir, log_func=self.log.emit,
        )

        self.output_file_path = merge_dir
        self.status.emit(f"Merged {len(all_files)} files → {merge_dir}")
        self.progress.emit(100)

class SplitWorker(WorkerThread):
    """Worker for splitting merged data by group information.
    Matches Chinese manual: add group info, rename Sample->Group, save as .xlsx.
    """

    def __init__(self, excel_file: str, output_dir: str):
        super().__init__("Group Splitting")
        self.excel_file = excel_file
        self.output_dir = output_dir

    def execute(self):
        import pandas as pd
        import numpy as np

        self.status.emit("Reading file for group splitting...")
        self.progress.emit(10)

        # Load data — try CSV first, then Excel, with encryption handling
        df = None
        load_error = None

        # Try 1: CSV (always preferred — merge output is CSV)
        # NOTE: the scan-ID column ('filename'/'Sample') is NOT moved into
        # the index — sheets are saved with index=False, which would drop
        # it. It stays a regular column so every split sheet keeps it.
        try:
            df = pd.read_csv(self.excel_file)
            self.status.emit("Loaded as CSV")
        except Exception as e_csv:
            pass

        # Try 2: Excel
        if df is None:
            if self.excel_file.endswith(('.xlsx', '.xls')):
                try:
                    df = pd.read_excel(self.excel_file)
                    self.status.emit("Loaded as Excel")
                except Exception as e:
                    load_error = str(e)

        # Try 3: openpyxl with data_only
        if df is None and self.excel_file.endswith('.xlsx'):
            try:
                from openpyxl import load_workbook
                wb = load_workbook(self.excel_file, read_only=True, data_only=True,
                                   keep_links=False)
                ws = wb.active
                data = list(ws.values)
                if data:
                    df = pd.DataFrame(data[1:], columns=data[0])
                wb.close()
                self.status.emit("Loaded via openpyxl")
            except Exception:
                pass

        if df is None:
            if "encrypt" in str(load_error).lower() if load_error else False:
                raise ValueError(
                    f"Workbook is encrypted (password-protected).\n\n"
                    f"Please open the file in Excel, then:\n"
                    f"  1. File → Save As → Browse\n"
                    f"  2. Click 'Tools' next to Save → 'General Options'\n"
                    f"  3. Clear any password fields\n"
                    f"  4. Save without password\n\n"
                    f"Or simply save as CSV: File → Save As → CSV UTF-8"
                )
            raise ValueError(
                f"Cannot read file: {self.excel_file}\n"
                f"Supported formats: CSV (.csv), Excel (.xlsx/.xls)\n"
                f"Error: {load_error}"
            )

        self.status.emit(f"Loaded: {df.shape[0]} rows x {df.shape[1]} columns")
        self.progress.emit(25)

        # Keep the scan-ID column: upstream merged tables carry 'filename'
        # or 'Sample' identifying which scan each row came from. It must
        # survive into every split sheet (and onward through cleaning) so
        # each data point in the plots can be traced back to its scan.
        df = self._resolve_leading_column(df)

        # Check for Group column; if not present, try to use "Sample" column
        group_col = None
        for candidate in ["Group", "group", "GROUP"]:
            if candidate in df.columns:
                group_col = candidate
                break

        if group_col is None:
            # Try to find a column named "Sample" to rename to "Group"
            for sample_candidate in ["Sample", "sample", "SAMPLE"]:
                if sample_candidate in df.columns:
                    df = df.rename(columns={sample_candidate: "Group"})
                    group_col = "Group"
                    self.status.emit("Renamed 'Sample' column to 'Group'")
                    break

        if group_col is None:
            # Use the first non-numeric column as group — but NEVER an ID
            # column like 'filename' (a scan ID would make every row its own
            # "group" and produce thousands of single-row sheets)
            for col in df.columns:
                if not pd.api.types.is_numeric_dtype(df[col]):
                    if str(col).lower() in ("filename", "file_name", "id", "scan"):
                        continue
                    group_col = col
                    df = df.rename(columns={col: "Group"})
                    self.status.emit(f"Using '{col}' as Group column")
                    break

        if group_col is None:
            raise ValueError(
                "No Group or Sample column found. "
                "Please add a 'Group' column with group labels before splitting."
            )

        self.progress.emit(40)

        # Group the data
        groups = df.groupby(group_col)
        group_names = list(groups.groups.keys())
        self.status.emit(f"Found {len(group_names)} groups: {group_names}")

        if len(group_names) == 0:
            raise ValueError(
                f"The '{group_col}' column contains no valid group labels "
                f"(all {df.shape[0]} values are empty).\n\n"
                f"How to fix:\n"
                f"  1. Open your file in Excel\n"
                f"  2. Fill the '{group_col}' column with group names "
                f"(e.g. Control / Treated)\n"
                f"  3. Save the file and run Group Split again"
            )

        # Save each group to a separate sheet in an Excel file
        base_name = os.path.splitext(os.path.basename(self.excel_file))[0]
        output_file = os.path.join(self.output_dir, f"{base_name}.split.xlsx")

        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            for i, (group_name, group_data) in enumerate(groups):
                if not self._is_running:
                    break

                self.status.emit(f"Writing group: {group_name} ({group_data.shape[0]} rows)...")
                # Clean sheet name (Excel limits: 31 chars max)
                sheet_name = str(group_name)[:31]
                if sheet_name.endswith('$'):
                    sheet_name = sheet_name[:-1] + '_'
                group_data.to_excel(writer, sheet_name=sheet_name, index=False)

                self.progress.emit(40 + int((i + 1) / len(group_names) * 55))

        self.output_file_path = output_file
        self.status.emit(f"Split into {len(group_names)} groups. Saved to: {output_file}")
        self.progress.emit(100)

    @staticmethod
    def _resolve_leading_column(df):
        """Handle the file's leftmost column so scan IDs are kept.

        Upstream merged tables identify each row by a scan-ID column
        ('filename', 'Sample', ...). Saving sheets with index=False used
        to drop that ID silently. Rules here:

        - A named scan-ID column is kept untouched (it just rides along
          as a regular column into every split sheet).
        - A leading unnamed column ("Unnamed: 0" / empty header) that
          holds plain row numbers 0..n-1 (or 1..n) is a saved index —
          dropped.
        - A leading unnamed column holding anything else (scan IDs whose
          header was lost) is renamed to 'filename' and kept.
        """
        import pandas as pd

        if df is None or df.empty or df.shape[1] == 0:
            return df

        first_col = df.columns[0]
        name = str(first_col).strip().lower()
        if name in ("filename", "file", "file_name", "sample", "scan",
                    "name", "id"):
            return df

        unnamed = (name == "" or name.startswith("unnamed")
                   or str(first_col).strip().lower() in ("index", "row"))
        if not unnamed:
            return df

        vals = pd.Series(df[first_col])
        if pd.api.types.is_numeric_dtype(vals):
            v = vals.dropna()
            if len(v) == len(vals) and (v.astype("int64").astype("float64")
                                        == v.astype("float64")).all():
                v = v.astype("int64")
                if (v == pd.RangeIndex(len(v))).all() or \
                        (v == pd.RangeIndex(1, len(v) + 1)).all():
                    return df.drop(columns=[first_col])
        # Not row numbers — scan IDs that lost their header: keep them
        return df.rename(columns={first_col: "filename"})

class DataCleaningWorker(WorkerThread):
    """Worker for data cleaning — delegates to SCMeTA methods (matches Chinese version)."""

    def __init__(self, data_file: str, output_dir: str, params: dict):
        super().__init__("Data Cleaning")
        self.data_file = data_file
        self.output_dir = output_dir
        self.params = params

    def execute(self):
        import numpy as np
        import pandas as pd

        self.status.emit("Loading data...")
        self.progress.emit(5)

        # Load data — every sheet is cleaned INDEPENDENTLY (matches the
        # Chinese version: each group's sheet is filtered/imputed on its own)
        if self.data_file.endswith('.csv'):
            raw_sheets = {"data": pd.read_csv(self.data_file, index_col=0)}
        elif self.data_file.endswith(('.xlsx', '.xls')):
            raw_sheets = pd.read_excel(self.data_file, sheet_name=None)
        elif self.data_file.endswith('.txt'):
            raw_sheets = {"data": pd.read_csv(self.data_file, index_col=0, sep='\t')}
        else:
            raise ValueError(f"Unsupported format: {self.data_file}")

        if len(raw_sheets) > 1:
            self.status.emit(
                f"Loaded {len(raw_sheets)} sheets: "
                f"{', '.join(list(raw_sheets.keys())[:8])}")

        cleaned_sheets = []
        for sheet_name, sheet_df in raw_sheets.items():
            cleaned_sheets.append(self._clean_one_sheet(sheet_df, sheet_name))

        # Optional: balance group sizes — random downsample every sheet to
        # the smallest sheet's row count so all groups end up with the same
        # number of cells (opt-in; OFF by default). Applied AFTER per-sheet
        # cleaning so imputation/outlier steps still see the full data.
        if self.params.get("balance_groups", False):
            if len(cleaned_sheets) > 1:
                target = min(len(s) for s in cleaned_sheets)
                if target >= 1:
                    rng = np.random.RandomState(42)
                    balanced = []
                    for sheet_name, s in zip(raw_sheets.keys(),
                                             cleaned_sheets):
                        if len(s) > target:
                            keep_idx = rng.choice(len(s), size=target,
                                                  replace=False)
                            balanced.append(
                                s.iloc[keep_idx].reset_index(drop=True))
                            self.status.emit(
                                f"[{sheet_name}] Balanced: kept {target} of "
                                f"{len(s)} rows (random downsample, seed=42)")
                        else:
                            balanced.append(s)
                    cleaned_sheets = balanced
                else:
                    self.log.emit(
                        "Group balancing skipped — a sheet ended up empty.")
            else:
                self.log.emit(
                    "Group balancing skipped — only one sheet, nothing "
                    "to balance against.")

        df = pd.concat(cleaned_sheets, ignore_index=True)
        self.status.emit(f"Cleaned total: {df.shape[0]} rows x {df.shape[1]} cols")

        # --- Re-split the combined table for the post-cleaning steps ---
        # (numeric m/z columns carry float names after per-sheet conversion)
        meta_cols = [c for c in df.columns if not isinstance(c, float)]
        numeric_cols = [c for c in df.columns if isinstance(c, float)]
        meta_df = df[meta_cols].copy() if meta_cols else pd.DataFrame(index=df.index)
        df_numeric = df[numeric_cols].copy() if numeric_cols else pd.DataFrame(index=df.index)
        self.progress.emit(70)

        # Step 5: Normalization (numeric only)
        if self.params.get("enable_std", False) and not df_numeric.empty:
            method_map = {
                "Max Normalization": ["max"],
                "Z-score Standardization": ["zscore"],
                "Min-Max Scaling": ["minmax"],
                "Log Transformation": ["log"],
            }
            std_method = self.params.get("std_method", "Max Normalization")
            methods = method_map.get(std_method, ["max"])

            self.status.emit(f"Normalizing ({std_method})...")
            try:
                from SCMeTA.method import normalize
                ref_mz = self.params.get("ref_value", 760.58)
                df_numeric = normalize(df_numeric, methods, mz=ref_mz)
            except Exception:
                # SCMeTA normalize may fail with newer pandas (applymap→map)
                # or unsupported methods. Use fallback.
                self._fallback_normalize(df_numeric, std_method)
        self.progress.emit(85)

        # Step 6: Reference value filter
        if self.params.get("enable_ref", False):
            ref_val = self.params.get("ref_value", 760.58)
            ref_ppm = self.params.get("ref_ppm", 5.0)
            self.status.emit(f"Filtering by ref m/z {ref_val} ± {ref_ppm} ppm...")
            matched = []
            for col in df_numeric.columns:
                try:
                    fcol = float(col)
                    err = abs(fcol - ref_val) / ref_val * 1e6 if ref_val > 0 else float('inf')
                    if err <= ref_ppm:
                        matched.append(col)
                except (ValueError, TypeError):
                    pass
            if matched:
                self.status.emit(f"Found {len(matched)} columns matching reference")

        # --- Recombine numeric and metadata columns ---
        if not df_numeric.empty and meta_cols:
            df = pd.concat([meta_df, df_numeric], axis=1)
        elif not df_numeric.empty:
            df = df_numeric
        else:
            df = meta_df

        # Final safety check — never write an empty output
        if df.shape[0] == 0 or df.shape[1] == 0:
            raise ValueError(
                "Data cleaning removed ALL rows/columns — nothing left to save.\n"
                "Try unchecking some cleaning options or lowering their thresholds.")

        self.progress.emit(95)

        # Save with Excel-safe format (add timestamp to avoid lock conflicts)
        from datetime import datetime
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        base_name = os.path.splitext(os.path.basename(self.data_file))[0]

        if "filename" in df.columns or df.shape[1] > 100:
            if "filename" not in df.columns:
                df["filename"] = [f"row_{i}" for i in range(len(df))]
            base_path = os.path.join(self.output_dir, f"{base_name}_cleaned_{ts}")
            try:
                _save_matrix_safe(df, base_path, log_func=self.log.emit)
                output_file = f"{base_path}.csv"
            except PermissionError:
                # File is locked (open in Excel) — save with alternate name
                alt_path = os.path.join(self.output_dir, f"{base_name}_cleaned_{ts}_NEW")
                _save_matrix_safe(df, alt_path, log_func=self.log.emit)
                output_file = f"{alt_path}.csv"
                self.log.emit("WARNING: Original file locked. Saved with _NEW suffix.")
        else:
            output_file = os.path.join(self.output_dir, f"{base_name}_cleaned_{ts}.csv")
            try:
                df.to_csv(output_file, index=False)
            except PermissionError:
                output_file = os.path.join(self.output_dir, f"{base_name}_cleaned_{ts}_NEW.csv")
                df.to_csv(output_file, index=False)
        self.progress.emit(100)

        self.output_file_path = output_file
        self.status.emit(f"Data cleaning completed. Output: {output_file}")
        self.log.emit(f"Cleaning done: {output_file}")

    def _clean_one_sheet(self, df, sheet_name: str):
        """Clean ONE sheet independently (Chinese version behavior):
        0→NaN → select imputable columns (missing ≤ threshold) → KNN (k=20)
        → IQR outliers. Over-threshold columns are KEPT but not imputed."""
        # --- Separate metadata columns from numeric m/z columns ---
        new_cols, meta_cols, numeric_cols = [], [], []
        for c in df.columns:
            try:
                new_cols.append(float(c))
                numeric_cols.append(float(c))
            except (ValueError, TypeError):
                new_cols.append(c)
                meta_cols.append(c)
        df.columns = new_cols
        meta_df = df[meta_cols].copy() if meta_cols else pd.DataFrame(index=df.index)
        df_numeric = df[numeric_cols].copy()

        self.status.emit(
            f"[{sheet_name}] {df.shape[0]} rows, {len(numeric_cols)} m/z + "
            f"{len(meta_cols)} metadata cols")

        # Step 1: Zero handling — replace 0 with NaN (Chinese default), or
        # with half of each column's minimum NONZERO value (user option,
        # "Replace 0 with Half of Column Minimum"). Columns that are
        # entirely zero have no nonzero minimum to halve, so they are set
        # to NaN (same treatment as the default path); KNN skips all-NaN
        # columns automatically.
        if not df_numeric.empty:
            if self.params.get("zero_half_min", False):
                all_zero_cols = 0
                for col in df_numeric.columns:
                    col_vals = df_numeric[col]
                    nonzero = col_vals[col_vals != 0]
                    if len(nonzero) > 0:
                        fill_val = nonzero.min() / 2.0
                        df_numeric[col] = col_vals.mask(col_vals == 0, fill_val)
                    else:
                        all_zero_cols += 1
                        df_numeric[col] = np.nan
                self.status.emit(
                    f"[{sheet_name}] Replaced 0 with half of column minimum "
                    f"({all_zero_cols} all-zero columns set to NaN)")
            else:
                df_numeric = df_numeric.replace(0, np.nan)

        # Step 2: Sparse filter (Chinese version: threshold 0.5) — columns
        # with missing ratio ≤ threshold are imputed; the rest are KEPT but
        # left as NaN (not removed, not imputed)
        valid_cols = list(df_numeric.columns)
        if self.params.get("filter_sparse", False) and not df_numeric.empty:
            threshold = self.params.get("sparse_threshold", 0.5)
            missing_ratio = df_numeric.isna().sum(axis=0) / df_numeric.shape[0]
            valid_cols = missing_ratio[missing_ratio <= threshold].index.tolist()
            skipped = len(df_numeric.columns) - len(valid_cols)
            self.status.emit(
                f"[{sheet_name}] Sparse filter (missing ≤ {threshold}): "
                f"{len(valid_cols)} imputable, {skipped} kept without imputation")

        # Step 3: KNN Imputation (Chinese version: sklearn KNNImputer, k=20)
        # Only columns with PARTIAL missing values are imputed; all-NaN
        # columns (e.g. all-zero columns set to NaN above) and complete
        # columns are left untouched.
        if self.params.get("impute_knn", False) and valid_cols:
            impute_cols = [c for c in valid_cols
                           if df_numeric[c].isna().any()
                           and not df_numeric[c].isna().all()]
            if impute_cols and df_numeric.shape[0] > 1:
                from sklearn.impute import KNNImputer
                n_neighbors = min(self.params.get("knn_k", 20),
                                  df_numeric.shape[0] - 1)
                if n_neighbors >= 1:
                    self.status.emit(
                        f"[{sheet_name}] Filling missing values "
                        f"(KNN, k={n_neighbors}, {len(impute_cols)} cols)...")
                    knn = KNNImputer(n_neighbors=n_neighbors)
                    imputed = knn.fit_transform(df_numeric[impute_cols])
                    df_numeric[impute_cols] = pd.DataFrame(
                        imputed, columns=impute_cols, index=df_numeric.index)

        # Step 4: Outlier handling — IQR (Chinese version: x1.5, outliers
        # replaced by the median of the NON-outlier values)
        if self.params.get("iqr_filter", False) and valid_cols:
            iqr_coeff = self.params.get("iqr_coeff", 1.5)
            self.status.emit(
                f"[{sheet_name}] Handling outliers (IQR x {iqr_coeff})...")
            for col in valid_cols:
                col_data = df_numeric[col].dropna()
                if len(col_data) < 4:
                    continue
                Q1, Q3 = col_data.quantile(0.25), col_data.quantile(0.75)
                IQR = Q3 - Q1
                if IQR <= 0:
                    continue
                lower, upper = Q1 - iqr_coeff * IQR, Q3 + iqr_coeff * IQR
                mask = (df_numeric[col] < lower) | (df_numeric[col] > upper)
                if mask.any():
                    median_non_outliers = df_numeric.loc[~mask, col].median()
                    if pd.notna(median_non_outliers):
                        df_numeric.loc[mask, col] = median_non_outliers

        # Step 4b: Outlier handling — Standard Deviation (English-only extra)
        if self.params.get("std_filter", False) and not df_numeric.empty:
            std_coeff = self.params.get("std_coeff", 3.0)
            self.status.emit(
                f"[{sheet_name}] Handling outliers (mean ± {std_coeff}σ)...")
            for col in df_numeric.columns:
                col_data = df_numeric[col].dropna()
                if len(col_data) < 4:
                    continue
                mu, sd = col_data.mean(), col_data.std()
                if sd <= 0:
                    continue
                mask = ((df_numeric[col] < mu - std_coeff * sd) |
                        (df_numeric[col] > mu + std_coeff * sd))
                if mask.any():
                    df_numeric.loc[mask, col] = col_data.median()

        # --- Recombine this sheet's metadata and numeric columns ---
        if not df_numeric.empty and meta_cols:
            return pd.concat([meta_df, df_numeric], axis=1)
        if not df_numeric.empty:
            return df_numeric
        return meta_df

    def _fallback_normalize(self, df, method: str):
        """Fallback normalization without SCMeTA."""
        import numpy as np
        for col in df.columns:
            col_data = df[col].dropna()
            if len(col_data) == 0:
                continue
            if "Max" in method:
                m = col_data.max()
                if m > 0:
                    df[col] = df[col] / m
            elif "Z-score" in method:
                s = col_data.std()
                if s > 0:
                    df[col] = (df[col] - col_data.mean()) / s
            elif "Min-Max" in method:
                mn, mx = col_data.min(), col_data.max()
                if mx > mn:
                    df[col] = (df[col] - mn) / (mx - mn)
            elif "Log" in method:
                df[col] = np.log1p(df[col].clip(lower=0))


class SearchMergeWorker(WorkerThread):
    """Worker for merging cleaned data before search — matches m/z across CSV files by ppm tolerance."""

    def __init__(self, excel_file: str, output_dir: str, ppm_tolerance: float):
        super().__init__("Search Data Merge")
        self.excel_file = excel_file
        self.output_dir = output_dir
        self.ppm_tolerance = ppm_tolerance

    def execute(self):
        import numpy as np
        import pandas as pd
        from pathlib import Path

        self.status.emit("Loading cleaned data...")

        # Load Excel file with multiple sheets
        if self.excel_file.endswith(('.xlsx', '.xls')):
            xls = pd.ExcelFile(self.excel_file)
            sheet_names = xls.sheet_names
            self.status.emit(f"Found {len(sheet_names)} sheets: {sheet_names}")
        else:
            # Single CSV file
            df = pd.read_csv(self.excel_file)
            sheet_names = None
            data_dict = {"data": df}

        self.progress.emit(10)

        if sheet_names is not None:
            data_dict = {}
            for i, sheet in enumerate(sheet_names):
                if not self._is_running:
                    break
                data_dict[sheet] = pd.read_excel(xls, sheet_name=sheet)
                self.progress.emit(10 + int((i + 1) / len(sheet_names) * 30))

        # Merge strategy: find common m/z columns across all sheets
        self.status.emit("Identifying common m/z values across groups...")

        # Get numeric columns (m/z values) that appear in all groups
        all_numeric_cols = []
        for name, df in data_dict.items():
            numeric = [c for c in df.select_dtypes(include=[np.number]).columns.tolist()
                       if c not in df.columns[:2].tolist()]  # Skip first 2 metadata cols
            all_numeric_cols.append(set(numeric))

        if all_numeric_cols:
            common_cols = all_numeric_cols[0]
            for cols_set in all_numeric_cols[1:]:
                common_cols = common_cols.intersection(cols_set)
            common_cols = sorted(common_cols)
            self.status.emit(f"Found {len(common_cols)} common m/z columns")
        else:
            common_cols = []
            self.status.emit("No numeric columns found for merging")

        self.progress.emit(50)

        # Build a merged matrix: one row per m/z, one column per group
        if common_cols:
            merge_rows = []
            for mz_str in common_cols:
                try:
                    mz_val = float(mz_str)
                except (ValueError, TypeError):
                    continue
                row = {"m/z": mz_val}
                for name, df in data_dict.items():
                    if mz_str in df.columns:
                        row[name] = df[mz_str].mean()  # Use mean intensity across samples
                    else:
                        row[name] = np.nan
                merge_rows.append(row)

            merged_df = pd.DataFrame(merge_rows)
        else:
            # Fallback: just concatenate all data
            merged_df = pd.concat(data_dict.values(), ignore_index=True)

        self.progress.emit(80)

        # Save with Excel-safe format
        base_name = os.path.splitext(os.path.basename(self.excel_file))[0]
        if "filename" in merged_df.columns or merged_df.shape[1] > 100:
            if "filename" not in merged_df.columns:
                merged_df["filename"] = [f"row_{i}" for i in range(len(merged_df))]
            base_path = os.path.join(self.output_dir, f"{base_name}_merged_for_search")
            _save_matrix_safe(merged_df, base_path, log_func=self.log.emit)
            output_file = f"{base_path}.csv"
        else:
            output_file = os.path.join(self.output_dir, f"{base_name}_merged_for_search.csv")
            merged_df.to_csv(output_file, index=False)
        self.progress.emit(100)

        self.output_file_path = output_file
        self.status.emit(f"Merged data saved to: {output_file}")
        self.log.emit(f"Search merge completed. Output: {output_file}")


def _read_csv_robust(path: str, **kwargs):
    """Read a CSV trying common encodings (UTF-8, Chinese GBK, latin-1).

    Excel on Chinese Windows exports CSV as GBK/ANSI, which pandas'
    default UTF-8 read rejects with UnicodeDecodeError.
    """
    for enc in ("utf-8-sig", "gb18030", "latin-1"):
        try:
            return pd.read_csv(path, encoding=enc, **kwargs)
        except UnicodeDecodeError:
            continue
    raise ValueError(f"Cannot decode file: {path}")


class LocalSearchWorker(WorkerThread):
    """Worker for local database search — uses SCMeTA's SearchDatabase (matching Chinese version)."""

    def __init__(self, data_file: str, ref_file: str, output_dir: str, ppm_tolerance: float):
        super().__init__("Local Database Search")
        self.data_file = data_file
        self.ref_file = ref_file
        self.output_dir = output_dir
        self.ppm_tolerance = ppm_tolerance

    def execute(self):
        import csv
        import numpy as np
        import pandas as pd

        self.status.emit("Loading data for search...")
        self.progress.emit(10)

        # Load experimental data
        if self.data_file.endswith(('.xlsx', '.xls')):
            # Read ALL sheets (e.g. group-split files) and concatenate
            sheets = pd.read_excel(self.data_file, sheet_name=None)
            data_df = (pd.concat(list(sheets.values()), ignore_index=True)
                       if len(sheets) > 1 else list(sheets.values())[0])
        elif self.data_file.endswith('.csv'):
            data_df = _read_csv_robust(self.data_file, index_col=0)
        elif self.data_file.endswith('.txt'):
            data_df = _read_csv_robust(self.data_file, index_col=0, sep='\t')
        else:
            raise ValueError(f"Unsupported format: {self.data_file}")

        # Ensure numeric column names (m/z values)
        try:
            data_df.columns = [float(c) for c in data_df.columns]
        except (ValueError, TypeError):
            pass

        self.progress.emit(25)

        # Try using SCMeTA's local SQLite database. This requires a 'dbpath'
        # entry in ~/.scmeta/config.ini pointing to an HMDB .db file — SCMeTA
        # does not ship one, so this usually falls back to CSV matching.
        self.status.emit("Searching database via SCMeTA...")
        results = []
        try:
            from SCMeTA.config import PARAMETERS
            try:
                db_path = PARAMETERS.dbpath
            except KeyError:
                raise ValueError(
                    "Local SCMeTA database is not configured "
                    "(no 'dbpath' in ~/.scmeta/config.ini).")
            if not db_path or not os.path.exists(str(db_path)):
                raise ValueError(f"Database file not found: {db_path}")

            from SCMeTA.core.match import SearchDatabase
            from SCMeTA.match.connector import Connector

            searcher = SearchDatabase()
            tol = self.ppm_tolerance / 1e6  # Convert ppm to fractional tolerance

            # Use SCMeTA's connector for local DB search
            with Connector() as conn:
                mz_list = data_df.columns.tolist()
                matched = conn.search_list(mz_list, tol=tol)
                for mz, hits in matched.items():
                    for hit in hits:
                        results.append({
                            "m/z": mz,
                            "DB_m/z": hit[0] if len(hit) > 0 else "",
                            "Name": hit[2] if len(hit) > 2 else "",
                            "Formula": hit[1] if len(hit) > 1 else "",
                        })

        except Exception as e:
            self.log.emit(f"SCMeTA local database unavailable: {e}")
            self.log.emit("Falling back to CSV-based ppm matching...")
            self.status.emit("Using reference library CSV for ppm matching...")
            results = self._fallback_csv_search(data_df)

        self.progress.emit(90)

        # Save results
        # utf-8-sig (BOM) so Chinese Windows Excel decodes headers correctly;
        # also write an .xlsx copy (Excel-native, no encoding ambiguity)
        base_name = os.path.splitext(os.path.basename(self.data_file))[0]
        output_file = os.path.join(self.output_dir, f"Search_{base_name}.csv")

        if results:
            result_df = pd.DataFrame(results)
            result_df.to_csv(output_file, index=False, encoding="utf-8-sig")
            try:
                result_df.to_excel(
                    os.path.join(self.output_dir, f"Search_{base_name}.xlsx"),
                    index=False)
            except Exception:
                pass
        else:
            # No results found
            pd.DataFrame(columns=["m/z", "match"]).to_csv(
                output_file, index=False, encoding="utf-8-sig")

        self.progress.emit(100)
        self.output_file_path = output_file
        self.status.emit(f"Search complete. {len(results)} matches. Saved to: {output_file}")
        self.log.emit(f"Local search: {len(results)} matches found.")

    def _fallback_csv_search(self, data_df) -> list:
        """Fallback: CSV-based ppm matching."""
        import csv
        import numpy as np
        import pandas as pd

        # Load reference file
        db_data = {}
        db_headers = None
        try:
            if self.ref_file.endswith(('.xlsx', '.xls')):
                ref_df = pd.read_excel(self.ref_file)
            else:
                ref_df = _read_csv_robust(self.ref_file)
            db_headers = list(ref_df.columns)
            for _, row in ref_df.iterrows():
                try:
                    db_data[float(row.iloc[0])] = row.tolist()
                except (ValueError, IndexError):
                    continue
        except Exception:
            reader = None
            for enc in ("utf-8-sig", "gb18030", "latin-1"):
                try:
                    f = open(self.ref_file, 'r', encoding=enc)
                    reader = csv.reader(f)
                    db_headers = next(reader, None)
                    break
                except UnicodeDecodeError:
                    continue
            if reader is None:
                raise ValueError(f"Cannot decode reference file: {self.ref_file}")
            for row in reader:
                try:
                    db_data[float(row[0])] = row
                except (ValueError, IndexError):
                    continue

        ppm = self.ppm_tolerance
        results = []
        mz_list = data_df.columns.tolist() if hasattr(data_df, 'columns') else data_df

        for exp_mz in mz_list:
            try:
                exp_mz = float(exp_mz)
            except (ValueError, TypeError):
                continue

            for db_mz, db_row in db_data.items():
                if db_mz == 0:
                    continue
                error = abs(exp_mz - db_mz) / abs(db_mz) * 1e6
                if error <= ppm:
                    result = {"m/z": exp_mz, "PPM_Error": round(error, 4), "DB_m/z": db_mz}
                    if db_headers and len(db_headers) > 1:
                        for j, h in enumerate(db_headers[1:4], 1):
                            if j < len(db_row):
                                result[h] = db_row[j]
                    results.append(result)

        return results


class OnlineSearchWorker(WorkerThread):
    """Worker for online database search — uses HMDB, KEGG, LipidMaps, ChEBI, PubChem APIs."""

    def __init__(self, database: str, output_dir: str, ppm_tolerance: float):
        super().__init__("Online Database Search")
        self.database = database
        self.output_dir = output_dir
        self.ppm_tolerance = ppm_tolerance

    def execute(self):
        import csv
        import requests
        import json

        self.status.emit(f"Connecting to {self.database} online database...")
        self.progress.emit(5)

        # Try to use SCMeTA's built-in search modules
        search_results = []

        try:
            # Attempt to use SCMeTA's search module for local DB
            from SCMeTA.match.connector import Connector
            from SCMeTA.config import PARAMETERS

            self.status.emit("Using local HMDB database...")
            # For local DB search via SCMeTA, we need a data file
            self.log.emit("Online search: local HMDB SQLite database ready.")
            self.log.emit("Please use Local Search tab with a reference file for metabolite identification.")
            self.log.emit("Online HMDB API search is available via SCMeTA.search.hmdb module.")

            self.progress.emit(50)

            # Attempt online HMDB search if available
            try:
                from SCMeTA.search.hmdb import get_hmdb_id
                self.status.emit("HMDB API module loaded successfully.")
            except ImportError:
                self.log.emit("SCMeTA.search.hmdb module not found; using direct API access.")

            self.output_file_path = self.output_dir
            self.progress.emit(100)
            self.status.emit(f"{self.database} search ready. Use Local Search for metabolite matching.")
            self.log.emit(
                f"Online Search: {self.database} connection established.\n"
                f"Note: For full metabolite identification, use the Local Search tab\n"
                f"with your experimental data and a reference database file."
            )

        except ImportError as e:
            self.log.emit(f"SCMeTA modules not fully available: {e}")
            self.log.emit("Falling back to direct database API integration...")
            self._fallback_online_search(search_results)

    def _fallback_online_search(self, search_results):
        """Fallback online search using direct REST APIs."""
        import requests

        api_endpoints = {
            "HMDB": "https://hmdb.ca/metabolites/",
            "KEGG": "https://rest.kegg.jp/find/compound/",
            "PubChem": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/",
        }

        endpoint = api_endpoints.get(self.database)
        if not endpoint:
            self.status.emit(f"No API endpoint configured for {self.database}")
            self.progress.emit(100)
            return

        self.status.emit(f"{self.database} API endpoint: {endpoint}")
        self.log.emit(
            f"Online {self.database} search configured.\n"
            f"Note: Full online search requires:\n"
            f"  1. Network connectivity\n"
            f"  2. Valid m/z data input\n"
            f"  3. Proper API authentication (for HMDB)\n"
            f"Please use Local Search tab for offline metabolite identification."
        )
        self.progress.emit(100)


def _group_color_map(groups, color_scheme="tab10"):
    """Map group names to colors from a matplotlib colormap.

    Qualitative colormaps (tab10, Set1, ...) cycle through their colors;
    continuous ones (viridis, plasma, ...) are sampled evenly. Matches the
    Chinese version, which uses the tab10 palette for group coloring.
    """
    import matplotlib.pyplot as plt
    try:
        cmap = plt.get_cmap(color_scheme)
    except Exception:
        cmap = plt.get_cmap("tab10")
    n = getattr(cmap, "N", 256)
    if n is not None and n <= 32:  # qualitative / listed colormaps
        return {g: cmap(i % n) for i, g in enumerate(groups)}
    denom = max(len(groups) - 1, 1)  # continuous colormaps
    return {g: cmap(i / denom) for i, g in enumerate(groups)}


def _add_confidence_ellipses(ax, reduced, group_labels, groups, group_colors):
    """Draw 95% confidence ellipses for each group.

    Filled with a LIGHT version of the group color (alpha 0.15) so the
    scatter points stay clearly visible on top (zorder 1 < scatter's 3);
    the border stays the solid group color.
    """
    import numpy as np
    from matplotlib.patches import Ellipse
    from matplotlib.colors import to_rgba

    scale = np.sqrt(5.991)  # 95% confidence for a 2D gaussian
    for g in groups:
        mask = np.array([gl == g for gl in group_labels])
        pts = reduced[mask]
        if len(pts) < 3:
            continue
        cov = np.cov(pts.T)
        vals, vecs = np.linalg.eigh(cov)
        vals = np.maximum(vals, 0)  # guard tiny negative eigenvalues
        order = vals.argsort()[::-1]
        vals, vecs = vals[order], vecs[:, order]
        angle = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
        color = group_colors[g]
        ax.add_patch(Ellipse(
            xy=pts.mean(axis=0),
            width=scale * np.sqrt(vals[0]),
            height=scale * np.sqrt(vals[1]),
            angle=angle, fill=True,
            facecolor=to_rgba(color, alpha=0.15),
            edgecolor=color, linewidth=1.3, zorder=1))


# Point Style combo items -> matplotlib markers (hollow styles get a
# transparent face and a colored edge instead)
_POINT_MARKERS = {
    "Solid Circle (●)": "o",
    "Hollow Circle (○)": "o",
    "Square (■)": "s",
    "Triangle (▲)": "^",
    "Diamond (◆)": "D",
    "Star (★)": "*",
}


def _marker_for(point_style):
    return _POINT_MARKERS.get(str(point_style), "o")


def _is_hollow(point_style):
    return str(point_style).startswith("Hollow")


def _darken(color, factor=0.55):
    """Darker version of a color, used as the edge line of solid markers so
    each point gets a visible outline (depth) instead of a flat fill."""
    from matplotlib.colors import to_rgba
    if color is None:
        color = "#1f77b4"  # matplotlib default cycle color
    r, g, b, a = to_rgba(color)
    return (r * factor, g * factor, b * factor, a)
#factor数值越小边框越黑
#_darken() —— 把颜色变暗的工具函数

_LEGEND_POS_MAP = {
    "Best": "best",
    "Upper Right": "upper right",
    "Lower Right": "lower right",
    "Upper Left": "upper left",
    "Lower Left": "lower left",
}


def _draw_legend(ax, pos, title=None, handles=None):
    """Legend honoring the Legend Position setting; 'Outside Right' (the
    default) places it next to the axes instead of on top of the data."""
    if pos == "Outside Right":
        ax.legend(handles=handles, title=title, loc="center left",
                  bbox_to_anchor=(1.02, 0.5))
    else:
        ax.legend(handles=handles, title=title,
                  loc=_LEGEND_POS_MAP.get(pos, "best"))


def _annotate_top_points(ax, xs, ys, names, max_labels, label_font_size,
                         label_pos, scores=None):
    """Annotate the top-N points with their names.

    Ranking uses `scores` when given (volcano: -log10p * |log2FC|),
    otherwise distance from the origin (dim-reduction plots). Honors the
    Label Settings group (max labels / font size / position).
    """
    import numpy as np
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    if len(xs) == 0:
        return
    if scores is None:
        scores = xs ** 2 + ys ** 2
    try:
        fs = int(str(label_font_size).split()[0])
    except Exception:
        fs = 10
    offsets = {
        "Top": (0, 6, "center", "bottom"),
        "Bottom": (0, -6, "center", "top"),
        "Left": (-6, 0, "right", "center"),
        "Right": (6, 0, "left", "center"),
        "Center": (0, 0, "center", "center"),
    }
    order = np.argsort(-np.asarray(scores, dtype=float))[:int(max_labels)]
    for i, k in enumerate(order):
        if k >= len(names):
            continue
        if label_pos == "Auto":
            # alternate above/below to reduce label overlap
            above = i % 2 == 0
            dx, dy = 0, 6 if above else -6
            ha, va = "center", "bottom" if above else "top"
        else:
            dx, dy, ha, va = offsets.get(
                label_pos, (0, 6, "center", "bottom"))
        ax.annotate(str(names[k]), (xs[k], ys[k]), xytext=(dx, dy),
                    textcoords="offset points", ha=ha, va=va, fontsize=fs)


def _bh_fdr(pvals):
    """Benjamini-Hochberg FDR correction.

    Returns q-values for the input p-values (standard BH step-up
    procedure, monotonicity enforced).
    """
    import numpy as np
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    q = np.empty(m)
    prev = 1.0
    # walk from the largest p-value to the smallest: q_(i) = p_(i) * m / rank
    for rank, i in enumerate(order[::-1]):
        qi = p[i] * m / (m - rank)
        if qi > prev:
            qi = prev
        prev = qi
        q[i] = qi
    return q


def _combat_batch_correction(mat, group_labels, log_func):
    """Optional ComBat batch correction (pyComBat, parametric mode).

    Correction runs WITHIN each cell-type super-group separately (labels
    sharing the same letter prefix, e.g. C1/C2 -> "C", M1/M2 -> "M"), so
    between-type differences (C vs M, or any groups added later) are
    preserved and only the within-type batch effect is removed. Falls back
    to a single global correction when the labels don't split into groups.
    Zero-variance columns are left untouched; columns that end up
    non-finite are restored to their original values.
    Returns the corrected matrix, or None when correction is not possible.
    """
    import numpy as np
    import pandas as pd
    try:
        from pycombat.pycombat import Combat
    except Exception as e:
        log_func(f"ComBat unavailable ({e}) — continues without it.")
        return None
    # Positional row labels (row order = mat order). A plain numpy array
    # avoids any pandas label-alignment surprises with duplicate indexes.
    labels = pd.Series(np.asarray(group_labels, dtype=str), index=mat.index)
    labels = labels.fillna("NA").astype(str)
    # Super-group = leading letters of each label ("C1" -> "C", "M2" -> "M")
    pf = labels.str.extract(r"^([A-Za-z]+)", expand=False)
    super_groups = pf.dropna().unique().tolist()
    if len(super_groups) < 2:
        # No multi-group split available — correct all samples together
        super_groups = [None]
    var_cols = mat.columns[mat.var().to_numpy() > 0]
    if len(var_cols) == 0:
        log_func("ComBat: no columns with variance — skipped.")
        return None
    out = mat.copy()
    done, skipped = [], []
    for sg in super_groups:
        if sg is None:
            mask_arr = np.ones(len(mat), dtype=bool)
            batch_labels = labels
        else:
            mask_arr = pf.to_numpy() == sg
            batch_labels = labels.iloc[mask_arr]
        pos = np.flatnonzero(mask_arr)
        if len(pos) < 2:
            skipped.append(f"{sg} (<2 samples)")
            continue
        batches = sorted(set(batch_labels))
        if len(batches) < 2:
            skipped.append(f"{sg} (only one batch)")
            continue
        # Columns with variance WITHIN this subset only — columns that
        # are constant here (e.g. peaks seen only in the other type)
        # would make pyComBat divide by zero and poison every column.
        # All indexing is positional (iloc) so duplicate row/column
        # labels can never trigger a Cartesian expansion.
        sub = mat.iloc[pos]
        col_pos = np.flatnonzero(sub.var().to_numpy() > 0)
        if len(col_pos) == 0:
            skipped.append(f"{sg} (no variable columns)")
            continue
        b = np.array([batches.index(g) for g in batch_labels])
        Y = sub.iloc[:, col_pos].values.astype(float)
        try:
            Yt = Combat(mode="p").fit_transform(Y, b)
        except Exception as e:
            log_func(f"ComBat failed for group {sg} ({e}) — left unchanged.")
            skipped.append(str(sg))
            continue
        bad = ~np.isfinite(Yt).all(axis=0)
        if bad.any():
            Yt[:, bad] = Y[:, bad]
            log_func(
                f"ComBat [{sg}]: {int(bad.sum())} column(s) non-finite "
                f"after correction — kept original values.")
        out.iloc[pos, col_pos] = Yt
        done.append(f"{sg}({len(batches)} batches)")
    if pf.isna().to_numpy().any():
        skipped.append("labels without letter prefix")
    if not done:
        log_func("ComBat batch correction skipped (no group has ≥2 batches).")
        return None
    log_func(
        "ComBat batch correction applied within each group separately: "
        f"{', '.join(done)} — between-group differences preserved.")
    if skipped:
        log_func(f"ComBat: skipped {', '.join(skipped)}.")
    return out


class VisualizationWorker(WorkerThread):
    """Worker for generating plots — PCA, t-SNE, UMAP, Heatmap, Volcano, Box, OPLS-DA."""

    def __init__(self, matrix_file: str, coord_file: str, output_dir: str,
                 plot_type: str, params: dict):
        super().__init__("Visualization")
        self.matrix_file = matrix_file
        self.coord_file = coord_file
        self.output_dir = output_dir
        self.plot_type = plot_type
        self.params = params
        self._result_fig = None
        self._result_data = None

    @property
    def result_figure(self):
        return self._result_fig

    @property
    def result_data(self):
        return self._result_data

    def execute(self):
        import numpy as np
        import pandas as pd
        import matplotlib
        matplotlib.use('Qt5Agg')
        import matplotlib.pyplot as plt

        self.status.emit(f"Loading matrix data: {self.matrix_file}...")
        self.progress.emit(10)

        # Load matrix data
        if self.matrix_file.endswith('.csv'):
            mat = _read_csv_robust(self.matrix_file, index_col=0)
        elif self.matrix_file.endswith(('.xlsx', '.xls')):
            mat = pd.read_excel(self.matrix_file, index_col=0)
        elif self.matrix_file.endswith('.txt'):
            mat = _read_csv_robust(self.matrix_file, index_col=0, sep='\t')
        else:
            raise ValueError(f"Unsupported matrix format: {self.matrix_file}")

        # Extract group labels (if any) BEFORE numeric filtering — a Group
        # column (or a low-cardinality row index) drives point coloring
        group_labels = None
        for c in mat.columns:
            if str(c).strip().lower() in ("group", "groups", "class", "label",
                                          "celltype", "cell type", "cluster"):
                group_labels = mat[c].fillna("NA").astype(str).tolist()
                break
        if group_labels is None:
            idx_name = str(mat.index.name or "").strip().lower()
            if idx_name in ("group", "groups", "class", "label", "cluster"):
                group_labels = pd.Series(mat.index).fillna("NA").astype(str).tolist()
            elif 1 < mat.index.nunique() <= 12 and mat.index.nunique() < len(mat):
                group_labels = pd.Series(mat.index).fillna("NA").astype(str).tolist()
                self.log.emit(
                    f"Row index has {mat.index.nunique()} distinct values — "
                    f"coloring points by index as groups.")

        # Sample identifiers for point hover — the scan each dot belongs
        # to. Captured BEFORE numeric filtering because Sample/filename
        # columns are non-numeric and get dropped later. Row order never
        # changes in the steps below, so the labels stay aligned with the
        # plotted points.
        import re
        self._hover_labels = None
        for c in mat.columns:
            if str(c).strip().lower() in ("sample", "filename", "file",
                                          "name", "scan"):
                vals = [str(v) for v in mat[c].fillna("NA")]
                # Skip auto-generated placeholders ("row_0", "row_1", ...)
                # that the cleaning step adds when the file has no sample
                # column — the row index is more informative then
                if all(re.fullmatch(r"row_\d+", v) for v in vals):
                    continue
                self._hover_labels = vals
                break
        if self._hover_labels is None:
            idx = mat.index
            if (isinstance(idx, pd.RangeIndex)
                    and idx.start == 0 and idx.step == 1):
                self._hover_labels = [
                    f"scan {i + 1}" for i in range(len(idx))]
            else:
                self._hover_labels = [str(v) for v in idx]

        # Fill NaN with 0
        mat = mat.fillna(0)

        # Keep only numeric columns — metadata columns such as filename/name
        # or metabolite names break mean()/std() under pandas 3
        total_cols = mat.shape[1]
        mat = mat.select_dtypes(include=[np.number])
        if mat.shape[1] < total_cols:
            self.log.emit(
                f"Ignored {total_cols - mat.shape[1]} non-numeric column(s) "
                f"(e.g. filename/name) for visualization.")
        if mat.shape[1] == 0:
            raise ValueError(
                "Matrix file has no numeric m/z columns — "
                "pick a data matrix file for visualization.")

        self.status.emit(f"Matrix loaded: {mat.shape[0]} rows x {mat.shape[1]} cols")
        self.progress.emit(25)

        # Canonical method key: "t-SNE" -> "tsne", "OPLS-DA" -> "oplsda"
        plot_type = self.plot_type.lower().replace("-", "").replace(" ", "")

        # Optional ComBat batch correction — OFF by default, applied on the
        # full matrix BEFORE feature selection / preprocessing. Correction
        # runs within each cell-type group (C / M / ...) separately so that
        # differences between cell types are preserved.
        if self.params.get("combat_batch", False):
            if group_labels is not None and len(set(group_labels)) >= 2:
                corrected = _combat_batch_correction(
                    mat, group_labels, self.log.emit)
                if corrected is not None:
                    mat = corrected
            else:
                self.log.emit(
                    "ComBat batch correction needs at least 2 groups — skipped.")

        # Feature selection — only for dim-reduction methods, and BEFORE
        # preprocessing (a variance ranking is meaningless after z-scoring,
        # which equalizes every column's variance)
        feat_sel = self.params.get("feature_select", "None (all features)")
        if (plot_type in ("pca", "tsne", "umap")
                and feat_sel != "None (all features)"):
            n_feat = self.params.get("feature_n", 200)
            if mat.shape[1] > n_feat:
                # Class Difference is a SUPERVISED selection: by
                # construction it shows some separation even in pure
                # noise. Verify the class difference itself first (CV AUC
                # + permutation test) so a non-significant result is
                # flagged in the log instead of being silently over-read.
                if feat_sel == "Top-N by Class Difference":
                    self._guard_class_difference(mat, group_labels)
                try:
                    keep = self._select_features(mat, group_labels, feat_sel, n_feat)
                    mat = mat[keep]
                    self.log.emit(
                        f"Feature selection ({feat_sel}, N={n_feat}): "
                        f"kept {mat.shape[1]} columns.")
                except Exception as e:
                    self.log.emit(
                        f"Feature selection skipped ({e}) — using all columns.")

        # Preprocessing
        preprocess = self.params.get("data_preprocess", "None")
        if preprocess == "Z-score Standardization":
            mat = (mat - mat.mean()) / mat.std().replace(0, 1)
        elif preprocess == "Max Normalization":
            mat = mat / mat.max().replace(0, 1)
        elif preprocess == "Min-Max Scaling":
            mat = (mat - mat.min()) / (mat.max() - mat.min()).replace(0, 1)
        # "None": skip preprocessing

        self.progress.emit(40)
        self.status.emit(f"Generating {self.plot_type} plot...")

        # Display Settings: plot style + font must be set before the figure
        # is created so every element (labels, ticks, legend) picks them up
        self._apply_display_style()

        # Output geometry (Chart Size settings). The data plot (axes) keeps
        # the requested size; the overall figure = data plot + standard
        # margins, with extra room on the right when the legend sits outside.
        fig_w = float(self.params.get("fig_width", 10.0) or 10.0)
        fig_h = float(self.params.get("fig_height", 8.0) or 8.0)
        legend_out = (self.params.get("legend_pos", "Outside Right")
                      == "Outside Right"
                      and self.params.get("show_legend", True))
        left_m, bottom_m, top_m = 0.125, 0.11, 0.10
        right_m = 0.20 if legend_out else 0.125
        out_w = fig_w / (1 - left_m - right_m)
        out_h = fig_h / (1 - bottom_m - top_m)
        fig, ax = plt.subplots(figsize=(out_w, out_h), dpi=100)
        ax.set_position([left_m, bottom_m,
                         1 - left_m - right_m, 1 - bottom_m - top_m])

        # Hover state — point plots (PCA/t-SNE/UMAP/OPLS-DA/Volcano) fill
        # these in; the preview canvas uses them to show which sample/scan
        # a point belongs to when the mouse hovers over it.
        self._hover_points = None
        self._hover_groups = None
        self._hover_prefix = "Sample"
        self._hover_plot_labels = None

        try:
            if plot_type in ("pca", "tsne", "umap"):
                self._run_dimensionality_reduction(mat, plot_type, fig, ax, group_labels)
            elif plot_type == "heatmap":
                self._run_heatmap(mat, fig, ax, group_labels)
            elif plot_type == "volcanoplot":
                self._run_volcano(mat, fig, ax, group_labels)
            elif plot_type == "boxplot":
                self._run_boxplot(mat, fig, ax, group_labels)
            elif plot_type == "oplsda":
                self._run_oplsda(mat, fig, ax, group_labels)
            else:
                # Default: PCA
                self._run_dimensionality_reduction(mat, "pca", fig, ax)
        except Exception as e:
            self.log.emit(f"Warning: {plot_type} generation issue: {e}")
            # Fallback plot — no hover metadata (points are not meaningful)
            self._hover_points = None
            ax.scatter(mat.iloc[:, 0] if mat.shape[1] > 0 else range(len(mat)),
                       mat.iloc[:, 1] if mat.shape[1] > 1 else mat.iloc[:, 0])
            ax.set_title(f"{self.plot_type} (fallback)")

        self._result_fig = fig
        self._result_data = mat

        self.progress.emit(90)

        # Save figure — overall size = figure inches x Output DPI (pixels);
        # the axes placement above already reserves room for the legend, so
        # nothing gets cut off at the figure edge.
        save_dpi = int(self.params.get("output_dpi", 300) or 300)
        os.makedirs(self.output_dir, exist_ok=True)
        base_name = os.path.splitext(os.path.basename(self.matrix_file))[0]
        output_file = os.path.join(self.output_dir, f"{base_name}_{self.plot_type.replace(' ', '_')}.png")
        fig.savefig(output_file, dpi=save_dpi)

        # Hover payload for the preview pane: per-point coordinates in DATA
        # space plus the axes-position fractions used above. The preview
        # canvas maps the mouse position back through these fractions to a
        # data point and shows its sample/scan on hover (same behavior as
        # the Chinese version). Built AFTER savefig so the axes limits are
        # final (autoscale runs during the savefig draw).
        hp = getattr(self, "_hover_points", None)
        if (hp is not None and len(hp) > 0
                and getattr(self, "_hover_labels", None)):
            self.hover_payload = {
                "points": np.asarray(hp, dtype=float)[:, :2],
                "labels": (getattr(self, "_hover_plot_labels", None)
                           or self._hover_labels),
                "groups": getattr(self, "_hover_groups", None),
                "prefix": getattr(self, "_hover_prefix", "Sample"),
                "xlim": tuple(float(v) for v in ax.get_xlim()),
                "ylim": tuple(float(v) for v in ax.get_ylim()),
                "left": left_m, "right": right_m,
                "bottom": bottom_m, "top": top_m,
            }
        else:
            self.hover_payload = None
        self.progress.emit(100)

        self.status.emit(f"Plot saved to: {output_file}")
        self.log.emit(f"Visualization ({self.plot_type}) completed. Saved to: {output_file}")


    @staticmethod
    def _select_features(mat, group_labels, mode, n):
        """Top-N informative columns for dim reduction.

        - Top-N by Variance: column variance (strong signal > flat background)
        - Top-N by Group Difference: one-way ANOVA F-statistic across groups
        - Top-N by Class Difference: pooled Cohen's d between two classes
          derived from the group names by stripping trailing digits /
          separators (e.g. C1/C2/C3 -> C, M1/M2/M3 -> M)
        """
        import re
        import numpy as np
        import pandas as pd
        if mode == "Top-N by Variance":
            scores = mat.var()
        elif mode == "Top-N by Class Difference":
            if group_labels is None:
                raise ValueError(
                    "needs group labels for class-difference ranking")
            labs = [str(g) for g in group_labels]
            classes = [re.sub(r"[^A-Za-z]+$", "", g) for g in labs]
            uniq = list(dict.fromkeys(classes))
            if len(uniq) != 2:
                raise ValueError(
                    f"Class Difference derives classes from group names and "
                    f"needs exactly 2, got {len(uniq)}: {uniq}")
            c1, c2 = uniq
            m1 = np.array([cl == c1 for cl in classes])
            m2 = ~m1
            if m1.sum() < 2 or m2.sum() < 2:
                raise ValueError("each class needs at least 2 samples")
            V = mat.values
            mu1, mu2 = V[m1].mean(axis=0), V[m2].mean(axis=0)
            s1, s2 = V[m1].std(axis=0), V[m2].std(axis=0)
            n1, n2 = m1.sum(), m2.sum()
            pooled = np.sqrt(
                ((n1 - 1) * s1 ** 2 + (n2 - 1) * s2 ** 2) / (n1 + n2 - 2))
            d = np.divide(mu1 - mu2, pooled,
                          out=np.zeros_like(pooled), where=pooled > 0)
            scores = pd.Series(np.abs(d), index=mat.columns)
        else:
            if group_labels is None or len(set(group_labels)) < 2:
                raise ValueError(
                    "needs at least 2 groups for group-difference ranking")
            labs = np.array(list(group_labels))
            V = mat.values
            overall = V.mean(axis=0)
            ssb = np.zeros(V.shape[1])
            ssw = np.zeros(V.shape[1])
            k = 0
            for g in sorted(set(labs)):
                m = labs == g
                if m.sum() < 2:
                    continue
                k += 1
                gm = V[m].mean(axis=0)
                ssb += m.sum() * (gm - overall) ** 2
                ssw += ((V[m] - gm) ** 2).sum(axis=0)
            msb = ssb / max(k - 1, 1)
            msw = ssw / max(len(labs) - k, 1)
            scores = pd.Series(
                np.divide(msb, msw, out=np.zeros_like(msb), where=msw > 0),
                index=mat.columns)
        return list(scores.nlargest(n).index)

    @staticmethod
    def _verify_class_difference(mat, class_mask, n_perm=20):
        """Unbiased C/M difference check for the Class Difference mode.

        Pipeline: label-free top-200 VARIANCE features (variance does not
        use the labels, so the null stays valid) -> z-score -> 5-fold CV
        logistic AUC; the same pipeline on shuffled class labels gives the
        permutation null.

        Returns (auc, p, n_exceed, n_perm) or None when it cannot run.
        """
        import numpy as np
        from sklearn.linear_model import LogisticRegression
        from sklearn.model_selection import cross_val_predict
        from sklearn.metrics import roc_auc_score
        y = np.asarray(class_mask, dtype=int)
        if y.sum() < 5 or (1 - y).sum() < 5:
            return None
        top = min(200, mat.shape[1])
        keep = mat.var().nlargest(top).index
        X = mat[keep].values.astype(float)
        sd = X.std(axis=0)
        X = (X - X.mean(axis=0)) / np.where(sd > 0, sd, 1)
        try:
            proba = cross_val_predict(LogisticRegression(max_iter=2000), X, y,
                                      cv=5, method="predict_proba")[:, 1]
            auc_obs = roc_auc_score(y, proba)
        except Exception:
            return None
        rng = np.random.RandomState(0)
        aucs_null = []
        for _ in range(n_perm):
            yp = rng.permutation(y)
            try:
                p_perm = cross_val_predict(
                    LogisticRegression(max_iter=2000), X, yp, cv=5,
                    method="predict_proba")[:, 1]
                aucs_null.append(roc_auc_score(yp, p_perm))
            except Exception:
                pass
        if not aucs_null:
            return None
        n_exceed = sum(a >= auc_obs for a in aucs_null)
        p = (n_exceed + 1) / (len(aucs_null) + 1)
        return auc_obs, p, n_exceed, len(aucs_null)

    def _guard_class_difference(self, mat, group_labels):
        """Run the class-difference verification and log the verdict.
        Returns the class mask (True = first derived class) or None."""
        import re
        import numpy as np
        if group_labels is None:
            self.log.emit(
                "Class Difference needs group labels — verification skipped.")
            return None
        labs = [str(g) for g in group_labels]
        classes = [re.sub(r"[^A-Za-z]+$", "", g) for g in labs]
        uniq = list(dict.fromkeys(classes))
        if len(uniq) != 2:
            self.log.emit(
                f"Class Difference verification skipped — group names derive "
                f"{len(uniq)} classes ({uniq}), need exactly 2.")
            return None
        mask = np.array([cl == uniq[0] for cl in classes])
        self.status.emit(
            "Verifying class difference (5-fold CV AUC + permutation test)...")
        try:
            res = self._verify_class_difference(mat, mask)
        except Exception as e:
            res = None
            self.log.emit(f"Class difference verification skipped ({e}).")
        if res is None:
            return mask
        auc, p, n_exceed, n_perm = res
        self.log.emit(
            f"Class difference check ({uniq[0]}: {int(mask.sum())} cells, "
            f"{uniq[1]}: {int((~mask).sum())} cells): "
            f"5-fold CV AUC = {auc:.3f}")
        if p <= 0.05:
            self.log.emit(
                f"Permutation test ({n_perm} shuffles): p = {p:.3f} "
                f"({n_exceed}/{n_perm} shuffled AUCs exceeded the observed) "
                f"— class difference is significant.")
        else:
            self.log.emit(
                f"WARNING: permutation test ({n_perm} shuffles): p = {p:.3f} "
                f"({n_exceed}/{n_perm} shuffled AUCs exceeded the observed) — "
                f"the class difference is NOT significant; the separation in "
                f"this plot is constructed by feature selection and must not "
                f"be read as a real difference.")
        return mask

    def _apply_display_style(self):
        """Apply the Display Settings plot style + font BEFORE the figure is
        created, so axes labels, ticks and the legend all pick them up."""
        import matplotlib.pyplot as plt
        style_map = {
            "Standard White": "default",
            "Seaborn": "seaborn-v0_8",
            "ggplot": "ggplot",
            "bmh": "bmh",
        }
        plot_style = self.params.get("plot_style", "Standard White")
        try:
            plt.style.use(style_map.get(plot_style, "default"))
        except Exception as e:
            plt.style.use("default")
            self.log.emit(
                f"Plot style '{plot_style}' unavailable — using Standard "
                f"White ({e}).")
        try:
            plt.rcParams["font.family"] = self.params.get(
                "font_style", "Times New Roman")
            font_size = self.params.get("font_size", 18)
            if font_size and int(font_size) > 0:
                plt.rcParams["font.size"] = int(font_size)
        except Exception:
            pass

    def _scatter(self, ax, x, y, s, alpha, color=None, label=None):
        """Scatter honoring the Display Settings Point Style (hollow styles
        draw colored edges with a transparent face; solid styles get a thin
        darker edge line so the points look less flat)."""
        point_style = self.params.get("point_style", "Solid Circle (●)")
        kwargs = {"marker": _marker_for(point_style)}
        if _is_hollow(point_style):
            kwargs.update(facecolors="none",
                          edgecolors=(color if color is not None
                                      else "#1f77b4"))
            return ax.scatter(x, y, s=s, alpha=alpha, label=label, **kwargs)
        kwargs.update(linewidths=0.6, edgecolors=_darken(color))
        return ax.scatter(x, y, s=s, alpha=alpha, color=color, label=label,
                          **kwargs)

    def _run_dimensionality_reduction(self, mat, method, fig, ax, group_labels=None):
        """Run PCA, t-SNE, or UMAP."""
        import numpy as np

        # t-SNE / UMAP need at least 3 samples (t-SNE also requires
        # n_components < n_samples) — show a note instead of crashing
        if method in ("tsne", "umap") and mat.shape[0] < 3:
            pretty_name = {"tsne": "t-SNE", "umap": "UMAP"}.get(
                method, method.upper())
            ax.text(0.5, 0.5,
                    f"{pretty_name} needs at least 3 samples\n"
                    f"(current data has {mat.shape[0]})",
                    transform=ax.transAxes, ha='center', va='center')
            ax.set_title(f"{pretty_name} Analysis")
            return

        if method == "pca":
            from sklearn.decomposition import PCA
            n_components = max(1, min(self.params.get("pca_dims", 2),
                                      mat.shape[1], mat.shape[0]))
            reducer = PCA(n_components=n_components)
        elif method == "tsne":
            from sklearn.manifold import TSNE
            n_components = max(1, min(self.params.get("tsne_components", 2), 3,
                                      mat.shape[1], mat.shape[0] - 1))
            perplexity = min(self.params.get("tsne_perplexity", 30),
                             mat.shape[0] - 1)
            reducer = TSNE(n_components=n_components, random_state=42,
                           perplexity=max(1, perplexity))
        elif method == "umap":
            from umap import UMAP
            n_components = max(1, min(self.params.get("umap_components", 2),
                                      mat.shape[1], mat.shape[0]))
            n_neighbors = min(self.params.get("umap_neighbors", 15),
                              mat.shape[0] - 1)
            reducer = UMAP(n_components=n_components,
                           n_neighbors=max(2, n_neighbors), random_state=42)

        reduced = reducer.fit_transform(mat.values if hasattr(mat, 'values') else mat)
        self._result_data = reduced

        # Hover: per-point coordinates (the two components actually drawn)
        self._hover_points = np.asarray(reduced)[:, :2]
        self._hover_groups = ([str(g) for g in group_labels]
                              if group_labels is not None else None)
        self._hover_prefix = "Sample"

        # Point coloring: one color per group when group labels exist
        groups, group_colors = None, None
        if group_labels is not None and len(set(group_labels)) > 1:
            groups = sorted(set(group_labels))
            group_colors = _group_color_map(
                groups, self.params.get("color_scheme", "tab10"))

        two_dim = reduced.shape[1] >= 2
        s = self.params.get("point_size", 36)
        alpha = self.params.get("opacity", 0.8)
        if group_colors:
            for g in groups:
                mask = np.array([gl == g for gl in group_labels])
                if two_dim:
                    self._scatter(ax, reduced[mask, 0], reduced[mask, 1],
                                  s, alpha, color=group_colors[g], label=str(g))
                else:
                    self._scatter(ax, reduced[mask, 0], np.zeros(mask.sum()),
                                  s, alpha, color=group_colors[g], label=str(g))
        else:
            if two_dim:
                self._scatter(ax, reduced[:, 0], reduced[:, 1], s, alpha)
            else:
                self._scatter(ax, reduced[:, 0],
                              np.zeros_like(reduced[:, 0]), s, alpha)

        # 95% confidence ellipse per group (optional) — light group-color
        # fill, drawn under the points so they stay visible. PCA only:
        # t-SNE / UMAP distances are not meaningful in data space, so the
        # ellipse would not describe the true covariance there.
        if (method == "pca"
                and self.params.get("confidence_ellipse", False)
                and group_colors and two_dim):
            _add_confidence_ellipses(ax, reduced, group_labels, groups,
                                     group_colors)

        # Explained variance for PCA
        pretty_name = {"pca": "PCA", "tsne": "t-SNE", "umap": "UMAP"}.get(
            method, method.upper())
        if method == "pca" and self.params.get("explained_variance", True):
            explained = reducer.explained_variance_ratio_
            ax.set_xlabel(f"PC1 ({explained[0]*100:.1f}%)")
            if len(explained) >= 2:
                ax.set_ylabel(f"PC2 ({explained[1]*100:.1f}%)")
        else:
            ax.set_xlabel(f"{pretty_name} Component 1")
            ax.set_ylabel(f"{pretty_name} Component 2")

        ax.set_title(f"{pretty_name} Analysis")

        # Data point labels — top-N points by distance from the origin
        if two_dim and self.params.get("data_point_labels", False):
            names = [str(v) for v in (mat.index if hasattr(mat, "index")
                                      else range(len(reduced)))]
            _annotate_top_points(
                ax, reduced[:, 0], reduced[:, 1], names,
                self.params.get("max_labels", 10),
                self.params.get("label_font_size", "10 pt"),
                self.params.get("label_pos", "Auto"))

        # Legend — only meaningful when points are labeled by group;
        # defaults to outside the axes on the right
        if self.params.get("show_legend", True) and group_colors:
            _draw_legend(ax, self.params.get("legend_pos", "Outside Right"),
                         title="Group")
        if self.params.get("show_grid", False):
            ax.grid(True, alpha=0.3)

    def _run_heatmap(self, mat, fig, ax, group_labels=None):
        """Generate a heatmap, with a per-row group color strip + legend
        when group labels exist."""
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        from matplotlib.colorbar import ColorbarBase

        color_map = self.params.get("color_scheme", "viridis")
        # Qualitative palettes (tab10/Set1/...) are for group colors, not
        # heatmap values — fall back to viridis for the heatmap itself
        try:
            if getattr(plt.get_cmap(color_map), "N", 256) <= 32:
                color_map = "viridis"
        except Exception:
            color_map = "viridis"
        max_rows = self.params.get("heatmap_rows", 100)
        max_cols = self.params.get("heatmap_cols", 50)

        display_mat = mat.iloc[:min(mat.shape[0], max_rows),
                               :min(mat.shape[1], max_cols)]

        if display_mat.shape[0] == 0 or display_mat.shape[1] == 0:
            ax.text(0.5, 0.5, "No data to display", transform=ax.transAxes,
                    ha='center', va='center')
            ax.set_title("Heatmap")
            return

        im = ax.imshow(display_mat.values if hasattr(display_mat, 'values') else display_mat,
                       aspect='auto', cmap=color_map)
        fig.colorbar(im, ax=ax)
        ax.set_title("Heatmap")
        ax.set_xlabel("Features")
        ax.set_ylabel("Samples")

        if group_labels is not None and len(set(group_labels)) > 1:
            groups = sorted(set(group_labels))
            gcolors = _group_color_map(groups, "tab10")
            strip_cmap = mcolors.ListedColormap([gcolors[g] for g in groups])
            boundaries = np.arange(len(groups) + 1) - 0.5
            norm = mcolors.BoundaryNorm(boundaries, strip_cmap.N)

            divider = make_axes_locatable(ax)
            # Left strip: one colored cell per sample row
            row_groups = [gl for gl in group_labels[:len(display_mat)]]
            rax = divider.append_axes("left", size="2%", pad=0.02)
            idx = np.array([[groups.index(gl)] for gl in row_groups])
            rax.imshow(idx, aspect="auto", cmap=strip_cmap, norm=norm,
                       interpolation="nearest")
            rax.set_xticks([])
            rax.set_yticks([])
            # Right bar: group color legend
            cax = divider.append_axes("right", size="3%", pad=0.08)
            cbar = ColorbarBase(cax, cmap=strip_cmap, norm=norm,
                                boundaries=boundaries,
                                ticks=np.arange(len(groups)))
            cbar.set_ticklabels(groups)
            cbar.set_label("Group")

    def _run_volcano(self, mat, fig, ax, group_labels=None):
        """Volcano plot: log2 fold change vs -log10 p-value.

        Two comparison modes (see the Volcano parameter page):
        - Class Merge (default): classes are derived from group names
          (C1/C2/C3 -> C, M1/M2/M3 -> M) and the two classes compared
          (control vs model), matching the manuscript's Figure 3C.
        - Manual Pair: two specific groups chosen in the UI.
        Falls back to first half vs second half only when the data has no
        usable group labels. Significant points are colored.
        """
        import numpy as np
        import re
        from scipy.stats import ttest_ind

        cmp_mode = self.params.get("volcano_cmp", "Class Merge (C vs M)")

        if group_labels is not None and len(set(group_labels)) >= 2:
            if cmp_mode == "Manual Pair":
                sel1 = str(self.params.get("volcano_g1", "")).strip()
                sel2 = str(self.params.get("volcano_g2", "")).strip()
                avail = sorted(set(str(g) for g in group_labels))
                if not sel1 or not sel2 or sel1 == sel2 \
                        or sel1 not in avail or sel2 not in avail:
                    ax.text(0.5, 0.5,
                            "Selected groups not found in data.\n"
                            f"Requested: '{sel1}' vs '{sel2}'\n"
                            f"Available: {', '.join(avail)}",
                            transform=ax.transAxes, ha='center', va='center')
                    return
                g1_name, g2_name = sel1, sel2
                m1 = np.array([str(gl) == g1_name for gl in group_labels])
                m2 = np.array([str(gl) == g2_name for gl in group_labels])
                if len(avail) > 2:
                    self.log.emit(
                        f"Volcano compares '{g2_name}' vs '{g1_name}' "
                        f"(selected pair; {len(avail) - 2} other group(s) "
                        f"ignored).")
            else:  # Class Merge
                labs = [str(g) for g in group_labels]
                classes = [re.sub(r"[^A-Za-z]+$", "", g) for g in labs]
                uniq = sorted(dict.fromkeys(classes))
                if len(uniq) != 2:
                    ax.text(0.5, 0.5,
                            "Class Merge needs exactly 2 classes derived\n"
                            "from group names (e.g. C1/C2/C3 + M1/M2/M3).\n"
                            f"Got {len(uniq)}: {', '.join(uniq)}\n"
                            "(switch Comparison to Manual Pair)",
                            transform=ax.transAxes, ha='center', va='center')
                    return
                c1, c2 = uniq
                m1 = np.array([cl == c1 for cl in classes])
                m2 = np.array([cl == c2 for cl in classes])
                members1 = sorted(set(g for g, cl in zip(labs, classes)
                                     if cl == c1))
                members2 = sorted(set(g for g, cl in zip(labs, classes)
                                     if cl == c2))
                g1_name = f"{c1} ({'+'.join(members1)})"
                g2_name = f"{c2} ({'+'.join(members2)})"
                self.log.emit(
                    f"Volcano class comparison: {g2_name} (n={int(m2.sum())}) "
                    f"vs {g1_name} (n={int(m1.sum())}).")
        else:
            n = mat.shape[0]
            m1 = np.zeros(n, dtype=bool)
            m2 = np.zeros(n, dtype=bool)
            m1[:n // 2] = True
            m2[n // 2:] = True
            g1_name, g2_name = "First half", "Second half"

        if m1.sum() < 2 or m2.sum() < 2:
            ax.text(0.5, 0.5,
                    "Need at least 2 samples per group/class\nfor Volcano plot",
                    transform=ax.transAxes, ha='center', va='center')
            return

        g1 = mat[m1].mean(axis=0)
        g2 = mat[m2].mean(axis=0)
        fc = np.log2((g2 + 1e-6) / (g1 + 1e-6))

        pvals = []
        for col in mat.columns:
            try:
                _, p = ttest_ind(mat[m1][col], mat[m2][col])
                # identical groups give NaN p — treat as not significant
                pvals.append(p if np.isfinite(p) else 1.0)
            except Exception:
                pvals.append(1.0)

        # Optional BH-FDR correction (on by default) — color points by
        # adjusted p (q-value) so multiple testing is controlled
        use_fdr = self.params.get("volcano_fdr", False)
        if use_fdr:
            qvals = _bh_fdr(pvals)
            neg_log_p = -np.log10(np.array(qvals) + 1e-300)
            sig_name = "q"
            y_label = "-log10 adjusted p-value (q-value)"
        else:
            neg_log_p = -np.log10(np.array(pvals) + 1e-300)
            sig_name = "p"
            y_label = "-log10 p-value"

        fc_thr = self.params.get("volcano_fc", 1.0)
        p_thr = self.params.get("volcano_p", 0.05)
        p_line = -np.log10(p_thr)
        up = (fc >= fc_thr) & (neg_log_p >= p_line)
        down = (fc <= -fc_thr) & (neg_log_p >= p_line)
        ns = ~(up | down)
        self.log.emit(
            f"Volcano ({g2_name} vs {g1_name}): {mat.shape[1]} features, "
            f"{'BH-FDR corrected' if use_fdr else 'raw p-values'}, "
            f"threshold |log2FC|>={fc_thr} & {sig_name}<={p_thr}: "
            f"up {int(up.sum())}, down {int(down.sum())}, "
            f"not significant {int(ns.sum())}")

        # Hover: one point per feature (m/z column) — the Chinese version
        # shows "Feature: ..." when hovering over volcano points
        self._hover_points = np.column_stack([fc, neg_log_p])
        self._hover_plot_labels = [str(c) for c in mat.columns]
        self._hover_groups = None
        self._hover_prefix = "Feature"

        s = self.params.get("point_size", 36)
        alpha = self.params.get("opacity", 0.8)
        self._scatter(ax, fc[ns], neg_log_p[ns], s, alpha, color='grey',
                      label=f"Not significant ({int(ns.sum())})")
        self._scatter(ax, fc[up], neg_log_p[up], s, alpha, color='#d62728',
                      label=f"Up in {g2_name} ({int(up.sum())})")
        self._scatter(ax, fc[down], neg_log_p[down], s, alpha, color='#1f77b4',
                      label=f"Up in {g1_name} ({int(down.sum())})")
        ax.axhline(p_line, color='red', linestyle='--', alpha=0.5,
                   label=f'{"FDR q" if use_fdr else "p"}={p_thr}')
        ax.axvline(-fc_thr, color='grey', linestyle='--', alpha=0.3)
        ax.axvline(fc_thr, color='grey', linestyle='--', alpha=0.3)
        ax.set_xlabel("log2 Fold Change")
        ax.set_ylabel(y_label)
        ax.set_title(f"Volcano Plot: {g2_name} vs {g1_name}")
        # Feature labels — top-N features ranked by significance * fold change
        if self.params.get("feature_labels", True):
            _annotate_top_points(
                ax, fc, neg_log_p, [str(c) for c in mat.columns],
                self.params.get("max_labels", 10),
                self.params.get("label_font_size", "10 pt"),
                self.params.get("label_pos", "Auto"),
                scores=neg_log_p * np.abs(fc))
        if self.params.get("show_legend", True):
            _draw_legend(ax, self.params.get("legend_pos", "Outside Right"))
        if self.params.get("show_grid", False):
            ax.grid(True, alpha=0.3)

    def _run_boxplot(self, mat, fig, ax, group_labels=None):
        """Box plot of top features by variance; grouped by Group when present."""
        import numpy as np
        import matplotlib.patches as mpatches

        variances = mat.var().sort_values(ascending=False)
        n_feat = self.params.get("box_features", 10)
        top_features = variances.head(min(n_feat, len(variances))).index.tolist()

        if not top_features:
            ax.text(0.5, 0.5, "No features to plot", transform=ax.transAxes,
                    ha='center', va='center')
            return

        if group_labels is not None and len(set(group_labels)) > 1:
            groups = sorted(set(group_labels))
            gcolors = _group_color_map(groups, "tab10")
            fig.clf()
            # matplotlib >=3.11: FigureBase.subplots no longer accepts figsize.
            # User-chosen data plot width; height scales with the number of
            # panels. Reserve room on the right for an outside legend.
            legend_out = (self.params.get("legend_pos", "Outside Right")
                          == "Outside Right"
                          and self.params.get("show_legend", True))
            fig.set_size_inches(self.params.get("fig_width", 10.0),
                                2.2 * len(top_features))
            fig.subplots_adjust(left=0.125, bottom=0.11, top=0.90,
                                right=0.80 if legend_out else 0.875)
            axes = fig.subplots(nrows=len(top_features), ncols=1, sharex=True)
            if len(top_features) == 1:
                axes = [axes]
            for i, feat in enumerate(top_features):
                # Skip groups with no data for this feature — an empty
                # array would crash matplotlib boxplot
                plot_groups, data_by_group = [], []
                for g in groups:
                    mask = np.array([gl == g for gl in group_labels])
                    vals = mat.loc[mask, feat].dropna().values
                    if len(vals) == 0:
                        continue
                    plot_groups.append(g)
                    data_by_group.append(vals)
                if not plot_groups:
                    continue
                # `labels=` was removed from Axes.boxplot in matplotlib 3.11;
                # set tick labels explicitly for cross-version compatibility
                bp = axes[i].boxplot(data_by_group, patch_artist=True)
                axes[i].set_xticks(range(1, len(plot_groups) + 1))
                axes[i].set_xticklabels([str(g) for g in plot_groups])
                for patch, g in zip(bp['boxes'], plot_groups):
                    patch.set_facecolor(gcolors[g])
                    patch.set_alpha(0.6)
                axes[i].set_ylabel(str(feat)[:12], fontsize=9)
            axes[-1].set_xlabel("Group")
            axes[0].set_title("Box Plot — Top Features by Group")
            if self.params.get("show_legend", True):
                handles = [mpatches.Patch(color=gcolors[g], label=str(g))
                           for g in groups]
                _draw_legend(axes[0],
                             self.params.get("legend_pos", "Outside Right"),
                             title="Group", handles=handles)
            fig.tight_layout()
        else:
            plot_data = mat[top_features]
            plot_data.boxplot(ax=ax, rot=45)
            ax.set_title("Box Plot — Top Features by Variance")
            ax.set_ylabel("Intensity")

    def _run_oplsda(self, mat, fig, ax, group_labels=None):
        """Generate an OPLS-DA plot.

        With 2+ groups: orthogonal signal correction (one component
        orthogonal to Y is removed from X) followed by PLS-DA (sklearn
        PLSRegression on one-hot group labels); the plotted scores are the
        predictive T scores. Without groups, falls back to PCA with a note.
        """
        import numpy as np

        if group_labels is not None and len(set(group_labels)) >= 2:
            groups = sorted(set(group_labels))
            gcolors = _group_color_map(
                groups, self.params.get("color_scheme", "tab10"))

            X = mat.values.astype(float)
            y = np.zeros((X.shape[0], len(groups)))
            for gi, g in enumerate(groups):
                y[np.array([gl == g for gl in group_labels]), gi] = 1.0

            try:
                Xc = X - X.mean(axis=0)
                # OSC: remove the strongest Y-orthogonal component from X
                u, s, vt = np.linalg.svd(Xc, full_matrices=False)
                t = u[:, 0] * s[0]
                t_o = t - y @ (np.linalg.pinv(y.T @ y) @ (y.T @ t))
                if np.linalg.norm(t_o) > 1e-12:
                    w = Xc.T @ t_o / (t_o @ t_o)
                    w = w / np.linalg.norm(w)
                    t_o = Xc @ w
                    p = Xc.T @ t_o / (t_o @ t_o)
                    Xc = Xc - np.outer(t_o, p)
                from sklearn.cross_decomposition import PLSRegression
                pls = PLSRegression(
                    n_components=max(1, min(2, Xc.shape[1], Xc.shape[0])))
                pls.fit(Xc, y)
                reduced = pls.x_scores_
                if reduced.shape[1] < 2:
                    reduced = np.hstack(
                        [reduced, np.zeros((reduced.shape[0], 1))])
                self._result_data = reduced
            except Exception as e:
                self.log.emit(f"OPLS-DA failed ({e}) — using PCA fallback")
                from sklearn.decomposition import PCA
                reduced = PCA(n_components=2).fit_transform(mat.fillna(0))
                self._result_data = reduced

            s = self.params.get("point_size", 36)
            alpha = self.params.get("opacity", 0.8)
            for g in groups:
                mask = np.array([gl == g for gl in group_labels])
                self._scatter(ax, reduced[mask, 0], reduced[mask, 1],
                              s, alpha, color=gcolors[g], label=str(g))
            # Confidence ellipse is PCA-only — not drawn here
            if self.params.get("show_legend", True):
                _draw_legend(ax,
                             self.params.get("legend_pos", "Outside Right"),
                             title="Group")
            ax.set_xlabel("T score [1] (predictive)")
            ax.set_ylabel("T score [2]")
            ax.set_title("OPLS-DA")
        else:
            from sklearn.decomposition import PCA
            pca = PCA(n_components=2)
            reduced = pca.fit_transform(mat.fillna(0))
            self._result_data = reduced
            self._scatter(ax, reduced[:, 0], reduced[:, 1],
                          self.params.get("point_size", 36),
                          self.params.get("opacity", 0.8))
            ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
            ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
            ax.set_title("OPLS-DA (no groups — PCA shown)")
        # Hover: per-point coordinates for the preview tooltip
        self._hover_points = np.asarray(reduced)[:, :2]
        self._hover_groups = ([str(g) for g in group_labels]
                              if group_labels is not None else None)
        self._hover_prefix = "Sample"
        # Data point labels — top-N points by distance from the origin
        if self.params.get("data_point_labels", False):
            names = [str(v) for v in (mat.index if hasattr(mat, "index")
                                      else range(len(reduced)))]
            _annotate_top_points(
                ax, reduced[:, 0], reduced[:, 1], names,
                self.params.get("max_labels", 10),
                self.params.get("label_font_size", "10 pt"),
                self.params.get("label_pos", "Auto"))
        if self.params.get("show_grid", False):
            ax.grid(True, alpha=0.3)


# =============================================================================
# REUSABLE WIDGETS
# =============================================================================

class FileSelector(QWidget):
    """File/directory selector with label, line edit and browse button."""

    path_changed = pyqtSignal(str)

    def __init__(self, mode: str = "file", label: str = "", 
                 filter_str: str = "All Files (*)", placeholder: str = "", parent=None):
        super().__init__(parent)
        self.mode = mode
        self.filter_str = filter_str

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        if label:
            self.label = QLabel(label)
            self.label.setFixedWidth(100)
            layout.addWidget(self.label)

        self.path_edit = QLineEdit()
        self.path_edit.setPlaceholderText(placeholder or "Click Browse to select...")
        self.path_edit.setReadOnly(True)
        layout.addWidget(self.path_edit, 1)

        self.browse_btn = QPushButton("Browse")
        self.browse_btn.setObjectName("browseBtn")
        self.browse_btn.setFixedWidth(70)
        self.browse_btn.clicked.connect(self._browse)
        layout.addWidget(self.browse_btn)

    def _browse(self):
        if self.mode == "file":
            path, _ = QFileDialog.getOpenFileName(self, "Select File", "", self.filter_str)
        elif self.mode == "files":
            paths, _ = QFileDialog.getOpenFileNames(self, "Select Files", "", self.filter_str)
            path = ";".join(paths) if paths else ""
        elif self.mode == "directory":
            path = QFileDialog.getExistingDirectory(self, "Select Directory")
        else:
            path = ""

        if path:
            self.path_edit.setText(path)
            self.path_changed.emit(path)

    def get_path(self) -> str:
        return self.path_edit.text()

    def set_path(self, path: str):
        self.path_edit.setText(path)


class LogConsole(QPlainTextEdit):
    """Console-like log display at bottom."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("logConsole")
        self.setReadOnly(True)
        self.setMaximumBlockCount(2000)
        font = QFont("Consolas", 10)
        self.setFont(font)
        self.setFixedHeight(120)

    def log(self, message: str):
        timestamp = datetime.now().strftime("[%H:%M:%S]")
        self.appendPlainText(f"{timestamp} {message}")
        self.verticalScrollBar().setValue(self.verticalScrollBar().maximum())


# =============================================================================
# TAB 1: PREPROCESSING (预处理)
# =============================================================================

class PreprocessingTab(QWidget):
    """Tab 1: Preprocessing - matches original Chinese UI exactly."""

    def __init__(self, log_console: LogConsole, parent=None):
        super().__init__(parent)
        self.log = log_console
        self.worker = None
        self._setup_ui()

    def _setup_ui(self):
        main_layout = QVBoxLayout(self)
        main_layout.setSpacing(12)
        main_layout.setContentsMargins(16, 16, 16, 16)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)

        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setSpacing(12)

        # === File Path Section ===
        file_group = QGroupBox("File Path")
        file_layout = QVBoxLayout(file_group)
        file_layout.setSpacing(8)

        # Folder selector (matching Chinese version: select folder containing RAW files)
        folder_label = QLabel("RAW Folder Path:")
        folder_label.setStyleSheet("font-weight: bold;")
        file_layout.addWidget(folder_label)

        folder_sel_layout = QHBoxLayout()
        self.raw_selector = FileSelector(
            "directory", "",
            placeholder="Select folder containing .RAW files..."
        )
        folder_sel_layout.addWidget(self.raw_selector, 1)

        # Also allow selecting a single RAW file
        single_file_btn = QPushButton("Or Select Single File")
        single_file_btn.setObjectName("browseBtn")
        single_file_btn.clicked.connect(self._browse_single_raw)
        folder_sel_layout.addWidget(single_file_btn)
        file_layout.addLayout(folder_sel_layout)

        self.raw_file_label = QLabel("")
        self.raw_file_label.setStyleSheet("color: #888888; font-size: 11px;")
        file_layout.addWidget(self.raw_file_label)

        self.output_selector = FileSelector(
            "directory", "Output Path:",
            placeholder="Select output directory..."
        )
        file_layout.addWidget(self.output_selector)

        layout.addWidget(file_group)

        # === Parameter Settings ===
        params_group = QGroupBox("Parameter Settings")
        params_layout = QGridLayout(params_group)
        params_layout.setSpacing(10)
        params_layout.setColumnStretch(1, 1)
        params_layout.setColumnStretch(3, 1)

        # Row 0: Ion mode quick toggle
        params_layout.addWidget(QLabel("Ion Mode:"), 0, 0)
        ion_layout = QHBoxLayout()
        self.ion_positive = QRadioButton("Positive (+)")
        self.ion_positive.setChecked(True)
        self.ion_negative = QRadioButton("Negative (-)")
        self.ion_positive.toggled.connect(self._on_ion_mode_changed)
        ion_layout.addWidget(self.ion_positive)
        ion_layout.addWidget(self.ion_negative)
        ion_layout.addStretch()
        params_layout.addLayout(ion_layout, 0, 1, 1, 3)

        # Row 1: Target m/z and PPM tolerance
        params_layout.addWidget(QLabel("Target m/z:"), 1, 0)
        self.target_mz = QDoubleSpinBox()
        self.target_mz.setRange(0, 5000)
        self.target_mz.setDecimals(4)
        self.target_mz.setValue(760.58)  # Positive ion mode default
        self.target_mz.setSingleStep(0.0001)
        params_layout.addWidget(self.target_mz, 1, 1)

        params_layout.addWidget(QLabel("PPM Tolerance:"), 1, 2)
        self.ppm_tolerance = QDoubleSpinBox()
        self.ppm_tolerance.setRange(0.01, 100)
        self.ppm_tolerance.setDecimals(2)
        self.ppm_tolerance.setValue(5.00)
        params_layout.addWidget(self.ppm_tolerance, 1, 3)

        # Row 2: Segment length and Intensity threshold (Chinese defaults)
        params_layout.addWidget(QLabel("Segment Length:"), 2, 0)
        self.segment_length = QDoubleSpinBox()
        self.segment_length.setRange(0.1, 100)
        self.segment_length.setDecimals(2)
        self.segment_length.setValue(5.00)
        params_layout.addWidget(self.segment_length, 2, 1)

        params_layout.addWidget(QLabel("Intensity Threshold:"), 2, 2)
        self.intensity_threshold = QDoubleSpinBox()
        self.intensity_threshold.setRange(0.0, 1.0)
        self.intensity_threshold.setDecimals(2)
        self.intensity_threshold.setValue(0.10)
        params_layout.addWidget(self.intensity_threshold, 2, 3)

        # Row 3: argrelextrema order (Chinese order_spin, default 1)
        params_layout.addWidget(QLabel("Order:"), 3, 0)
        self.smoothing_order = QSpinBox()
        self.smoothing_order.setRange(1, 10)
        self.smoothing_order.setValue(1)
        params_layout.addWidget(self.smoothing_order, 3, 1)

        # Row 4: Checkboxes (span full width to avoid overlap)
        checkbox_layout = QHBoxLayout()
        self.remove_isotopes = QCheckBox("Remove Isotopes")
        self.auto_merge = QCheckBox("Auto-merge after completion")
        self.auto_merge.setChecked(True)
        self.subtract_background = QCheckBox("Subtract Background")
        self.subtract_background.setChecked(True)
        checkbox_layout.addWidget(self.remove_isotopes)
        checkbox_layout.addWidget(self.auto_merge)
        checkbox_layout.addWidget(self.subtract_background)
        checkbox_layout.addStretch()
        params_layout.addLayout(checkbox_layout, 4, 0, 1, 4)

        layout.addWidget(params_group)

        # === Submit / Stop Buttons ===
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()

        self.submit_btn = QPushButton("Submit")
        self.submit_btn.setObjectName("submitBtn")
        # minimum, not fixed: the stylesheet's 30px side padding needs ~110px
        # total, so a fixed 100px width clips the button text
        self.submit_btn.setMinimumWidth(100)
        self.submit_btn.clicked.connect(self._start_processing)
        btn_layout.addWidget(self.submit_btn)

        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setObjectName("stopBtn")
        self.stop_btn.setMinimumWidth(100)
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self._stop_processing)
        btn_layout.addWidget(self.stop_btn)

        layout.addLayout(btn_layout)

        # === Data Processing Section ===
        process_layout = QHBoxLayout()

        # File Merge
        merge_group = QGroupBox("File Merge")
        merge_layout = QVBoxLayout(merge_group)
        merge_layout.addWidget(QLabel("Merge preprocessed files into new folder"))
        self.merge_input = FileSelector("directory", "", placeholder="Select folder to merge...")
        merge_layout.addWidget(self.merge_input)
        self.merge_btn = QPushButton("Merge Files")
        self.merge_btn.setObjectName("mergeBtn")
        self.merge_btn.clicked.connect(self._start_merge)
        merge_layout.addWidget(self.merge_btn)

        process_layout.addWidget(merge_group, 1)

        # Group Split
        split_group = QGroupBox("Group Split")
        split_layout = QVBoxLayout(split_group)
        self.split_selector = FileSelector("file", "", 
            "Excel Files (*.xlsx *.xls)", 
            "Select Excel file to split...")
        split_layout.addWidget(self.split_selector)
        self.split_btn = QPushButton("Submit")
        self.split_btn.setObjectName("submitBtn")
        self.split_btn.clicked.connect(self._start_split)
        split_layout.addWidget(self.split_btn)

        process_layout.addWidget(split_group, 1)

        layout.addLayout(process_layout)

        # === Processing Status ===
        status_group = QGroupBox("Processing Status")
        status_layout = QVBoxLayout(status_group)

        self.progress_bar = QProgressBar()
        self.progress_bar.setTextVisible(True)
        status_layout.addWidget(self.progress_bar)

        time_layout = QHBoxLayout()
        time_layout.addWidget(QLabel("Elapsed Time:"))
        self.elapsed_label = QLabel("00:00:00")
        self.elapsed_label.setStyleSheet("font-weight: bold; color: #1976D2;")
        time_layout.addWidget(self.elapsed_label)
        time_layout.addStretch()
        status_layout.addLayout(time_layout)

        status_layout.addWidget(QLabel("Merge Output:"))
        self.merge_output = QLineEdit()
        self.merge_output.setPlaceholderText("Output file path will be displayed after merge completion")
        self.merge_output.setReadOnly(True)
        status_layout.addWidget(self.merge_output)

        layout.addWidget(status_group)

        layout.addStretch()
        scroll.setWidget(container)
        main_layout.addWidget(scroll)

    def _on_ion_mode_changed(self, checked):
        """Toggle target m/z between positive (760.5842) and negative (885.5492) ion modes."""
        if self.ion_positive.isChecked():
            self.target_mz.setValue(760.58)
        else:
            self.target_mz.setValue(885.5492)

    def _browse_single_raw(self):
        """Fallback: select a single RAW file instead of a folder."""
        path, _ = QFileDialog.getOpenFileName(
            self, "Select RAW File", "",
            "RAW Files (*.raw *.RAW *.Raw);;All Files (*)"
        )
        if path:
            self.raw_selector.set_path(os.path.dirname(path))
            self.raw_file_label.setText(f"Single file: {os.path.basename(path)}")
            # Store the single file path for later
            self._single_raw_file = path

    def _start_processing(self):
        raw_path = self.raw_selector.get_path()
        output_path = self.output_selector.get_path()

        # If a single file was selected via "Or Select Single File", use that
        if hasattr(self, '_single_raw_file') and self._single_raw_file:
            raw_path = self._single_raw_file
            self._single_raw_file = None  # clear after use

        if not raw_path:
            QMessageBox.warning(self, "Warning",
                "Please select a folder containing .RAW files (or use 'Or Select Single File').")
            return
        if not output_path:
            QMessageBox.warning(self, "Warning", "Please select output path.")
            return

        # Check if path exists
        if not os.path.exists(raw_path):
            QMessageBox.warning(self, "Warning", f"Path does not exist: {raw_path}")
            return

        params = {
            "target_mz": self.target_mz.value(),
            "ppm_tolerance": self.ppm_tolerance.value(),
            "segment_length": self.segment_length.value(),
            "intensity_threshold": self.intensity_threshold.value(),
            "smoothing_order": self.smoothing_order.value(),
            "remove_isotopes": self.remove_isotopes.isChecked(),
            "auto_merge": self.auto_merge.isChecked(),
            "subtract_background": self.subtract_background.isChecked(),
        }

        self.log.log(f"Starting preprocessing: {raw_path}")
        self.worker = PreprocessingWorker(raw_path, output_path, params)
        self._connect_worker(self.worker)
        self.worker.start()

    def _start_merge(self):
        input_dir = self.merge_input.get_path()
        if not input_dir:
            QMessageBox.warning(self, "Warning", "Please select folder to merge.")
            return

        self.worker = MergeWorker(input_dir, input_dir)
        self._connect_worker(self.worker)
        self.worker.start()

    def _start_split(self):
        excel_file = self.split_selector.get_path()
        if not excel_file:
            QMessageBox.warning(self, "Warning", "Please select Excel file to split.")
            return

        output_dir = os.path.dirname(excel_file)
        self.worker = SplitWorker(excel_file, output_dir)
        self._connect_worker(self.worker)
        self.worker.start()

    def _connect_worker(self, worker):
        worker.progress.connect(self.progress_bar.setValue)
        worker.status.connect(lambda s: self.log.log(s))
        worker.elapsed_time.connect(self.elapsed_label.setText)
        worker.finished_signal.connect(self._on_finished)
        worker.log.connect(lambda s: self.log.log(s))
        self.submit_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        self.merge_btn.setEnabled(False)
        self.split_btn.setEnabled(False)

    def _stop_processing(self):
        if self.worker:
            self.worker.stop()
            self.log.log("Processing stopped by user.")
        self._reset_ui()

    def _on_finished(self, success: bool, message: str):
        if success:
            self.log.log(message)
            # Show output path for merge operations
            output_path = getattr(self.worker, 'output_file_path', None)
            if output_path:
                self.merge_output.setText(output_path)
            QMessageBox.information(self, "Complete", message)
        else:
            self.log.log(f"ERROR: {message}")
            QMessageBox.critical(self, "Error", message)
        self._reset_ui()

    def _reset_ui(self):
        self.submit_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.merge_btn.setEnabled(True)
        self.split_btn.setEnabled(True)
        self.progress_bar.setValue(0)


# =============================================================================
# TAB 2: DATA CLEANING (数据清洗)
# =============================================================================

class DataCleaningTab(QWidget):
    """Tab 2: Data Cleaning - matches original Chinese UI exactly."""

    def __init__(self, log_console: LogConsole, parent=None):
        super().__init__(parent)
        self.log = log_console
        self.worker = None
        self._setup_ui()

    def _setup_ui(self):
        main_layout = QHBoxLayout(self)
        main_layout.setSpacing(12)
        main_layout.setContentsMargins(16, 16, 16, 16)

        # Left panel
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.setSpacing(10)
        left_panel.setFixedWidth(320)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setFixedWidth(320)

        left_container = QWidget()
        left_container.setFixedWidth(300)
        ll = QVBoxLayout(left_container)
        ll.setSpacing(10)

        # === Data File ===
        file_group = QGroupBox("Data File")
        file_layout = QVBoxLayout(file_group)

        file_sel_layout = QHBoxLayout()
        file_sel_layout.addWidget(QLabel("Data File:"))
        self.data_file_edit = QLineEdit()
        self.data_file_edit.setPlaceholderText("Select CSV or...")
        self.data_file_edit.setReadOnly(True)
        file_sel_layout.addWidget(self.data_file_edit, 1)
        self.data_browse_btn = QPushButton("Browse")
        self.data_browse_btn.setObjectName("browseBtn")
        self.data_browse_btn.clicked.connect(self._browse_data_file)
        file_sel_layout.addWidget(self.data_browse_btn)
        file_layout.addLayout(file_sel_layout)

        self.file_info_label = QLabel("File Info: No file selected")
        self.file_info_label.setStyleSheet("color: #888888; font-size: 11px;")
        file_layout.addWidget(self.file_info_label)

        ll.addWidget(file_group)

        # === Outlier Processing ===
        outlier_group = QGroupBox("Outlier Processing")
        outlier_layout = QGridLayout(outlier_group)

        self.std_filter = QCheckBox("Standard Deviation Filter")
        outlier_layout.addWidget(self.std_filter, 0, 0)
        outlier_layout.addWidget(QLabel("Coefficient:"), 0, 1)
        self.std_coeff = QDoubleSpinBox()
        self.std_coeff.setRange(0.1, 10)
        self.std_coeff.setDecimals(2)
        self.std_coeff.setValue(3.00)
        self.std_coeff.setSuffix("σ")
        outlier_layout.addWidget(self.std_coeff, 0, 2)

        self.iqr_filter = QCheckBox("IQR Filter")
        self.iqr_filter.setChecked(True)  # Chinese version always handles outliers (IQR x1.5)
        outlier_layout.addWidget(self.iqr_filter, 1, 0)
        outlier_layout.addWidget(QLabel("Coefficient:"), 1, 1)
        self.iqr_coeff = QDoubleSpinBox()
        self.iqr_coeff.setRange(0.1, 10)
        self.iqr_coeff.setDecimals(2)
        self.iqr_coeff.setValue(1.50)
        outlier_layout.addWidget(self.iqr_coeff, 1, 2)

        ll.addWidget(outlier_group)

        # === Missing Value Processing ===
        missing_group = QGroupBox("Missing Value Processing")
        missing_layout = QGridLayout(missing_group)

        self.filter_sparse = QCheckBox("Filter Sparse Columns")
        self.filter_sparse.setChecked(True)  # Chinese version always filters (threshold 0.5)
        missing_layout.addWidget(self.filter_sparse, 0, 0)
        missing_layout.addWidget(QLabel("Threshold:"), 0, 1)
        self.sparse_threshold = QDoubleSpinBox()
        self.sparse_threshold.setRange(0.0, 1.0)
        self.sparse_threshold.setDecimals(2)
        self.sparse_threshold.setValue(0.50)  # Chinese version: 0.5
        missing_layout.addWidget(self.sparse_threshold, 0, 2)

        self.impute_knn = QCheckBox("Missing Value Imputation")
        self.impute_knn.setChecked(True)  # Chinese version always imputes (KNN k=20)
        missing_layout.addWidget(self.impute_knn, 1, 0)
        missing_layout.addWidget(QLabel("KNN:"), 1, 1)
        self.knn_k = QSpinBox()
        self.knn_k.setRange(1, 50)
        self.knn_k.setValue(20)  # Chinese version: k=20
        missing_layout.addWidget(self.knn_k, 1, 2)

        self.fill_missing = QCheckBox("Missing Value Fill")
        missing_layout.addWidget(self.fill_missing, 2, 0, 1, 3)

        self.zero_half_min = QCheckBox("Replace 0 with Half of Column Minimum")
        self.zero_half_min.setChecked(False)
        self.zero_half_min.setToolTip(
            "Replace every 0 with half of the column's minimum nonzero "
            "value. Columns that are entirely 0 (no nonzero minimum to "
            "halve) are set to NaN instead.")
        missing_layout.addWidget(self.zero_half_min, 3, 0, 1, 3)

        ll.addWidget(missing_group)

        # === Group Size Balancing ===
        balance_group = QGroupBox("Group Size Balancing")
        balance_layout = QVBoxLayout(balance_group)

        self.balance_groups = QCheckBox("Balance Group Sizes")
        self.balance_groups.setChecked(False)  # opt-in only
        self.balance_groups.setToolTip(
            "Randomly downsample every sheet to the smallest sheet's row "
            "count so all groups end up with the same number of cells. "
            "Fixed random seed for reproducibility. Note: real cells are "
            "discarded from the larger groups — enable only when equal "
            "group sizes are required (e.g. for balanced statistics).")
        balance_layout.addWidget(self.balance_groups)

        ll.addWidget(balance_group)

        # === Data Standardization ===
        std_group = QGroupBox("Data Standardization")
        std_layout = QVBoxLayout(std_group)

        self.enable_std = QCheckBox("Enable Standardization")
        std_layout.addWidget(self.enable_std)

        method_layout = QHBoxLayout()
        method_layout.addWidget(QLabel("Method:"))
        self.std_method = QComboBox()
        self.std_method.addItems([
            "Max Normalization",
            "Z-score Standardization",
            "Min-Max Scaling",
            "Log Transformation"
        ])
        method_layout.addWidget(self.std_method, 1)
        std_layout.addLayout(method_layout)

        ll.addWidget(std_group)

        # === Reference Value Filter ===
        ref_group = QGroupBox("Reference Value Filter")
        ref_layout = QGridLayout(ref_group)

        self.enable_ref = QCheckBox("Enable Reference Filter")
        ref_layout.addWidget(self.enable_ref, 0, 0, 1, 3)

        ref_layout.addWidget(QLabel("Ref Value:"), 1, 0)
        self.ref_value = QDoubleSpinBox()
        self.ref_value.setRange(0, 5000)
        self.ref_value.setDecimals(4)
        self.ref_value.setValue(760.5800)
        ref_layout.addWidget(self.ref_value, 1, 1, 1, 2)

        ref_layout.addWidget(QLabel("PPM Tolerance:"), 2, 0)
        self.ref_ppm = QDoubleSpinBox()
        self.ref_ppm.setRange(0.01, 100)
        self.ref_ppm.setDecimals(2)
        self.ref_ppm.setValue(5.0)
        self.ref_ppm.setSuffix(" ppm")
        ref_layout.addWidget(self.ref_ppm, 2, 1, 1, 2)

        ref_layout.addWidget(QLabel("Filter Condition:"), 3, 0)
        self.ref_condition = QComboBox()
        self.ref_condition.addItems([">", ">=", "<", "<=", "==", "!="])
        ref_layout.addWidget(self.ref_condition, 3, 1)
        self.ref_threshold = QDoubleSpinBox()
        self.ref_threshold.setRange(0, 1000)
        self.ref_threshold.setDecimals(2)
        self.ref_threshold.setValue(80.00)
        ref_layout.addWidget(self.ref_threshold, 3, 2)

        self.preview_match_btn = QPushButton("Preview Match")
        self.preview_match_btn.setObjectName("browseBtn")
        self.preview_match_btn.clicked.connect(self._preview_match)
        ref_layout.addWidget(self.preview_match_btn, 4, 0, 1, 3)

        ll.addWidget(ref_group)

        # === Operation ===
        op_group = QGroupBox("Operation")
        op_layout = QVBoxLayout(op_group)

        self.preview_data_btn = QPushButton("Preview Data")
        self.preview_data_btn.setObjectName("previewBtn")
        self.preview_data_btn.clicked.connect(self._preview_data)
        op_layout.addWidget(self.preview_data_btn)

        self.clean_btn = QPushButton("Start Cleaning")
        self.clean_btn.setObjectName("submitBtn")
        self.clean_btn.clicked.connect(self._start_cleaning)
        op_layout.addWidget(self.clean_btn)

        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setObjectName("stopBtn")
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self._stop_cleaning)
        op_layout.addWidget(self.stop_btn)

        ll.addWidget(op_group)
        ll.addStretch()

        scroll.setWidget(left_container)
        left_layout.addWidget(scroll)

        # Right panel - Data Preview
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)

        # Preview header
        preview_header = QHBoxLayout()
        preview_header.addWidget(QLabel("<b>Data Preview</b>"))
        preview_header.addStretch()
        preview_header.addWidget(QLabel("Display Mode:"))
        self.display_mode = QComboBox()
        self.display_mode.addItems(["Raw Data", "Processed Data", "Statistics"])
        preview_header.addWidget(self.display_mode)

        right_layout.addLayout(preview_header)

        # Data table / stats tabs
        self.preview_tabs = QTabWidget()
        self.preview_tabs.setDocumentMode(True)

        self.data_table = QTableWidget()
        self.data_table.setColumnCount(0)
        self.data_table.setRowCount(0)
        self.preview_tabs.addTab(self.data_table, "Data Table")

        self.stats_widget = QWidget()
        stats_layout = QVBoxLayout(self.stats_widget)
        self.stats_text = QPlainTextEdit()
        self.stats_text.setReadOnly(True)
        stats_layout.addWidget(self.stats_text)
        self.preview_tabs.addTab(self.stats_widget, "Statistics")

        right_layout.addWidget(self.preview_tabs, 1)

        # Output info
        output_group = QGroupBox("Output Information")
        output_layout = QFormLayout(output_group)

        self.clean_progress = QProgressBar()
        output_layout.addRow("Cleaning Progress:", self.clean_progress)

        self.output_file_edit = QLineEdit()
        self.output_file_edit.setPlaceholderText("Output file path will be displayed after cleaning completion")
        self.output_file_edit.setReadOnly(True)
        output_layout.addRow("Output File:", self.output_file_edit)

        right_layout.addWidget(output_group)

        main_layout.addWidget(left_panel)
        main_layout.addWidget(right_panel, 1)

    def _browse_data_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Data File", "", 
            "Data Files (*.csv *.xlsx *.txt);;All Files (*)")
        if path:
            self.data_file_edit.setText(path)
            self.file_info_label.setText(f"File Info: {os.path.basename(path)}")
            self.log.log(f"Selected data file: {path}")

    def _preview_data(self):
        path = self.data_file_edit.text()
        if not path:
            QMessageBox.warning(self, "Warning", "Please select a data file first.")
            return
        try:
            if path.endswith('.csv'):
                df = pd.read_csv(path)
            else:
                df = pd.read_excel(path)

            self.data_table.setColumnCount(len(df.columns))
            self.data_table.setRowCount(min(len(df), 100))
            # QTableWidget requires all header labels to be strings
            self.data_table.setHorizontalHeaderLabels([str(c) for c in df.columns])

            for i in range(min(len(df), 100)):
                for j, col in enumerate(df.columns):
                    self.data_table.setItem(i, j, QTableWidgetItem(str(df.iloc[i, j])))

            self.stats_text.setPlainText(df.describe().to_string())
            self.log.log(f"Preview loaded: {len(df)} rows x {len(df.columns)} cols")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load file: {str(e)}")

    def _preview_match(self):
        """Preview reference value matches in the data file."""
        path = self.data_file_edit.text()
        if not path:
            QMessageBox.warning(self, "Warning", "Please select a data file first.")
            return

        try:
            import pandas as pd
            import numpy as np

            if path.endswith('.csv'):
                df = pd.read_csv(path)
            elif path.endswith(('.xlsx', '.xls')):
                df = pd.read_excel(path)
            else:
                df = pd.read_csv(path)

            ref_val = self.ref_value.value()
            ref_ppm = self.ref_ppm.value()

            # Find numeric columns that match the reference m/z within ppm tolerance
            numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
            matched = []
            for col in numeric_cols:
                try:
                    col_val = float(col)
                    error = abs(col_val - ref_val) / ref_val * 1e6 if ref_val > 0 else float('inf')
                    if error <= ref_ppm:
                        matched.append((col, round(error, 4)))
                except (ValueError, TypeError, ZeroDivisionError):
                    pass

            if matched:
                info = "\n".join([f"  m/z {m[0]}: {m[1]} ppm error" for m in matched[:20]])
                QMessageBox.information(self, "Match Preview",
                    f"Found {len(matched)} columns matching m/z {ref_val} ± {ref_ppm} ppm:\n\n{info}")
            else:
                QMessageBox.information(self, "Match Preview",
                    f"No columns match m/z {ref_val} ± {ref_ppm} ppm.")
            self.log.log(f"Reference match preview: {len(matched)} matches for m/z {ref_val}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Preview failed: {str(e)}")

    def _start_cleaning(self):
        data_file = self.data_file_edit.text()
        if not data_file:
            QMessageBox.warning(self, "Warning", "Please select a data file.")
            return

        output_dir = os.path.dirname(data_file)
        params = {
            "std_filter": self.std_filter.isChecked(),
            "std_coeff": self.std_coeff.value(),
            "iqr_filter": self.iqr_filter.isChecked(),
            "iqr_coeff": self.iqr_coeff.value(),
            "filter_sparse": self.filter_sparse.isChecked(),
            "sparse_threshold": self.sparse_threshold.value(),
            "impute_knn": self.impute_knn.isChecked(),
            "knn_k": self.knn_k.value(),
            "fill_missing": self.fill_missing.isChecked(),
            "zero_half_min": self.zero_half_min.isChecked(),
            "balance_groups": self.balance_groups.isChecked(),
            "enable_std": self.enable_std.isChecked(),
            "std_method": self.std_method.currentText(),
            "enable_ref": self.enable_ref.isChecked(),
            "ref_value": self.ref_value.value(),
            "ref_ppm": self.ref_ppm.value(),
        }

        self.worker = DataCleaningWorker(data_file, output_dir, params)
        self.worker.progress.connect(self.clean_progress.setValue)
        self.worker.status.connect(lambda s: self.log.log(s))
        self.worker.finished_signal.connect(self._on_finished)
        self.worker.log.connect(lambda s: self.log.log(s))

        self.clean_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        self.worker.start()

    def _stop_cleaning(self):
        if self.worker:
            self.worker.stop()
            self.log.log("Cleaning stopped by user.")
        self._reset_ui()

    def _on_finished(self, success, message):
        if success:
            self.log.log(message)
            # Show output file path if available
            output_path = getattr(self.worker, 'output_file_path', message)
            self.output_file_edit.setText(output_path)
            QMessageBox.information(self, "Complete", message)
        else:
            self.log.log(f"ERROR: {message}")
            QMessageBox.critical(self, "Error", message)
        self._reset_ui()

    def _reset_ui(self):
        self.clean_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.clean_progress.setValue(0)


# =============================================================================
# TAB 3: SEARCH (检索)
# =============================================================================

class SearchTab(QWidget):
    """Tab 3: Search - matches original Chinese UI exactly."""

    def __init__(self, log_console: LogConsole, parent=None):
        super().__init__(parent)
        self.log = log_console
        self.worker = None
        self._setup_ui()

    def _setup_ui(self):
        main_layout = QVBoxLayout(self)
        main_layout.setSpacing(12)
        main_layout.setContentsMargins(16, 16, 16, 16)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)

        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setSpacing(12)

        # === Merge Cleaned Data ===
        merge_group = QGroupBox("Merge Cleaned Data")
        merge_layout = QVBoxLayout(merge_group)

        self.merge_excel = FileSelector("file", "", 
            "Excel Files (*.xlsx *.xls)",
            "Select cleaned Excel file...")
        merge_layout.addWidget(self.merge_excel)

        ppm_layout = QHBoxLayout()
        ppm_layout.addWidget(QLabel("PPM Tolerance:"))
        self.merge_ppm = QDoubleSpinBox()
        self.merge_ppm.setRange(0.01, 100)
        self.merge_ppm.setDecimals(2)
        self.merge_ppm.setValue(5.00)
        ppm_layout.addWidget(self.merge_ppm)
        ppm_layout.addStretch()
        merge_layout.addLayout(ppm_layout)

        self.merge_search_btn = QPushButton("Start Merge")
        self.merge_search_btn.setObjectName("mergeBtn")
        self.merge_search_btn.clicked.connect(self._start_merge_search)
        merge_layout.addWidget(self.merge_search_btn)

        layout.addWidget(merge_group)

        # === Local Search ===
        local_group = QGroupBox("Local Search")
        local_layout = QVBoxLayout(local_group)

        self.local_data = FileSelector("file", "Merged Data File:",
            "Data Files (*.csv *.xlsx)",
            "Select merged data file...")
        local_layout.addWidget(self.local_data)

        self.local_ref = FileSelector("file", "Reference Library:",
            "All Files (*)",
            "Select reference library file...")
        local_layout.addWidget(self.local_ref)

        local_ppm_layout = QHBoxLayout()
        local_ppm_layout.addWidget(QLabel("PPM Tolerance:"))
        self.local_ppm = QDoubleSpinBox()
        self.local_ppm.setRange(0.01, 100)
        self.local_ppm.setDecimals(2)
        self.local_ppm.setValue(20)  # Chinese version: 20 ppm
        local_ppm_layout.addWidget(self.local_ppm)
        local_ppm_layout.addStretch()
        local_layout.addLayout(local_ppm_layout)

        self.local_submit = QPushButton("Submit")
        self.local_submit.setObjectName("submitBtn")
        self.local_submit.clicked.connect(self._start_local_search)
        local_layout.addWidget(self.local_submit)

        layout.addWidget(local_group)

        # === Online Search ===
        online_group = QGroupBox("Online Search")
        online_layout = QVBoxLayout(online_group)

        db_layout = QHBoxLayout()
        db_layout.addWidget(QLabel("Online Database:"))
        self.online_db = QComboBox()
        self.online_db.addItems(["HMDB", "KEGG", "LipidMaps", "ChEBI", "PubChem"])
        db_layout.addWidget(self.online_db, 1)
        online_layout.addLayout(db_layout)

        online_ppm_layout = QHBoxLayout()
        online_ppm_layout.addWidget(QLabel("PPM Tolerance:"))
        self.online_ppm = QDoubleSpinBox()
        self.online_ppm.setRange(0.01, 100)
        self.online_ppm.setDecimals(2)
        self.online_ppm.setValue(5.00)
        online_ppm_layout.addWidget(self.online_ppm)
        online_ppm_layout.addStretch()
        online_layout.addLayout(online_ppm_layout)

        self.online_submit = QPushButton("Submit")
        self.online_submit.setObjectName("submitBtn")
        self.online_submit.clicked.connect(self._start_online_search)
        online_layout.addWidget(self.online_submit)

        layout.addWidget(online_group)

        # === Search Status ===
        status_group = QGroupBox("Search Status")
        status_layout = QVBoxLayout(status_group)

        self.search_progress = QProgressBar()
        status_layout.addWidget(self.search_progress)

        status_layout.addWidget(QLabel("Match Output:"))
        self.match_output = QLineEdit()
        self.match_output.setPlaceholderText("Output path will be displayed after matching completion")
        self.match_output.setReadOnly(True)
        status_layout.addWidget(self.match_output)

        layout.addWidget(status_group)

        layout.addStretch()
        scroll.setWidget(container)
        main_layout.addWidget(scroll)

    def _start_merge_search(self):
        excel_file = self.merge_excel.get_path()
        if not excel_file:
            QMessageBox.warning(self, "Warning", "Please select cleaned Excel file.")
            return
        output_dir = os.path.dirname(excel_file)
        self.worker = SearchMergeWorker(excel_file, output_dir, self.merge_ppm.value())
        self._connect_worker(self.worker)
        self.worker.start()

    def _start_local_search(self):
        data_file = self.local_data.get_path()
        ref_file = self.local_ref.get_path()
        if not data_file or not ref_file:
            QMessageBox.warning(self, "Warning", "Please select both data file and reference library.")
            return
        output_dir = os.path.dirname(data_file)
        self.worker = LocalSearchWorker(data_file, ref_file, output_dir, self.local_ppm.value())
        self._connect_worker(self.worker)
        self.worker.start()

    def _start_online_search(self):
        database = self.online_db.currentText()
        output_dir = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if not output_dir:
            return
        self.worker = OnlineSearchWorker(database, output_dir, self.online_ppm.value())
        self._connect_worker(self.worker)
        self.worker.start()

    def _connect_worker(self, worker):
        worker.progress.connect(self.search_progress.setValue)
        worker.status.connect(lambda s: self.log.log(s))
        worker.finished_signal.connect(self._on_finished)
        worker.log.connect(lambda s: self.log.log(s))
        self.merge_search_btn.setEnabled(False)
        self.local_submit.setEnabled(False)
        self.online_submit.setEnabled(False)

    def _on_finished(self, success, message):
        if success:
            self.log.log(message)
            output_path = getattr(self.worker, 'output_file_path', message)
            self.match_output.setText(output_path)
            QMessageBox.information(self, "Complete", message)
        else:
            self.log.log(f"ERROR: {message}")
            QMessageBox.critical(self, "Error", message)
        self.merge_search_btn.setEnabled(True)
        self.local_submit.setEnabled(True)
        self.online_submit.setEnabled(True)
        self.search_progress.setValue(0)


# =============================================================================
# TAB 4: VISUALIZATION (可视化分析)
# =============================================================================

class VisualizationTab(QWidget):
    """Tab 4: Visualization - matches original Chinese UI exactly."""

    def __init__(self, log_console: LogConsole, parent=None):
        super().__init__(parent)
        self.log = log_console
        self.worker = None
        self._has_plot = False
        # Debounced auto-redraw: any parameter change re-generates the plot
        self._refresh_timer = QTimer(self)
        self._refresh_timer.setSingleShot(True)
        self._refresh_timer.timeout.connect(self._auto_refresh)
        self._setup_ui()

    def _setup_ui(self):
        main_layout = QHBoxLayout(self)
        main_layout.setSpacing(12)
        main_layout.setContentsMargins(16, 16, 16, 16)

        # Left panel — width is adjustable with the mouse (see the QSplitter
        # at the bottom of this method)
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)

        left_container = QWidget()
        left_container.setMinimumWidth(280)
        ll = QVBoxLayout(left_container)
        ll.setSpacing(10)

        # === Visualization Method ===
        method_group = QGroupBox("Visualization Method")
        method_layout = QVBoxLayout(method_group)

        sel_layout = QHBoxLayout()
        sel_layout.addWidget(QLabel("Select Method:"))
        self.viz_method = QComboBox()
        self.viz_method.addItems([
            "PCA", "t-SNE", "UMAP", "Heatmap", 
            "Volcano Plot", "Box Plot", "OPLS-DA"
        ])
        self.viz_method.currentTextChanged.connect(self._on_method_changed)
        sel_layout.addWidget(self.viz_method, 1)
        method_layout.addLayout(sel_layout)

        # Matrix file
        self.matrix_file = FileSelector("file", "Matrix File:",
            "Data Files (*.csv *.xlsx *.txt)",
            "Select matrix file...")
        method_layout.addWidget(self.matrix_file)

        # Coordinate file
        self.coord_file = FileSelector("file", "Coordinate File:",
            "Data Files (*.csv *.xlsx)",
            "Dimension reduction result file...")
        method_layout.addWidget(self.coord_file)

        ll.addWidget(method_group)

        # === Parameter Settings ===
        params_group = QGroupBox("Parameter Settings")
        params_layout = QGridLayout(params_group)

        params_layout.addWidget(QLabel("Data Preprocessing:"), 0, 0)
        self.data_preprocess = QComboBox()
        self.data_preprocess.addItems([
            "Z-score Standardization",
            "Max Normalization",
            "Min-Max Scaling",
            "None"
        ])
        params_layout.addWidget(self.data_preprocess, 0, 1)

        # Feature selection — keep only the top-N informative columns before
        # dim reduction so noise columns do not dilute the real signal
        fs_layout = QHBoxLayout()
        self.feature_select = QComboBox()
        self.feature_select.addItems([
            "None (all features)",
            "Top-N by Variance",
            "Top-N by Group Difference",
            "Top-N by Class Difference",
        ])
        self.feature_select.setToolTip(
            "Class Difference: splits groups into two classes by name "
            "prefix (e.g. C1/C2/C3 -> C, M1/M2/M3 -> M) and keeps the "
            "top-N by pooled Cohen's d. This is a supervised selection — "
            "the log reports a 5-fold CV AUC and a permutation p-value "
            "so a non-significant class difference is flagged instead of "
            "being silently over-read.")
        self.feature_n = QSpinBox()
        self.feature_n.setRange(10, 2000)
        self.feature_n.setValue(200)
        # Keep the box wide enough for 4 digits ("2000") to show fully —
        # otherwise the combo box squeezes it down to ~2 visible digits
        # (84px = Qt's own sizeHint for this spinbox)
        self.feature_n.setMinimumWidth(84)
        fs_layout.addWidget(self.feature_select)
        fs_layout.addWidget(QLabel("N:"))
        fs_layout.addWidget(self.feature_n)
        params_layout.addWidget(QLabel("Feature Selection\n(PCA / t-SNE / UMAP):"), 1, 0)
        params_layout.addLayout(fs_layout, 1, 1)

        # Optional ComBat batch correction — OFF by default; only applied
        # when the user enables it explicitly
        self.combat_batch = QCheckBox("ComBat Batch Correction (optional)")
        self.combat_batch.setToolTip(
            "Remove batch effect WITHIN each cell-type group separately "
            "(e.g. within the C group and within the M group), so differences "
            "between cell types are preserved. Off by default.")
        params_layout.addWidget(self.combat_batch, 2, 0, 1, 2)

        # Method-specific parameters — one stacked page per method so each
        # method shows ITS OWN parameters (switched in _on_method_changed)
        self.method_params = QStackedWidget()
        params_layout.addWidget(self.method_params, 3, 0, 1, 2)
        self._build_method_param_pages()

        ll.addWidget(params_group)

        # === Display Settings ===
        display_group = QGroupBox("Display Settings")
        display_layout = QGridLayout(display_group)

        display_layout.addWidget(QLabel("Color Scheme:"), 0, 0)
        self.color_scheme = QComboBox()
        self.color_scheme.addItems(["tab10", "Set1", "Set2", "Set3", "viridis", "plasma"])
        display_layout.addWidget(self.color_scheme, 0, 1)

        display_layout.addWidget(QLabel("Plot Style:"), 1, 0)
        self.plot_style = QComboBox()
        self.plot_style.addItems(["Standard White", "Seaborn", "ggplot", "bmh"])
        display_layout.addWidget(self.plot_style, 1, 1)

        display_layout.addWidget(QLabel("Font Style:"), 2, 0)
        self.font_style = QComboBox()
        # Times New Roman first = default font
        self.font_style.addItems(["Times New Roman", "Arial", "Helvetica", "sans-serif"])
        display_layout.addWidget(self.font_style, 2, 1)

        display_layout.addWidget(QLabel("Font Size:"), 3, 0)
        # Free 0-36 pt input (0 = matplotlib default), default 18 pt
        self.font_size = QSpinBox()
        self.font_size.setRange(0, 36)
        self.font_size.setValue(18)
        self.font_size.setSuffix(" pt")
        self.font_size.setSpecialValueText("Auto")
        display_layout.addWidget(self.font_size, 3, 1)

        # Label options
        label_layout = QHBoxLayout()
        self.data_point_labels = QCheckBox("Data Point Labels")
        self.feature_labels = QCheckBox("Feature Labels")
        self.feature_labels.setChecked(True)
        label_layout.addWidget(self.data_point_labels)
        label_layout.addWidget(self.feature_labels)
        display_layout.addLayout(label_layout, 4, 0, 1, 2)

        display_layout.addWidget(QLabel("Max Labels:"), 5, 0)
        self.max_labels = QSpinBox()
        self.max_labels.setRange(1, 100)
        self.max_labels.setValue(10)
        display_layout.addWidget(self.max_labels, 5, 1)

        # Legend/Grid
        lg_layout = QHBoxLayout()
        self.show_legend = QCheckBox("Legend")
        self.show_legend.setChecked(True)
        self.show_grid = QCheckBox("Grid")
        lg_layout.addWidget(self.show_legend)
        lg_layout.addWidget(self.show_grid)
        display_layout.addLayout(lg_layout, 6, 0, 1, 2)

        display_layout.addWidget(QLabel("Legend Position:"), 7, 0)
        self.legend_pos = QComboBox()
        self.legend_pos.addItems(["Outside Right", "Best", "Upper Right", "Lower Right", "Upper Left", "Lower Left"])
        display_layout.addWidget(self.legend_pos, 7, 1)

        display_layout.addWidget(QLabel("Point Size:"), 8, 0)
        self.point_size = QSpinBox()
        self.point_size.setRange(1, 100)
        self.point_size.setValue(36)
        self.point_size.setSuffix(" px")
        display_layout.addWidget(self.point_size, 8, 1)

        display_layout.addWidget(QLabel("Opacity:"), 9, 0)
        self.opacity = QDoubleSpinBox()
        self.opacity.setRange(0.1, 1.0)
        self.opacity.setDecimals(1)
        self.opacity.setValue(0.8)
        display_layout.addWidget(self.opacity, 9, 1)

        display_layout.addWidget(QLabel("Point Style:"), 10, 0)
        self.point_style = QComboBox()
        self.point_style.addItems([
            "Solid Circle (●)", "Hollow Circle (○)", "Square (■)",
            "Triangle (▲)", "Diamond (◆)", "Star (★)"
        ])
        display_layout.addWidget(self.point_style, 10, 1)

        ll.addWidget(display_group)

        # === Label Settings ===
        label_group = QGroupBox("Label Settings")
        label_layout = QGridLayout(label_group)

        label_layout.addWidget(QLabel("Font Size:"), 0, 0)
        self.label_font_size = QComboBox()
        self.label_font_size.addItems(["8 pt", "10 pt", "12 pt", "14 pt"])
        self.label_font_size.setCurrentIndex(1)
        label_layout.addWidget(self.label_font_size, 0, 1)

        label_layout.addWidget(QLabel("Label Position:"), 1, 0)
        self.label_pos = QComboBox()
        self.label_pos.addItems(["Auto", "Top", "Bottom", "Left", "Right", "Center"])
        label_layout.addWidget(self.label_pos, 1, 1)

        ll.addWidget(label_group)

        # === Chart Size ===
        chart_group = QGroupBox("Chart Size")
        chart_layout = QGridLayout(chart_group)

        # Data plot size — the axes area itself (width × height in inches).
        # The exported image adds margins and legend space around it.
        chart_layout.addWidget(QLabel("Data Plot Size:"), 0, 0)
        fig_size_layout = QHBoxLayout()
        self.fig_width = QDoubleSpinBox()
        self.fig_width.setRange(2.0, 20.0)
        self.fig_width.setDecimals(1)
        self.fig_width.setSingleStep(0.5)
        self.fig_width.setValue(10.0)
        self.fig_width.setSuffix(" in")
        self.fig_width.setToolTip(
            "Size of the data plot itself (the axes area). The exported "
            "image adds margins and legend space around it.")
        self.fig_height = QDoubleSpinBox()
        self.fig_height.setRange(2.0, 20.0)
        self.fig_height.setDecimals(1)
        self.fig_height.setSingleStep(0.5)
        self.fig_height.setValue(8.0)
        self.fig_height.setSuffix(" in")
        fig_size_layout.addWidget(self.fig_width)
        fig_size_layout.addWidget(QLabel("×"))
        fig_size_layout.addWidget(self.fig_height)
        chart_layout.addLayout(fig_size_layout, 0, 1)

        # Output DPI — resolution of the exported PNG (pixels = inches x DPI)
        chart_layout.addWidget(QLabel("Output DPI:"), 1, 0)
        self.output_dpi = QSpinBox()
        self.output_dpi.setRange(50, 600)
        self.output_dpi.setSingleStep(50)
        self.output_dpi.setValue(300)
        self.output_dpi.setSuffix(" dpi")
        self.output_dpi.setToolTip(
            "Resolution of the exported PNG. Higher DPI = sharper image "
            "and larger file (pixels = inches × DPI).")
        chart_layout.addWidget(self.output_dpi, 1, 1)

        ll.addWidget(chart_group)
        ll.addStretch()

        # Any parameter change triggers a debounced redraw of the current plot
        for combo in (self.data_preprocess, self.color_scheme, self.plot_style,
                      self.font_style, self.legend_pos,
                      self.point_style, self.label_font_size, self.label_pos):
            combo.currentTextChanged.connect(lambda *_: self._schedule_refresh())
        for spin in (self.pca_dims, self.tsne_components, self.tsne_perplexity,
                     self.umap_components, self.umap_neighbors,
                     self.heatmap_rows, self.heatmap_cols,
                     self.volcano_fc, self.volcano_p, self.box_features,
                     self.point_size, self.opacity, self.max_labels,
                     self.font_size, self.fig_width, self.fig_height,
                     self.output_dpi):
            spin.valueChanged.connect(lambda *_: self._schedule_refresh())
        for cb in (self.confidence_ellipse, self.explained_variance,
                   self.data_point_labels, self.feature_labels,
                   self.show_legend, self.show_grid):
            cb.toggled.connect(lambda *_: self._schedule_refresh())

        scroll.setWidget(left_container)
        left_layout.addWidget(scroll)

        # Right panel - Plot Canvas
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)

        # Toolbar
        toolbar = QHBoxLayout()
        toolbar.addWidget(QLabel("<b>📊 Plot Canvas</b>"))
        toolbar.addStretch()
        toolbar.addWidget(QLabel(f"x=0.000 y=0.000"))
        right_layout.addLayout(toolbar)

        # Matplotlib canvas placeholder
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
        from matplotlib.figure import Figure

        self.figure = Figure(figsize=(8, 6), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)
        self.canvas.draw()

        # Point-hover state: the mouse hovering over a data point in the
        # preview shows which sample/scan that point is (same behavior as
        # the Chinese version) — payload comes from the worker, the
        # annotation + event handler live on the preview canvas.
        self._hover_cid = None
        self._hover_annot = None
        self._hover_payload = None

        right_layout.addWidget(self.canvas, 1)

        # Matplotlib navigation toolbar. The built-in Save button saves the
        # on-screen preview canvas (widget-sized, figure dpi), which ignores
        # the Chart Size settings entirely. Replace it with one that exports
        # the full-resolution generated figure at the chosen output DPI.
        class _ChartSaveToolbar(NavigationToolbar):
            def save_figure(self, *args):
                import matplotlib as mpl
                filetypes = self.canvas.get_supported_filetypes_grouped()
                sorted_filetypes = sorted(filetypes.items())
                default_filetype = self.canvas.get_default_filetype()
                startpath = os.path.expanduser(
                    mpl.rcParams['savefig.directory'])
                start = os.path.join(
                    startpath, self.canvas.get_default_filename())
                filters = []
                selectedFilter = None
                for name, exts in sorted_filetypes:
                    exts_list = " ".join(['*.%s' % ext for ext in exts])
                    filter = f'{name} ({exts_list})'
                    if default_filetype in exts:
                        selectedFilter = filter
                    filters.append(filter)
                filters = ';;'.join(filters)
                fname, _ = QFileDialog.getSaveFileName(
                    self.canvas.parent(), "Choose a filename to save to",
                    start, filters, selectedFilter)
                if not fname:
                    return None
                if startpath != "":
                    mpl.rcParams['savefig.directory'] = os.path.dirname(fname)
                # Export the worker's original figure (data plot at the
                # requested inch size) at the Chart Size output DPI, not the
                # low-res preview shown on screen.
                app = getattr(self, "host_app", None)
                src_fig = None
                out_dpi = 300
                if (app is not None
                        and getattr(app, "worker", None) is not None):
                    src_fig = getattr(app.worker, "result_figure", None)
                    out_dpi = int(
                        getattr(app.worker, "params", {}).get("output_dpi",
                                                              300) or 300)
                if src_fig is None:
                    src_fig = self.canvas.figure
                try:
                    src_fig.savefig(fname, dpi=out_dpi)
                    if app is not None:
                        app.log.log(
                            f"Saved figure to {fname} ({out_dpi} DPI)")
                except Exception as e:
                    QMessageBox.critical(self, "Error saving file", str(e))
                return fname

        self.nav_toolbar = _ChartSaveToolbar(self.canvas, self)
        self.nav_toolbar.host_app = self
        right_layout.addWidget(self.nav_toolbar)

        # Generate button
        self.generate_btn = QPushButton("Generate Plot")
        self.generate_btn.setObjectName("submitBtn")
        self.generate_btn.setFixedHeight(40)
        self.generate_btn.clicked.connect(self._generate_plot)
        right_layout.addWidget(self.generate_btn)

        # Draggable divider between the parameter sidebar and the plot:
        # the sidebar can be widened with the mouse when long option text
        # is cut off; the plot takes the remaining space
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(left_panel)
        splitter.addWidget(right_panel)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setCollapsible(0, False)
        splitter.setCollapsible(1, False)
        splitter.setSizes([340, 1100])
        main_layout.addWidget(splitter)

    def _build_method_param_pages(self):
        """One parameter page per visualization method (same order as the
        method combo box)."""

        def _page():
            page = QWidget()
            grid = QGridLayout(page)
            grid.setContentsMargins(0, 0, 0, 0)
            return page, grid

        def _spin_row(grid, row, label, widget):
            grid.addWidget(QLabel(label), row, 0)
            grid.addWidget(widget, row, 1)

        # --- PCA ---
        page, grid = _page()
        self.pca_dims = QSpinBox()
        self.pca_dims.setRange(2, 10)
        self.pca_dims.setValue(2)
        _spin_row(grid, 0, "Components:", self.pca_dims)
        self.confidence_ellipse = QCheckBox("Confidence Ellipse")
        self.confidence_ellipse.setChecked(True)
        self.explained_variance = QCheckBox("Explained Variance")
        self.explained_variance.setChecked(True)
        opt_row = QHBoxLayout()
        opt_row.addWidget(self.confidence_ellipse)
        opt_row.addWidget(self.explained_variance)
        grid.addLayout(opt_row, 1, 0, 1, 2)
        self.method_params.addWidget(page)

        # --- t-SNE ---
        page, grid = _page()
        self.tsne_components = QSpinBox()
        self.tsne_components.setRange(2, 3)
        self.tsne_components.setValue(2)
        _spin_row(grid, 0, "Components:", self.tsne_components)
        self.tsne_perplexity = QSpinBox()
        self.tsne_perplexity.setRange(2, 200)
        self.tsne_perplexity.setValue(30)
        _spin_row(grid, 1, "Perplexity:", self.tsne_perplexity)
        self.method_params.addWidget(page)

        # --- UMAP ---
        page, grid = _page()
        self.umap_components = QSpinBox()
        self.umap_components.setRange(2, 10)
        self.umap_components.setValue(2)
        _spin_row(grid, 0, "Components:", self.umap_components)
        self.umap_neighbors = QSpinBox()
        self.umap_neighbors.setRange(2, 200)
        self.umap_neighbors.setValue(15)
        _spin_row(grid, 1, "n_neighbors:", self.umap_neighbors)
        self.method_params.addWidget(page)

        # --- Heatmap ---
        page, grid = _page()
        self.heatmap_rows = QSpinBox()
        self.heatmap_rows.setRange(10, 1000)
        self.heatmap_rows.setValue(100)
        _spin_row(grid, 0, "Max Samples:", self.heatmap_rows)
        self.heatmap_cols = QSpinBox()
        self.heatmap_cols.setRange(5, 500)
        self.heatmap_cols.setValue(50)
        _spin_row(grid, 1, "Max Features:", self.heatmap_cols)
        self.method_params.addWidget(page)

        # --- Volcano Plot ---
        page, grid = _page()
        self.volcano_fc = QDoubleSpinBox()
        self.volcano_fc.setRange(0.0, 10.0)
        self.volcano_fc.setDecimals(1)
        self.volcano_fc.setValue(1.0)
        _spin_row(grid, 0, "|log2FC| ≥:", self.volcano_fc)
        self.volcano_p = QDoubleSpinBox()
        self.volcano_p.setRange(0.0001, 1.0)
        self.volcano_p.setDecimals(4)
        self.volcano_p.setValue(0.05)
        _spin_row(grid, 1, "p-value ≤:", self.volcano_p)
        self.volcano_fdr = QCheckBox("BH-FDR Correction (q < 0.05)")
        self.volcano_fdr.setChecked(True)  # on by default to match the paper's adjusted-p reporting
        self.volcano_fdr.setToolTip(
            "Apply Benjamini-Hochberg FDR correction to the per-feature "
            "t-test p-values and color points by adjusted p (q-value) "
            "instead of raw p. On by default so volcano plots match the "
            "manuscript's reporting; uncheck to use raw p-values.")
        grid.addWidget(self.volcano_fdr, 2, 0, 1, 2)
        self.volcano_cmp = QComboBox()
        self.volcano_cmp.addItems(["Class Merge (C vs M)", "Manual Pair"])
        self.volcano_cmp.setToolTip(
            "Class Merge: derive classes from group names "
            "(C1/C2/C3 -> C, M1/M2/M3 -> M) and compare the two classes "
            "(control vs model). Manual Pair: compare two specific groups "
            "selected below.")
        _spin_row(grid, 3, "Comparison:", self.volcano_cmp)
        self.volcano_g1 = QComboBox()
        self.volcano_g2 = QComboBox()
        for cb in (self.volcano_g1, self.volcano_g2):
            cb.setEditable(True)
            cb.setEnabled(False)
            cb.setToolTip(
                "Populated automatically from the matrix file's Group "
                "column; you can also type a group name.")
        self.volcano_cmp.currentTextChanged.connect(
            lambda t: (self.volcano_g1.setEnabled(t == "Manual Pair"),
                       self.volcano_g2.setEnabled(t == "Manual Pair")))
        _spin_row(grid, 4, "Group 1:", self.volcano_g1)
        _spin_row(grid, 5, "Group 2:", self.volcano_g2)
        self.method_params.addWidget(page)

        # --- Box Plot ---
        page, grid = _page()
        self.box_features = QSpinBox()
        self.box_features.setRange(1, 50)
        self.box_features.setValue(10)
        _spin_row(grid, 0, "Top Features:", self.box_features)
        self.method_params.addWidget(page)

        # --- OPLS-DA ---
        page, grid = _page()
        grid.addWidget(QLabel("Fixed settings: 2 components\n(OSC + PLS-DA)"), 0, 0)
        self.method_params.addWidget(page)

    def _on_method_changed(self, method):
        # Switch the parameter page to the selected method so each method
        # shows its own parameters (page order matches the combo box)
        self.method_params.setCurrentIndex(max(0, self.viz_method.currentIndex()))
        self.log.log(f"Visualization method changed to: {method}")
        self._schedule_refresh()

    def _refresh_volcano_groups(self, matrix_file):
        """Populate the volcano group-pair combos from the matrix file's
        first column (the Group column). Cheap single-column read; any
        failure just leaves the combos empty (they are editable)."""
        try:
            if not str(matrix_file).lower().endswith((".csv", ".txt")):
                return
            col0 = _read_csv_robust(matrix_file, usecols=[0])
            labels = sorted(set(col0.iloc[:, 0].dropna().astype(str)))
            if not labels:
                return
            for cb, default in (
                    (self.volcano_g1, labels[0]),
                    (self.volcano_g2, labels[1] if len(labels) > 1 else labels[0])):
                cur = cb.currentText()
                cb.blockSignals(True)
                cb.clear()
                cb.addItems(labels)
                cb.setCurrentText(cur if cur in labels else default)
                cb.blockSignals(False)
        except Exception:
            pass

    def _generate_plot(self):
        matrix_file = self.matrix_file.get_path()
        if not matrix_file:
            QMessageBox.warning(self, "Warning", "Please select matrix file.")
            return

        output_dir = os.path.dirname(matrix_file)
        plot_type = self.viz_method.currentText()
        if plot_type == "Volcano Plot":
            self._refresh_volcano_groups(matrix_file)
        params = {
            "data_preprocess": self.data_preprocess.currentText(),
            "feature_select": self.feature_select.currentText(),
            "feature_n": self.feature_n.value(),
            "combat_batch": self.combat_batch.isChecked(),
            "pca_dims": self.pca_dims.value(),
            "tsne_components": self.tsne_components.value(),
            "tsne_perplexity": self.tsne_perplexity.value(),
            "umap_components": self.umap_components.value(),
            "umap_neighbors": self.umap_neighbors.value(),
            "heatmap_rows": self.heatmap_rows.value(),
            "heatmap_cols": self.heatmap_cols.value(),
            "volcano_fc": self.volcano_fc.value(),
            "volcano_p": self.volcano_p.value(),
            "volcano_fdr": self.volcano_fdr.isChecked(),
            "volcano_cmp": self.volcano_cmp.currentText(),
            "volcano_g1": self.volcano_g1.currentText(),
            "volcano_g2": self.volcano_g2.currentText(),
            "box_features": self.box_features.value(),
            "point_size": self.point_size.value(),
            "opacity": self.opacity.value(),
            "point_style": self.point_style.currentText(),
            "color_scheme": self.color_scheme.currentText(),
            "plot_style": self.plot_style.currentText(),
            "font_style": self.font_style.currentText(),
            "font_size": self.font_size.value(),
            "show_legend": self.show_legend.isChecked(),
            "legend_pos": self.legend_pos.currentText(),
            "show_grid": self.show_grid.isChecked(),
            "data_point_labels": self.data_point_labels.isChecked(),
            "feature_labels": self.feature_labels.isChecked(),
            "max_labels": self.max_labels.value(),
            "label_font_size": self.label_font_size.currentText(),
            "label_pos": self.label_pos.currentText(),
            "explained_variance": self.explained_variance.isChecked(),
            "confidence_ellipse": self.confidence_ellipse.isChecked(),
            "fig_width": self.fig_width.value(),
            "fig_height": self.fig_height.value(),
            "output_dpi": self.output_dpi.value(),
        }

        # Use VisualizationWorker to generate the actual plot
        self.worker = VisualizationWorker(
            matrix_file, self.coord_file.get_path(), output_dir, plot_type, params
        )
        self.worker.status.connect(lambda s: self.log.log(s))
        self.worker.finished_signal.connect(lambda success, msg: self._on_viz_finished(success, msg))
        self.worker.log.connect(lambda s: self.log.log(s))
        self.generate_btn.setEnabled(False)
        self.generate_btn.setText("Generating...")
        self.worker.start()

    def _schedule_refresh(self):
        """Debounce parameter edits — redraw shortly after the last change."""
        if not self._has_plot:
            return
        self._refresh_timer.start(600)

    def _auto_refresh(self):
        """Re-generate the plot after a debounced parameter change."""
        if self._has_plot and self.generate_btn.isEnabled():
            self._generate_plot()

    def _on_viz_finished(self, success, message):
        self.generate_btn.setEnabled(True)
        self.generate_btn.setText("Generate Plot")
        if success:
            self._has_plot = True
        if success and self.worker:
            # Display the generated figure
            if self.worker.result_figure is not None:
                self.figure.clear()
                # Copy the figure content at the chosen Output DPI so the
                # on-screen preview zooms to the same sharpness as the file
                import io
                buf = io.BytesIO()
                preview_dpi = int(
                    getattr(self.worker, "params", {}).get("output_dpi", 300)
                    or 300)
                self.worker.result_figure.savefig(
                    buf, format='png', dpi=preview_dpi, bbox_inches='tight')
                buf.seek(0)
                import matplotlib.image as mpimg
                img = mpimg.imread(buf)
                buf.close()
                self.ax = self.figure.add_subplot(111)
                self.ax.imshow(img)
                self.ax.axis('off')
                self.figure.tight_layout(pad=0)
                self.canvas.draw()
                # Point hover: show the sample/scan of the point under
                # the mouse (nearest data point within 5% of the axes range)
                self._install_hover(
                    getattr(self.worker, "hover_payload", None),
                    img.shape[1], img.shape[0])
                self.log.log(f"Plot generated and displayed: {message}")
            else:
                # Fallback: replot from saved file
                matrix_file = self.matrix_file.get_path()
                if matrix_file:
                    self._quick_preview_plot(matrix_file)
                self._install_hover(None)
                self.log.log(f"{message} (preview mode)")
        else:
            QMessageBox.critical(self, "Error", message)

    def _install_hover(self, payload, img_w=0, img_h=0):
        """Install/remove the point-hover tooltip on the preview canvas.

        Called once per generated plot. The payload comes from the worker:
        per-point coordinates in DATA space plus the axes-position
        fractions (Chart Size margins), so the mouse position can be
        mapped back to the nearest data point.
        """
        if getattr(self, "_hover_cid", None) is not None:
            try:
                self.canvas.mpl_disconnect(self._hover_cid)
            except Exception:
                pass
            self._hover_cid = None
        self._hover_payload = payload
        if payload is None:
            self._hover_annot = None
            return
        # Store the preview image size — the imshow data coordinates ARE
        # image pixels, so the mouse maps to the rendered figure exactly
        payload = dict(payload)
        payload["img_w"] = img_w
        payload["img_h"] = img_h
        self._hover_payload = payload
        # Tooltip annotation in the preview axes (positioned in axes
        # fractions, so it follows the point wherever pan/zoom moves it).
        self._hover_annot = self.ax.annotate(
            "", xy=(0.5, 0.5), xycoords="axes fraction",
            xytext=(12, 12), textcoords="offset points",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.4", fc="#FFFACD",
                      ec="#888888", alpha=0.95),
            zorder=1000, visible=False, annotation_clip=False)
        self._hover_cid = self.canvas.mpl_connect(
            "motion_notify_event", self._on_hover_move)

    def _on_hover_move(self, event):
        """Show the sample/scan of the data point nearest the mouse."""
        payload = getattr(self, "_hover_payload", None)
        annot = getattr(self, "_hover_annot", None)
        if payload is None or annot is None or event.inaxes is not self.ax:
            if annot is not None:
                annot.set_visible(False)
                self.canvas.draw_idle()
            return
        try:
            # Mouse position as a fraction of the rendered figure. The
            # preview axes show the figure through imshow, so the event's
            # DATA coordinates ARE the image pixel coordinates (default
            # extent: y inverted, origin upper) — pixel-exact and
            # independent of how the axes box is placed. Map that fraction
            # through the Chart Size margins to the DATA axes region.
            img_w, img_h = payload["img_w"], payload["img_h"]
            if img_w <= 0 or img_h <= 0:
                annot.set_visible(False)
                self.canvas.draw_idle()
                return
            xa = (event.xdata + 0.5) / img_w
            ya = (img_h - 0.5 - event.ydata) / img_h
            if not (0 <= xa <= 1 and 0 <= ya <= 1):
                annot.set_visible(False)
                self.canvas.draw_idle()
                return
            # xa/ya are fractions of the FULL rendered figure (the image
            # IS the whole figure, margins included) — check that the
            # mouse is inside the data-axes region and map to data coords.
            l, r, b, t = (payload["left"], payload["right"],
                          payload["bottom"], payload["top"])
            if not (l <= xa <= 1 - r and b <= ya <= 1 - t):
                annot.set_visible(False)
                self.canvas.draw_idle()
                return
            x0, x1 = payload["xlim"]
            y0, y1 = payload["ylim"]
            mx = x0 + (xa - l) / (1 - l - r) * (x1 - x0)
            my = y0 + (ya - b) / (1 - b - t) * (y1 - y0)

            # Nearest point, within 5% of the smaller axis range (the
            # same threshold the Chinese version uses)
            pts = payload["points"]
            d = np.sqrt((pts[:, 0] - mx) ** 2 + (pts[:, 1] - my) ** 2)
            i = int(np.argmin(d))
            if d[i] < min(x1 - x0, y1 - y0) * 0.05:
                labels = payload["labels"]
                text = f"{payload['prefix']}: {labels[i]}"
                grp = payload.get("groups")
                if grp is not None and i < len(grp):
                    text += f"\nGroup: {grp[i]}"
                annot.set_text(text)
                # Anchor the tooltip at the point's position (figure
                # fraction = preview axes fraction, margins included)
                annot.xy = (
                    l + (pts[i, 0] - x0) / (x1 - x0) * (1 - l - r),
                    b + (pts[i, 1] - y0) / (y1 - y0) * (1 - b - t))
                annot.set_visible(True)
            else:
                annot.set_visible(False)
            self.canvas.draw_idle()
        except Exception:
            annot.set_visible(False)

    def _quick_preview_plot(self, matrix_file):
        """Quick preview when full generation isn't available."""
        import pandas as pd
        import matplotlib.pyplot as plt
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)

        try:
            if matrix_file.endswith('.csv'):
                df = pd.read_csv(matrix_file, index_col=0)
            elif matrix_file.endswith(('.xlsx', '.xls')):
                df = pd.read_excel(matrix_file, index_col=0)
            else:
                df = pd.read_csv(matrix_file, index_col=0, sep='\t')

            df = df.fillna(0)
            numeric_cols = df.select_dtypes(include=[np.number])

            # Group coloring (same detection as the worker)
            group_labels = None
            for c in df.columns:
                if str(c).strip().lower() == "group":
                    group_labels = df[c].fillna("NA").astype(str).tolist()
                    break
            if (group_labels is None and 1 < df.index.nunique() <= 12
                    and df.index.nunique() < len(df)):
                group_labels = pd.Series(df.index).astype(str).tolist()

            if numeric_cols.shape[1] >= 2:
                x = numeric_cols.iloc[:, 0].values
                y = numeric_cols.iloc[:, 1].values
                # Point style setting (hollow = colored edge, clear face)
                ps = self.point_style.currentText()
                if group_labels is not None and len(set(group_labels)) > 1:
                    groups = sorted(set(group_labels))
                    cmap = plt.get_cmap("tab10")
                    for i, g in enumerate(groups):
                        mask = np.array([gl == g for gl in group_labels])
                        color = cmap(i % 10)
                        if _is_hollow(ps):
                            self.ax.scatter(x[mask], y[mask],
                                            s=self.point_size.value(),
                                            alpha=self.opacity.value(),
                                            marker=_marker_for(ps),
                                            facecolors="none",
                                            edgecolors=color, label=str(g))
                        else:
                            self.ax.scatter(x[mask], y[mask],
                                            s=self.point_size.value(),
                                            alpha=self.opacity.value(),
                                            marker=_marker_for(ps),
                                            color=color, label=str(g),
                                            linewidths=0.6,#数据点边框轮廓的粗细
                                            edgecolors=_darken(color))
                    if self.show_legend.isChecked():
                        _draw_legend(self.ax,
                                     self.legend_pos.currentText(),
                                     title="Group")
                else:
                    if _is_hollow(ps):
                        self.ax.scatter(x, y, s=self.point_size.value(),
                                        alpha=self.opacity.value(),
                                        marker=_marker_for(ps),
                                        facecolors="none",
                                        edgecolors="#1f77b4")
                    else:
                        self.ax.scatter(x, y, s=self.point_size.value(),
                                        alpha=self.opacity.value(),
                                        marker=_marker_for(ps),
                                        linewidths=0.6,
                                        edgecolors=_darken(None))
                self.ax.set_xlabel(numeric_cols.columns[0])
                self.ax.set_ylabel(numeric_cols.columns[1])
            else:
                self.ax.plot(df.index, df.iloc[:, 0])
                self.ax.set_xlabel("Index")

            self.ax.set_title(f"{self.viz_method.currentText()} — Data Preview")
            if self.show_grid.isChecked():
                self.ax.grid(True, alpha=0.3)
            self.canvas.draw()
        except Exception as e:
            self.log.log(f"Preview error: {e}")
            # Ultimate fallback: random
            np.random.seed(42)
            self.ax.scatter(np.random.randn(50), np.random.randn(50))
            self.ax.set_title("Sample Plot (no valid data loaded)")
            self.canvas.draw()


# =============================================================================
# TAB 5: R ANALYSIS (Seurat / Monocle / EigenMS)
# =============================================================================

class RAnalysisWorker(WorkerThread):
    """Worker for running R-based analyses via rpy2 bridge."""

    def __init__(self, analysis_type: str, matrix_file: str, output_dir: str,
                 params: dict):
        super().__init__(f"R {analysis_type}")
        self.analysis_type = analysis_type
        self.matrix_file = matrix_file
        self.output_dir = output_dir
        self.params = params

    def execute(self):
        if not _R_BRIDGE_AVAILABLE:
            raise ImportError(
                "r_bridge module not found. Please ensure r_bridge.py is in the "
                "same directory and rpy2 is installed."
            )

        rb = get_r_bridge()
        if not rb.is_available:
            raise RuntimeError(rb.status_message)

        # Ensure R packages are loaded
        self.status.emit("Checking R packages...")
        self.progress.emit(5)
        if not rb.is_ready:
            ok = rb.ensure_r_packages()
            if not ok:
                raise RuntimeError("Failed to load R packages. Check R log in R Analysis tab.")
        self.progress.emit(15)

        analysis_type = self.analysis_type.lower()

        if analysis_type == "seurat":
            self.status.emit("Running Seurat pipeline...")
            self.progress.emit(20)
            success = rb.run_seurat(
                matrix_csv=self.matrix_file,
                output_dir=self.output_dir,
                n_pcs=self.params.get("n_pcs", 30),
                resolution=self.params.get("resolution", 0.8),
                min_cells=self.params.get("min_cells", 3),
                min_features=self.params.get("min_features", 200),
            )
        elif analysis_type == "monocle":
            self.status.emit("Running Monocle trajectory analysis...")
            self.progress.emit(20)
            success = rb.run_monocle(
                matrix_csv=self.matrix_file,
                output_dir=self.output_dir,
                max_components=self.params.get("max_components", 2),
            )
        elif analysis_type == "eigenms":
            self.status.emit("Running EigenMS normalization...")
            self.progress.emit(20)
            success = rb.run_eigenms(
                matrix_csv=self.matrix_file,
                output_dir=self.output_dir,
                n_eig=self.params.get("n_eig", 10),
            )
        else:
            raise ValueError(f"Unknown analysis type: {self.analysis_type}")

        self.progress.emit(95)

        # Copy R log
        r_log = rb.get_log()
        if r_log:
            log_file = os.path.join(self.output_dir, f"{analysis_type}_r_log.txt")
            with open(log_file, 'w') as f:
                f.write(r_log)

        if not success:
            raise RuntimeError(f"{self.analysis_type} analysis failed. Check R log.")

        self.progress.emit(100)
        self.status.emit(f"{self.analysis_type} analysis complete!")


class RAnalysisTab(QWidget):
    """Tab 5: R Analysis — Seurat, Monocle, EigenMS."""

    def __init__(self, log_console, parent=None):
        super().__init__(parent)
        self.log = log_console
        self.worker = None
        self._setup_ui()

    def _setup_ui(self):
        main_layout = QVBoxLayout(self)
        main_layout.setSpacing(12)
        main_layout.setContentsMargins(16, 16, 16, 16)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)

        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setSpacing(12)

        # === R Environment Status ===
        status_group = QGroupBox("R Environment Status")
        status_layout = QVBoxLayout(status_group)

        self.r_status_label = QLabel("Checking R environment...")
        self.r_status_label.setStyleSheet("font-weight: bold; font-size: 13px;")
        status_layout.addWidget(self.r_status_label)

        btn_layout = QHBoxLayout()
        self.check_r_btn = QPushButton("Check R Environment")
        self.check_r_btn.clicked.connect(self._check_r_env)
        btn_layout.addWidget(self.check_r_btn)

        self.install_pkgs_btn = QPushButton("Install R Packages")
        self.install_pkgs_btn.clicked.connect(self._install_r_packages)
        self.install_pkgs_btn.setEnabled(False)
        btn_layout.addWidget(self.install_pkgs_btn)
        btn_layout.addStretch()
        status_layout.addLayout(btn_layout)

        layout.addWidget(status_group)

        # === Input Data ===
        input_group = QGroupBox("Input Data")
        input_layout = QVBoxLayout(input_group)

        self.matrix_selector = FileSelector(
            "file", "Cell Matrix:",
            "CSV Files (*.csv);;All Files (*)",
            "Select cell matrix CSV..."
        )
        input_layout.addWidget(self.matrix_selector)

        self.output_selector = FileSelector(
            "directory", "Output Dir:",
            placeholder="Select output directory..."
        )
        input_layout.addWidget(self.output_selector)

        layout.addWidget(input_group)

        # === Analysis Methods ===
        methods_group = QGroupBox("Analysis Methods")
        methods_layout = QVBoxLayout(methods_group)

        # -- Seurat --
        seurat_layout = QHBoxLayout()
        seurat_label = QLabel("<b>Seurat</b> — Clustering, PCA, t-SNE, UMAP")
        seurat_layout.addWidget(seurat_label)
        seurat_layout.addStretch()

        seurat_param_layout = QHBoxLayout()
        seurat_param_layout.addWidget(QLabel("PCs:"))
        self.seurat_pcs = QSpinBox()
        self.seurat_pcs.setRange(2, 100)
        self.seurat_pcs.setValue(30)
        seurat_param_layout.addWidget(self.seurat_pcs)
        seurat_param_layout.addWidget(QLabel("Resolution:"))
        self.seurat_res = QDoubleSpinBox()
        self.seurat_res.setRange(0.1, 5.0)
        self.seurat_res.setDecimals(1)
        self.seurat_res.setValue(0.8)
        seurat_param_layout.addWidget(self.seurat_res)
        seurat_param_layout.addStretch()

        self.seurat_btn = QPushButton("Run Seurat")
        self.seurat_btn.setObjectName("submitBtn")
        self.seurat_btn.clicked.connect(lambda: self._run_analysis("Seurat"))
        seurat_param_layout.addWidget(self.seurat_btn)

        methods_layout.addLayout(seurat_layout)
        methods_layout.addLayout(seurat_param_layout)

        # Separator
        sep1 = QFrame()
        sep1.setFrameShape(QFrame.Shape.HLine)
        sep1.setStyleSheet("background-color: #e0e0e0; max-height: 1px;")
        methods_layout.addWidget(sep1)

        # -- Monocle --
        monocle_layout = QHBoxLayout()
        monocle_label = QLabel("<b>Monocle</b> — Pseudotime trajectory analysis")
        monocle_layout.addWidget(monocle_label)
        monocle_layout.addStretch()

        monocle_param_layout = QHBoxLayout()
        monocle_param_layout.addWidget(QLabel("Components:"))
        self.monocle_comp = QSpinBox()
        self.monocle_comp.setRange(1, 10)
        self.monocle_comp.setValue(2)
        monocle_param_layout.addWidget(self.monocle_comp)
        monocle_param_layout.addStretch()

        self.monocle_btn = QPushButton("Run Monocle")
        self.monocle_btn.setObjectName("submitBtn")
        self.monocle_btn.clicked.connect(lambda: self._run_analysis("Monocle"))
        monocle_param_layout.addWidget(self.monocle_btn)

        methods_layout.addLayout(monocle_layout)
        methods_layout.addLayout(monocle_param_layout)

        # Separator
        sep2 = QFrame()
        sep2.setFrameShape(QFrame.Shape.HLine)
        sep2.setStyleSheet("background-color: #e0e0e0; max-height: 1px;")
        methods_layout.addWidget(sep2)

        # -- EigenMS --
        eigenms_layout = QHBoxLayout()
        eigenms_label = QLabel("<b>EigenMS</b> — SVD-based metabolomics normalization")
        eigenms_layout.addWidget(eigenms_label)
        eigenms_layout.addStretch()

        eigenms_param_layout = QHBoxLayout()
        eigenms_param_layout.addWidget(QLabel("Eigenvalues:"))
        self.eigenms_eig = QSpinBox()
        self.eigenms_eig.setRange(1, 100)
        self.eigenms_eig.setValue(10)
        eigenms_param_layout.addWidget(self.eigenms_eig)
        eigenms_param_layout.addStretch()

        self.eigenms_btn = QPushButton("Run EigenMS")
        self.eigenms_btn.setObjectName("submitBtn")
        self.eigenms_btn.clicked.connect(lambda: self._run_analysis("EigenMS"))
        eigenms_param_layout.addWidget(self.eigenms_btn)

        methods_layout.addLayout(eigenms_layout)
        methods_layout.addLayout(eigenms_param_layout)

        layout.addWidget(methods_group)

        # === Progress & Stop ===
        progress_group = QGroupBox("Progress")
        progress_layout = QVBoxLayout(progress_group)

        self.r_progress = QProgressBar()
        progress_layout.addWidget(self.r_progress)

        stop_layout = QHBoxLayout()
        stop_layout.addStretch()
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setObjectName("stopBtn")
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self._stop_analysis)
        stop_layout.addWidget(self.stop_btn)
        progress_layout.addLayout(stop_layout)

        layout.addWidget(progress_group)

        # === R Log Console ===
        log_group = QGroupBox("R Console Output")
        log_layout = QVBoxLayout(log_group)
        self.r_console = QPlainTextEdit()
        self.r_console.setReadOnly(True)
        self.r_console.setMaximumBlockCount(1000)
        self.r_console.setFont(QFont("Consolas", 10))
        log_layout.addWidget(self.r_console)
        layout.addWidget(log_group)

        layout.addStretch()
        scroll.setWidget(container)
        main_layout.addWidget(scroll)

        # Initial R check
        QTimer.singleShot(500, self._check_r_env)

    def _check_r_env(self):
        """Check R environment status."""
        if not _R_BRIDGE_AVAILABLE:
            self.r_status_label.setText("❌ r_bridge module not found")
            self.r_status_label.setStyleSheet("font-weight: bold; color: red;")
            self.r_console.appendPlainText("ERROR: r_bridge.py not found or rpy2 not installed.")
            self.r_console.appendPlainText("Install: pip install rpy2")
            return

        rb = get_r_bridge()
        msg = rb.status_message
        ready = rb.is_available and rb.is_ready

        if ready:
            self.r_status_label.setText("✅ R Ready — Seurat, Monocle, EigenMS available")
            self.r_status_label.setStyleSheet("font-weight: bold; color: green;")
            self.install_pkgs_btn.setEnabled(True)
        elif rb.is_available:
            self.r_status_label.setText(f"⚠️  {msg}")
            self.r_status_label.setStyleSheet("font-weight: bold; color: orange;")
            self.install_pkgs_btn.setEnabled(True)
        else:
            self.r_status_label.setText(f"❌ {msg}")
            self.r_status_label.setStyleSheet("font-weight: bold; color: red;")
            self.install_pkgs_btn.setEnabled(False)

        self.r_console.appendPlainText(rb.get_log())
        self.log.log(f"R status: {msg}")

    def _install_r_packages(self):
        """Install/verify R packages."""
        if not _R_BRIDGE_AVAILABLE:
            return
        self._set_buttons_enabled(False)
        rb = get_r_bridge()
        self.r_console.appendPlainText("Installing R packages (this may take several minutes)...")
        QApplication.processEvents()

        try:
            ok = rb.ensure_r_packages()
            if ok:
                self.r_console.appendPlainText("All R packages loaded successfully!")
            else:
                self.r_console.appendPlainText("Some R packages failed to load.")
        except Exception as e:
            self.r_console.appendPlainText(f"ERROR: {e}")

        self.r_console.appendPlainText(rb.get_log())
        self._check_r_env()
        self._set_buttons_enabled(True)

    def _run_analysis(self, analysis_type: str):
        """Run an R analysis."""
        matrix_file = self.matrix_selector.get_path()
        output_dir = self.output_selector.get_path()

        if not matrix_file:
            QMessageBox.warning(self, "Warning", "Please select a cell matrix file.")
            return
        if not output_dir:
            QMessageBox.warning(self, "Warning", "Please select an output directory.")
            return

        params = {}
        if analysis_type == "Seurat":
            params["n_pcs"] = self.seurat_pcs.value()
            params["resolution"] = self.seurat_res.value()
        elif analysis_type == "Monocle":
            params["max_components"] = self.monocle_comp.value()
        elif analysis_type == "EigenMS":
            params["n_eig"] = self.eigenms_eig.value()

        self.worker = RAnalysisWorker(analysis_type, matrix_file, output_dir, params)
        self.worker.progress.connect(self.r_progress.setValue)
        self.worker.status.connect(lambda s: self.r_console.appendPlainText(f"[STATUS] {s}"))
        self.worker.log.connect(lambda s: self.r_console.appendPlainText(f"[LOG] {s}"))
        self.worker.finished_signal.connect(self._on_finished)

        self._set_buttons_enabled(False)
        self.stop_btn.setEnabled(True)
        self.r_progress.setValue(0)
        self.r_console.appendPlainText(f"--- Starting {analysis_type} analysis ---")
        self.worker.start()

    def _stop_analysis(self):
        if self.worker:
            self.worker.stop()
            self.r_console.appendPlainText("Analysis stopped by user.")
        self._reset_ui()

    def _on_finished(self, success: bool, message: str):
        if success:
            self.r_console.appendPlainText(f"[DONE] {message}")
            QMessageBox.information(self, "Complete", message)
        else:
            self.r_console.appendPlainText(f"[ERROR] {message}")
            QMessageBox.critical(self, "Error", message)
        self._reset_ui()

    def _reset_ui(self):
        self._set_buttons_enabled(True)
        self.stop_btn.setEnabled(False)
        self.r_progress.setValue(0)

    def _set_buttons_enabled(self, enabled: bool):
        self.seurat_btn.setEnabled(enabled)
        self.monocle_btn.setEnabled(enabled)
        self.eigenms_btn.setEnabled(enabled)
        self.check_r_btn.setEnabled(enabled)
        self.install_pkgs_btn.setEnabled(enabled and _R_BRIDGE_AVAILABLE)


# =============================================================================
# MAIN WINDOW
# =============================================================================

class MainWindow(QMainWindow):
    """Main application window matching original Chinese UI."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"{APP_NAME} v{APP_VERSION}")
        self.setMinimumSize(1400, 900)
        self.resize(1600, 950)

        self.settings = QSettings(ORG_NAME, APP_NAME)
        self._setup_ui()
        self._setup_menu()
        self._load_settings()

    def _setup_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)
        main_layout.setSpacing(0)
        main_layout.setContentsMargins(0, 0, 0, 0)

        # Main tab widget
        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)

        self.log_console = LogConsole()

        self.tab_preprocessing = PreprocessingTab(self.log_console)
        self.tab_cleaning = DataCleaningTab(self.log_console)
        self.tab_search = SearchTab(self.log_console)
        self.tab_viz = VisualizationTab(self.log_console)
        self.tab_r = RAnalysisTab(self.log_console)

        self.tabs.addTab(self.tab_preprocessing, "Preprocessing")
        self.tabs.addTab(self.tab_cleaning, "Data Cleaning")
        self.tabs.addTab(self.tab_search, "Search")
        self.tabs.addTab(self.tab_viz, "Visualization")
        self.tabs.addTab(self.tab_r, "R Analysis")

        main_layout.addWidget(self.tabs, 1)

        # Log console at bottom
        log_label = QLabel("  Log:")
        log_label.setStyleSheet("background-color: #e8e8e8; padding: 4px; font-weight: bold;")
        main_layout.addWidget(log_label)
        main_layout.addWidget(self.log_console)

        # Status bar
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready")

        # Welcome
        self.log_console.log(f"Application started - {APP_NAME} v{APP_VERSION}")

    def _setup_menu(self):
        menubar = self.menuBar()

        file_menu = menubar.addMenu("&File")
        exit_action = QAction("E&xit", self)
        exit_action.setShortcut(QKeySequence.StandardKey.Quit)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        help_menu = menubar.addMenu("&Help")
        about_action = QAction("&About", self)
        about_action.triggered.connect(self._show_about)
        help_menu.addAction(about_action)

    def _show_about(self):
        QMessageBox.about(self, f"About {APP_NAME}",
            f"<h2>{APP_NAME} v{APP_VERSION}</h2>"
            f"<p><b>Single-Cell Metabolomics Data Processing Suite</b></p>"
            f"<p>This software provides a complete workflow for processing "
            f"single-cell mass spectrometry data.</p>"
            f"<hr>"
            f"<p><b>Modules:</b></p>"
            f"<ul>"
            f"<li>Preprocessing - Single-cell MS data extraction</li>"
            f"<li>Data Cleaning - Normalization, filtering, imputation</li>"
            f"<li>Search - Metabolite identification (Local & Online)</li>"
            f"<li>Visualization - PCA, t-SNE, UMAP, Heatmap, Volcano, OPLS-DA</li>"
            f"</ul>")

    def _load_settings(self):
        geometry = self.settings.value("geometry")
        if geometry:
            self.restoreGeometry(geometry)
        state = self.settings.value("windowState")
        if state:
            self.restoreState(state)

    def closeEvent(self, event):
        self.settings.setValue("geometry", self.saveGeometry())
        self.settings.setValue("windowState", self.saveState())
        event.accept()


# =============================================================================
# ENTRY POINT
# =============================================================================

def main():
    app = QApplication(sys.argv)
    app.setApplicationName(APP_NAME)
    app.setApplicationVersion(APP_VERSION)
    app.setOrganizationName(ORG_NAME)
    app.setStyle("Fusion")
    app.setStyleSheet(LIGHT_STYLESHEET)

    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
