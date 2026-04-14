#!/usr/bin/env python3
"""
MS Peak Inspector — Interactive mass spectrometry peak review tool.
"""

import os
import sys
import json
import csv
import itertools
import logging
import time
import traceback
from datetime import datetime
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends import _backend_tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tksheet import Sheet
from auto_selection import analyze_window


APP_DIR = os.path.dirname(os.path.abspath(__file__))
LOG_DIR = os.path.join(APP_DIR, "logs")
LOG_FILE_PATH = os.path.join(
    LOG_DIR,
    f"peak_inspector_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log",
)


def _setup_logging():
    os.makedirs(LOG_DIR, exist_ok=True)

    logger = logging.getLogger("peak_inspector")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    file_handler = logging.FileHandler(LOG_FILE_PATH, encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    logger.info("Logging started")
    logger.info("Log file: %s", LOG_FILE_PATH)
    return logger


LOGGER = _setup_logging()
MMU_PER_DA = 1000.0
SUPERSCRIPT_TRANSLATION = str.maketrans("0123456789+-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")
SUBSCRIPT_TRANSLATION = str.maketrans("0123456789+-", "₀₁₂₃₄₅₆₇₈₉₊₋")


def _log_exception(message, exc_info=None):
    if exc_info is None:
        LOGGER.exception(message)
        return

    exc_type, exc_value, exc_tb = exc_info
    formatted = "".join(traceback.format_exception(exc_type, exc_value, exc_tb))
    LOGGER.error("%s\n%s", message, formatted)


def _global_excepthook(exc_type, exc_value, exc_tb):
    _log_exception("Unhandled exception", (exc_type, exc_value, exc_tb))
    sys.__excepthook__(exc_type, exc_value, exc_tb)


sys.excepthook = _global_excepthook


def _format_duration_ms(start_time):
    return f"{(time.perf_counter() - start_time) * 1000.0:.1f} ms"


def _da_to_mmu(error_da):
    if error_da is None:
        return None
    return float(error_da) * MMU_PER_DA


def _to_superscript(value):
    return str(value).translate(SUPERSCRIPT_TRANSLATION)


def _to_subscript(value):
    return str(value).translate(SUBSCRIPT_TRANSLATION)


def _is_transient_tkagg_error(exc):
    """Return True for TkAgg redraw errors that occur during widget teardown."""
    msg = str(exc).lower()
    tokens = (
        "application has been destroyed",
        "invalid command name",
        "main thread is not in main loop",
        "pyaggimagephoto",
        "pyimage",
        "can't invoke",
    )
    return any(token in msg for token in tokens)


class SafeFigureCanvasTkAgg(FigureCanvasTkAgg):
    """TkAgg canvas that suppresses backend redraw failures during shutdown."""

    def _can_draw(self):
        try:
            widget = self.get_tk_widget()
            return bool(widget.winfo_exists())
        except tk.TclError:
            return False
        except RuntimeError:
            return False

    def draw(self):
        if not self._can_draw():
            return
        try:
            super().draw()
        except (tk.TclError, RuntimeError) as exc:
            if self._can_draw() and not _is_transient_tkagg_error(exc):
                raise

    def draw_idle(self):
        if not self._can_draw():
            return
        try:
            super().draw_idle()
        except (tk.TclError, RuntimeError) as exc:
            if self._can_draw() and not _is_transient_tkagg_error(exc):
                raise

    def blit(self, bbox=None):
        if not self._can_draw():
            return
        try:
            super().blit(bbox=bbox)
        except (tk.TclError, RuntimeError) as exc:
            if self._can_draw() and not _is_transient_tkagg_error(exc):
                raise


_ORIGINAL_TK_BLIT = _backend_tk._blit


def _safe_tk_blit(argsid):
    try:
        _ORIGINAL_TK_BLIT(argsid)
    except (tk.TclError, RuntimeError) as exc:
        if not _is_transient_tkagg_error(exc):
            raise


_backend_tk._blit = _safe_tk_blit

# ── Defaults ─────────────────────────────────────────────────────────────────
DEFAULT_CSV       = "data/samples/neodynium.csv"
DEFAULT_TARGETS   = [
157.901,
158.903,
159.904,
160.906,
161.907,
163.91,
165.914
]
DEFAULT_HALF_WIN  = 0.15
DEFAULT_VIEW_WIN  = 0.01
DEFAULT_USE_SG_SMOOTHING = True
DEFAULT_SG_WINDOW = 0
DEFAULT_SG_POLY   = 4
DEFAULT_INTENSITY_FRACTION = 0.01
DEFAULT_ISOTOPE_JSON = "data/reference/isotopes.json"
DEFAULT_MIN_ISOTOPIC_COMPOSITION = 1e-5
DEFAULT_MAX_COMBO_SIZE = 2
DEFAULT_MAX_LABEL_ROWS = 60
SUMMARY_METHODS = [
    ("manual", "Manual", "manual"),
    ("auto_active", "Auto", "auto"),
    ("derivative", "Deriv", "deriv"),
    ("gradient", "Grad", "gradient"),
    ("threshold", "Thresh", "thresh"),
    ("minima", "Minima", "minima"),
    ("combined", "Union", "combined"),
]

PERIODIC_TABLE_LAYOUT = [
    ["H", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "He"],
    ["Li", "Be", "", "", "", "", "", "", "", "", "", "", "B", "C", "N", "O", "F", "Ne"],
    ["Na", "Mg", "", "", "", "", "", "", "", "", "", "", "Al", "Si", "P", "S", "Cl", "Ar"],
    ["K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr"],
    ["Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe"],
    ["Cs", "Ba", "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn"],
    ["Fr", "Ra", "Ac", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"],
    ["", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""],
    ["", "", "", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", ""],
    ["", "", "", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", ""],
]

# ── Signal processing ─────────────────────────────────────────────────────────

def load_spectrum(path):
    """Load a two-column CSV (m/z, intensity). Handles optional header row."""
    for skip in (0, 1):
        try:
            df = pd.read_csv(path, skiprows=skip, header=None,
                             names=["x", "y"], usecols=[0, 1])
            df["x"] = pd.to_numeric(df["x"], errors="coerce")
            df["y"] = pd.to_numeric(df["y"], errors="coerce")
            df.dropna(inplace=True)
            if len(df) > 5:
                df.sort_values("x", inplace=True)
                df.reset_index(drop=True, inplace=True)
                return df
        except Exception:
            continue
    raise RuntimeError(
        "Could not parse CSV. Expected two numeric columns (m/z, intensity)."
    )


def sum_intensity_in_bounds(sub_df, left, right):
    # add y in range
    if sub_df is None or left is None or right is None:
        return None
    lo = min(float(left), float(right))
    hi = max(float(left), float(right))
    in_range = sub_df[(sub_df["x"] >= lo) & (sub_df["x"] <= hi)]
    if in_range.empty:
        return 0.0
    return float(in_range["y"].sum())


def percent_error_vs_manual(auto_val, manual_val):
    # compare to manual
    if auto_val is None or manual_val is None:
        return None
    if float(manual_val) == 0:
        return None
    return ((float(auto_val) - float(manual_val)) / float(manual_val)) * 100.0


def range_width(left, right):
    if left is None or right is None:
        return None
    return round(float(right) - float(left), 8)


def load_isotope_database(path, min_isotopic_composition=DEFAULT_MIN_ISOTOPIC_COMPOSITION):
    """Load isotope records keyed by element symbol for selection/matching."""
    if not os.path.exists(path):
        return {"by_element": {}, "elements": []}

    try:
        with open(path, "r", encoding="utf-8") as fh:
            records = json.load(fh)
    except Exception:
        return {"by_element": {}, "elements": []}

    by_element = {}
    for rec in records:
        symbol = str(rec.get("atomic_symbol", "")).strip()
        mass_number = rec.get("mass_number")
        mass_value = rec.get("relative_atomic_mass_value")
        comp_value = rec.get("isotopic_composition_value")

        try:
            comp_value = float(comp_value)
        except (TypeError, ValueError):
            comp_value = None

        if (
            not symbol
            or mass_number is None
            or mass_value is None
            or comp_value is None
            or comp_value <= float(min_isotopic_composition)
        ):
            continue

        label = rec.get("isotope_label") or f"{symbol}-{mass_number}"
        iso = {
            "symbol": symbol,
            "mass_number": int(mass_number),
            "label": label,
            "mass": float(mass_value),
            "uncertainty": rec.get("relative_atomic_mass_uncertainty"),
        }
        by_element.setdefault(symbol, []).append(iso)

    for symbol, lst in by_element.items():
        lst.sort(key=lambda x: x["mass_number"])

    return {
        "by_element": by_element,
        "elements": sorted(by_element.keys()),
    }


def format_isotope_combo_label(items):
    """Render isotope combinations using isotope-style chemical notation."""
    if not items:
        return ""

    grouped = []
    for item in items:
        symbol = str(item.get("symbol", "")).strip()
        mass_number = item.get("mass_number")
        if not symbol or mass_number is None:
            continue

        if (
            grouped
            and grouped[-1]["symbol"] == symbol
            and grouped[-1]["mass_number"] == int(mass_number)
        ):
            grouped[-1]["count"] += 1
            continue

        grouped.append({
            "symbol": symbol,
            "mass_number": int(mass_number),
            "count": 1,
        })

    if not grouped:
        return ""

    parts = []
    for group in grouped:
        symbol = group["symbol"]
        mass_number = group["mass_number"]
        count = group["count"]

        # Oxide peaks are conventionally shown as NdO rather than 16O.
        show_mass_number = not (symbol == "O" and mass_number == 16)
        part = f"{_to_superscript(mass_number)}{symbol}" if show_mass_number else symbol
        if count > 1:
            part += _to_subscript(count)
        parts.append(part)

    return "".join(parts)


def compute_isotope_label_matches(selected_isotopes, observed_mz, max_combo_size=2,
                                  max_results=60, max_combinations=80000):
    """Build candidate isotope-combination labels ranked by absolute error."""
    if not selected_isotopes or observed_mz is None:
        return [], 0, False

    n = len(selected_isotopes)
    max_combo_size = max(1, int(max_combo_size))
    tested = 0
    truncated = False
    matches = []

    for combo_size in range(1, max_combo_size + 1):
        for idx_combo in itertools.combinations_with_replacement(range(n), combo_size):
            tested += 1
            if tested > max_combinations:
                truncated = True
                break

            items = [selected_isotopes[i] for i in idx_combo]
            calc_mz = float(sum(it["mass"] for it in items))
            error_da = float(observed_mz - calc_mz)
            error_ppm = (error_da / float(observed_mz)) * 1e6 if observed_mz else None
            raw_label = " + ".join(it["label"] for it in items)
            label = format_isotope_combo_label(items) or raw_label

            matches.append({
                "label": label,
                "raw_label": raw_label,
                "calc_mz": calc_mz,
                "error_da": error_da,
                "error_mmu": _da_to_mmu(error_da),
                "error_ppm": error_ppm,
                "combo_size": combo_size,
            })

        if truncated:
            break

    matches.sort(key=lambda m: (abs(m["error_da"]), m["combo_size"], m["label"]))
    return matches[:max_results], tested, truncated


# ── Main application ──────────────────────────────────────────────────────────

class PeakInspector(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("MS Peak Inspector")
        self.geometry("1580x820")
        self.minsize(1200, 600)
        self.protocol("WM_DELETE_WINDOW", self._on_close)
        LOGGER.info("PeakInspector window initialized")

        # app state
        self.spectrum_df  = None
        self.targets      = []       # sorted list[float]
        self.current_idx  = 0
        self.results      = {}       # target -> analysis dict
        self.manual       = {}       # target -> {left, peak, right}
        self.selected_x   = None
        self._suppress_target_select = False
        self._suppress_table_select = False
        self.analysis     = None
        self.label_matches = []
        self.isotope_db = load_isotope_database(DEFAULT_ISOTOPE_JSON)
        self.selected_isotopes = []
        self._selected_isotope_keys = set()
        self.selector_window = None
        self.analysis_pref_window = None
        self.analysis_pref_listbox = None
        self.analysis_pref_summary_var = tk.StringVar(value="Select an analysis type to view/edit parameters.")
        self.analysis_pref_editor_frame = None
        self._analysis_pref_hover_name = None
        self._default_load_after_id = None
        self._startup_after_id = None
        self._closing = False
        self.selector_element_var = tk.StringVar(value="")
        self.selected_isotope_count_var = tk.StringVar(value="Selected isotopes: 0")
        self.ratio_source_var = tk.StringVar(value="combined")
        self.ratio_status_var = tk.StringVar(value="Select one or more principal peaks.")
        self.ratio_active_principal_var = tk.StringVar(value="")
        self.ratio_principals = set()
        self.ratio_secondary_map = {}
        self._ratio_active_principal_lookup = {}
        self._ratio_secondary_targets_view = []
        self._suppress_ratio_primary_select = False
        self._suppress_ratio_secondary_select = False
        self._sheet_widgets = []

        # tk variables
        self.spectrum_path_var = tk.StringVar(value=DEFAULT_CSV)
        self.half_win_var      = tk.DoubleVar(value=DEFAULT_HALF_WIN)
        self.view_win_var      = tk.DoubleVar(value=DEFAULT_VIEW_WIN)
        self.use_sg_smoothing_var = tk.BooleanVar(value=DEFAULT_USE_SG_SMOOTHING)
        self.sg_window_var     = tk.IntVar(value=DEFAULT_SG_WINDOW)
        self.sg_poly_var       = tk.IntVar(value=DEFAULT_SG_POLY)
        self.use_derivative_var = tk.BooleanVar(value=True)
        self.use_gradient_var   = tk.BooleanVar(value=True)
        self.use_threshold_var  = tk.BooleanVar(value=True)
        self.intensity_fraction_var = tk.DoubleVar(value=DEFAULT_INTENSITY_FRACTION)
        self.use_minima_var     = tk.BooleanVar(value=True)
        self.mode_var           = tk.StringVar(value="intersection")
        self.show_deriv_var    = tk.BooleanVar(value=True)
        self.new_target_var    = tk.StringVar()
        self.max_combo_var     = tk.IntVar(value=DEFAULT_MAX_COMBO_SIZE)
        self.max_label_rows_var = tk.IntVar(value=DEFAULT_MAX_LABEL_ROWS)
        self.status_var        = tk.StringVar(value="Load spectrum and targets to begin.")
        self.sel_label_var     = tk.StringVar(value="Selected x: —")

        self._build_menu()
        self._build_ui()
        self._bind_keys()
        self._set_startup_window_state()
        self._queue_default_startup_load()

    def _queue_default_startup_load(self):
        self._default_load_after_id = self.after(75, self._load_default_startup_data)

    def _load_default_startup_data(self):
        start_time = time.perf_counter()
        self._default_load_after_id = None
        if self._closing or not self.winfo_exists():
            return

        LOGGER.info("Startup defaults begin")
        self._status("Loading defaults...")
        if os.path.exists(DEFAULT_CSV):
            self._load_spectrum_file(DEFAULT_CSV)
            LOGGER.info("Startup defaults: spectrum loaded in %s", _format_duration_ms(start_time))
        for target in DEFAULT_TARGETS:
            self._add_target(target)
        LOGGER.info(
            "Startup defaults: added %d targets in %s",
            len(DEFAULT_TARGETS),
            _format_duration_ms(start_time),
        )
        if self.targets:
            self._navigate_to(0)
        LOGGER.info("Startup defaults complete in %s", _format_duration_ms(start_time))

    # ── Menu ─────────────────────────────────────────────────────────────────

    def _build_menu(self):
        mb = tk.Menu(self)
        self.config(menu=mb)
        fm = tk.Menu(mb, tearoff=0)
        mb.add_cascade(label="File", menu=fm)
        fm.add_command(label="Load Spectrum CSV…",  command=self._open_spectrum)
        fm.add_command(label="Load Targets File…",  command=self._open_targets)
        fm.add_separator()
        fm.add_command(label="Export Results CSV…", command=self._export_results)
        fm.add_separator()
        fm.add_command(label="Exit", command=self.quit)

        self.analysis_menu = tk.Menu(mb, tearoff=0)
        mb.add_cascade(label="Analysis", menu=self.analysis_menu)
        self._rebuild_analysis_menu()

    def _rebuild_analysis_menu(self):
        m = self.analysis_menu
        m.delete(0, tk.END)

        m.add_checkbutton(
            label="Savitzky-Golay smoothing",
            variable=self.use_sg_smoothing_var,
            command=lambda: self._toggle_analysis_method("Savitzky-Golay smoothing"),
        )
        m.add_separator()

        m.add_checkbutton(
            label="Derivative bounds",
            variable=self.use_derivative_var,
            command=lambda: self._toggle_analysis_method("Derivative bounds"),
        )
        m.add_checkbutton(
            label="Gradient bounds",
            variable=self.use_gradient_var,
            command=lambda: self._toggle_analysis_method("Gradient bounds"),
        )
        m.add_checkbutton(
            label="Threshold bounds",
            variable=self.use_threshold_var,
            command=lambda: self._toggle_analysis_method("Threshold bounds"),
        )
        m.add_checkbutton(
            label="Minima bounds",
            variable=self.use_minima_var,
            command=lambda: self._toggle_analysis_method("Minima bounds"),
        )
        m.add_separator()
        m.add_radiobutton(
            label="Mode: intersection",
            variable=self.mode_var,
            value="intersection",
            command=self._set_analysis_mode,
        )
        m.add_radiobutton(
            label="Mode: union",
            variable=self.mode_var,
            value="union",
            command=self._set_analysis_mode,
        )
        m.add_separator()
        m.add_command(
            label="Edit Smoothing…",
            command=lambda: self._open_analysis_preferences("Smoothing"),
        )
        m.add_command(
            label="Edit Derivative…",
            command=lambda: self._open_analysis_preferences("Derivative"),
        )
        m.add_command(
            label="Edit Gradient…",
            command=lambda: self._open_analysis_preferences("Gradient"),
        )
        m.add_command(
            label="Edit Threshold…",
            command=lambda: self._open_analysis_preferences("Threshold"),
        )
        m.add_command(
            label="Edit Minima…",
            command=lambda: self._open_analysis_preferences("Minima"),
        )
        m.add_separator()
        m.add_command(
            label="Analysis Preferences…",
            command=lambda: self._open_analysis_preferences("General"),
        )

        m.bind("<<MenuSelect>>", self._on_analysis_menu_select)

    def _analysis_preference_summaries(self):
        return {
            "General": (
                f"Win ±={self.half_win_var.get():.5g}, View ±={self.view_win_var.get():.5g}, "
                f"mode={self.mode_var.get()}"
            ),
            "Smoothing": (
                f"enabled={bool(self.use_sg_smoothing_var.get())}, "
                f"SG win={self.sg_window_var.get()}, poly={self.sg_poly_var.get()}"
            ),
            "Derivative": f"enabled={bool(self.use_derivative_var.get())}",
            "Gradient": f"enabled={bool(self.use_gradient_var.get())}",
            "Threshold": (
                f"enabled={bool(self.use_threshold_var.get())}, "
                f"fraction={self.intensity_fraction_var.get():.5g}"
            ),
            "Minima": f"enabled={bool(self.use_minima_var.get())}",
        }

    def _toggle_analysis_method(self, method_label):
        self._reanalyze()
        summaries = self._analysis_preference_summaries()
        name = {
            "Savitzky-Golay smoothing": "Smoothing",
            "Derivative bounds": "Derivative",
            "Gradient bounds": "Gradient",
            "Threshold bounds": "Threshold",
            "Minima bounds": "Minima",
        }.get(method_label)
        if name:
            self._status(f"{method_label} -> {summaries[name]}")

    def _set_analysis_mode(self):
        self._reanalyze()
        self._status(f"Mode set to {self.mode_var.get()}")

    def _on_analysis_menu_select(self, _event):
        try:
            idx = self.analysis_menu.index("active")
        except Exception:
            idx = None
        if idx is None:
            return

        try:
            label = self.analysis_menu.entrycget(idx, "label")
        except Exception:
            return

        summaries = self._analysis_preference_summaries()
        if label == "Savitzky-Golay smoothing":
            self._status(f"Savitzky-Golay smoothing: {summaries['Smoothing']}")
        elif label == "Derivative bounds":
            self._status(f"Derivative bounds: {summaries['Derivative']}")
        elif label == "Gradient bounds":
            self._status(f"Gradient bounds: {summaries['Gradient']}")
        elif label == "Threshold bounds":
            self._status(f"Threshold bounds: {summaries['Threshold']}")
        elif label == "Minima bounds":
            self._status(f"Minima bounds: {summaries['Minima']}")
        elif label.startswith("Mode:"):
            self._status(f"Current mode: {summaries['General']}")

    def _open_analysis_preferences(self, preselect="General"):
        if self.analysis_pref_window is not None and self.analysis_pref_window.winfo_exists():
            self.analysis_pref_window.lift()
            self.analysis_pref_window.focus_set()
            self._select_analysis_pref_item(preselect)
            return

        win = tk.Toplevel(self)
        win.title("Analysis Preferences")
        win.geometry("700x430")
        win.minsize(620, 360)
        self.analysis_pref_window = win
        win.protocol("WM_DELETE_WINDOW", self._close_analysis_preferences)

        root = ttk.Frame(win, padding=8)
        root.pack(fill=tk.BOTH, expand=True)

        content = ttk.PanedWindow(root, orient=tk.HORIZONTAL)
        content.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(content, padding=(4, 4))
        content.add(left, weight=0)
        ttk.Label(left, text="Analysis Types", font=("", 10, "bold")).pack(anchor=tk.W)

        lb = tk.Listbox(left, height=8, exportselection=False, activestyle="dotbox")
        for name in ("General", "Smoothing", "Derivative", "Gradient", "Threshold", "Minima"):
            lb.insert(tk.END, name)
        lb.pack(fill=tk.BOTH, expand=True, pady=(6, 0))
        lb.bind("<<ListboxSelect>>", self._on_analysis_pref_select)
        lb.bind("<Motion>", self._on_analysis_pref_hover)
        self.analysis_pref_listbox = lb

        right = ttk.Frame(content, padding=(8, 4))
        content.add(right, weight=1)
        ttk.Label(right, text="Parameters", font=("", 10, "bold")).pack(anchor=tk.W)

        summary_lbl = ttk.Label(
            right,
            textvariable=self.analysis_pref_summary_var,
            anchor=tk.W,
            justify=tk.LEFT,
            wraplength=420,
            foreground="#404040",
        )
        summary_lbl.pack(fill=tk.X, pady=(4, 8))

        editor = ttk.LabelFrame(right, text="Edit", padding=10)
        editor.pack(fill=tk.BOTH, expand=True)
        self.analysis_pref_editor_frame = editor

        btns = ttk.Frame(root)
        btns.pack(fill=tk.X, pady=(8, 0))
        ttk.Button(btns, text="Apply + Reanalyze", command=self._apply_analysis_preferences).pack(side=tk.RIGHT)
        ttk.Button(btns, text="Close", command=self._close_analysis_preferences).pack(side=tk.RIGHT, padx=(0, 6))

        self._select_analysis_pref_item(preselect)

    def _close_analysis_preferences(self):
        if self.analysis_pref_window is not None and self.analysis_pref_window.winfo_exists():
            self.analysis_pref_window.destroy()
        self.analysis_pref_window = None
        self.analysis_pref_listbox = None
        self.analysis_pref_editor_frame = None
        self._analysis_pref_hover_name = None

    def _on_analysis_pref_hover(self, event):
        if self.analysis_pref_listbox is None:
            return
        idx = self.analysis_pref_listbox.nearest(event.y)
        if idx < 0:
            return
        try:
            name = self.analysis_pref_listbox.get(idx)
        except Exception:
            return
        if name == self._analysis_pref_hover_name:
            return
        self._analysis_pref_hover_name = name
        summaries = self._analysis_preference_summaries()
        self.analysis_pref_summary_var.set(f"{name}: {summaries.get(name, '')}")

    def _on_analysis_pref_select(self, _event):
        if self.analysis_pref_listbox is None:
            return
        sel = self.analysis_pref_listbox.curselection()
        if not sel:
            return
        name = self.analysis_pref_listbox.get(sel[0])
        self._render_analysis_pref_editor(name)

    def _select_analysis_pref_item(self, name):
        if self.analysis_pref_listbox is None:
            return
        entries = [self.analysis_pref_listbox.get(i) for i in range(self.analysis_pref_listbox.size())]
        if name not in entries:
            name = "General"
        idx = entries.index(name)
        self.analysis_pref_listbox.selection_clear(0, tk.END)
        self.analysis_pref_listbox.selection_set(idx)
        self.analysis_pref_listbox.see(idx)
        self._render_analysis_pref_editor(name)

    def _render_analysis_pref_editor(self, name):
        if self.analysis_pref_editor_frame is None:
            return

        for child in self.analysis_pref_editor_frame.winfo_children():
            child.destroy()

        summaries = self._analysis_preference_summaries()
        self.analysis_pref_summary_var.set(f"{name}: {summaries.get(name, '')}")

        f = self.analysis_pref_editor_frame
        if name == "General":
            ttk.Label(f, text="Half window (target ±)").grid(row=0, column=0, sticky="w", pady=2)
            ttk.Entry(f, textvariable=self.half_win_var, width=10).grid(row=0, column=1, sticky="w", padx=8)

            ttk.Label(f, text="View window (plot ±)").grid(row=1, column=0, sticky="w", pady=2)
            ttk.Entry(f, textvariable=self.view_win_var, width=10).grid(row=1, column=1, sticky="w", padx=8)

            ttk.Label(f, text="Selection mode").grid(row=2, column=0, sticky="w", pady=2)
            mode_box = ttk.Combobox(
                f,
                textvariable=self.mode_var,
                width=12,
                values=("intersection", "union"),
                state="readonly",
            )
            mode_box.grid(row=2, column=1, sticky="w", padx=8)
            mode_box.bind("<<ComboboxSelected>>", lambda _e: self._rebuild_analysis_menu())

        elif name == "Smoothing":
            ttk.Checkbutton(
                f,
                text="Enable Savitzky-Golay smoothing",
                variable=self.use_sg_smoothing_var,
                command=self._rebuild_analysis_menu,
            ).grid(row=0, column=0, sticky="w", pady=2)

            ttk.Label(f, text="Savitzky-Golay window").grid(row=1, column=0, sticky="w", pady=(8, 2))
            ttk.Entry(f, textvariable=self.sg_window_var, width=10).grid(row=1, column=1, sticky="w", padx=8)

            ttk.Label(f, text="Savitzky-Golay polynomial").grid(row=2, column=0, sticky="w", pady=2)
            ttk.Entry(f, textvariable=self.sg_poly_var, width=10).grid(row=2, column=1, sticky="w", padx=8)

            ttk.Label(
                f,
                text="When disabled, raw y-values are used directly for peak and bound detection.",
                foreground="#555555",
                wraplength=430,
            ).grid(row=3, column=0, columnspan=2, sticky="w", pady=(8, 0))

        elif name == "Derivative":
            ttk.Checkbutton(
                f,
                text="Enable derivative bounds",
                variable=self.use_derivative_var,
                command=self._rebuild_analysis_menu,
            ).grid(row=0, column=0, sticky="w", pady=2)
            ttk.Label(
                f,
                text="Derivative bounds use nearest d2y/dx2 zero-crossings around auto peak.",
                foreground="#555555",
                wraplength=430,
            ).grid(row=1, column=0, sticky="w", pady=(4, 0))

        elif name == "Gradient":
            ttk.Checkbutton(
                f,
                text="Enable gradient bounds",
                variable=self.use_gradient_var,
                command=self._rebuild_analysis_menu,
            ).grid(row=0, column=0, sticky="w", pady=2)
            ttk.Label(
                f,
                text="Gradient bounds follow the dy/dx sign pattern around the peak and ignore isolated noisy sign flips.",
                foreground="#555555",
                wraplength=430,
            ).grid(row=1, column=0, sticky="w", pady=(4, 0))

        elif name == "Threshold":
            ttk.Checkbutton(
                f,
                text="Enable threshold bounds",
                variable=self.use_threshold_var,
                command=self._rebuild_analysis_menu,
            ).grid(row=0, column=0, sticky="w", pady=2)
            ttk.Label(f, text="Intensity fraction of auto-peak height").grid(row=1, column=0, sticky="w", pady=(8, 2))
            ttk.Entry(f, textvariable=self.intensity_fraction_var, width=10).grid(row=1, column=1, sticky="w", padx=8)

        elif name == "Minima":
            ttk.Checkbutton(
                f,
                text="Enable minima bounds",
                variable=self.use_minima_var,
                command=self._rebuild_analysis_menu,
            ).grid(row=0, column=0, sticky="w", pady=2)
            ttk.Label(
                f,
                text="Minima bounds use nearest local minima on each side of the auto peak.",
                foreground="#555555",
                wraplength=430,
            ).grid(row=1, column=0, sticky="w", pady=(4, 0))

    def _apply_analysis_preferences(self):
        self._rebuild_analysis_menu()
        self._reanalyze()
        self._status("Applied analysis preferences and reanalyzed all targets.")

    def _build_sheet(self, parent, headers, *, show_row_index=False, default_column_width=110,
                     height=None, width=None, readonly=True):
        sheet = Sheet(
            parent,
            headers=list(headers),
            data=[],
            width=width,
            height=height,
            show_row_index=show_row_index,
            show_top_left=show_row_index,
            show_header=True,
            theme="light blue",
            default_column_width=default_column_width,
            auto_resize_columns=0,
            auto_resize_row_index="empty",
            table_wrap="",
            header_wrap="c",
            paste_can_expand_x=False,
            paste_can_expand_y=False,
        )
        sheet.enable_bindings(
            "single_select",
            "drag_select",
            "ctrl_select",
            "column_select",
            "row_select",
            "copy",
            "select_all",
        )
        if readonly:
            sheet.readonly()
        sheet.pack(fill=tk.BOTH, expand=True)
        self._sheet_widgets.append(sheet)
        return sheet

    def _set_sheet_data(self, sheet, rows, headers=None):
        if headers is not None:
            sheet.headers(list(headers), redraw=False)
        sheet.set_sheet_data([list(row) for row in rows], redraw=True)

    def _set_sheet_current_row(self, sheet, row, column=0, *, as_cell=False):
        sheet.deselect("all", redraw=False)
        if row is None or row < 0 or row >= sheet.get_total_rows():
            sheet.set_refresh_timer()
            return
        if as_cell:
            sheet.select_cell(row, column, redraw=False, run_binding_func=False)
        else:
            sheet.select_row(row, redraw=False, run_binding_func=False)
        sheet.see(row, column, redraw=True)

    def _get_sheet_current_row(self, sheet):
        selected = sheet.get_currently_selected()
        if not selected:
            return None
        return selected.row

    def _get_sheet_selected_rows(self, sheet):
        rows = set(sheet.get_selected_rows(return_tuple=True))
        rows.update(row for row, _column in sheet.get_selected_cells())
        current = sheet.get_currently_selected()
        if current and current.row is not None and current.type_ != "column":
            rows.add(current.row)
        return sorted(rows)

    def _widget_belongs_to_sheet(self, widget, sheet):
        current = widget
        while current is not None:
            if current is sheet:
                return True
            current = getattr(current, "master", None)
        return False

    def _is_text_input_widget(self, widget):
        classes = {"Entry", "TEntry", "Text", "TCombobox", "Spinbox"}
        try:
            return widget.winfo_class() in classes
        except Exception:
            return False

    def _bind_sheet_arrow_keys(self, sheet, left_right_handler, up_down_handler=None):
        for widget in (sheet, sheet.MT, sheet.RI, sheet.CH, sheet.TL):
            widget.bind("<Left>", left_right_handler)
            widget.bind("<Right>", left_right_handler)
            if up_down_handler is not None:
                widget.bind("<Up>", up_down_handler)
                widget.bind("<Down>", up_down_handler)

    def _move_target_index(self, delta):
        if not self.targets:
            return None
        new_idx = max(0, min(len(self.targets) - 1, self.current_idx + int(delta)))
        if new_idx != self.current_idx:
            self._navigate_to(new_idx)
        return new_idx

    def _on_data_prev_next_key(self, event):
        if self._is_text_input_widget(event.widget):
            return None
        delta = 1 if event.keysym == "Right" else -1
        self._move_target_index(delta)
        return "break"

    def _on_list_prev_next_key(self, event):
        if self._is_text_input_widget(event.widget):
            return None
        origin_sheet = self.target_sheet if self._widget_belongs_to_sheet(event.widget, self.target_sheet) else self.summary_sheet
        delta = 1 if event.keysym == "Right" else -1
        new_idx = self._move_target_index(delta)
        if new_idx is not None:
            self._set_sheet_current_row(origin_sheet, new_idx)
            origin_sheet.focus_set()
        return "break"

    def _on_target_summary_up_down_key(self, event):
        if self._is_text_input_widget(event.widget):
            return None
        origin_sheet = self.target_sheet if self._widget_belongs_to_sheet(event.widget, self.target_sheet) else self.summary_sheet
        delta = -1 if event.keysym == "Up" else 1
        new_idx = self._move_target_index(delta)
        if new_idx is not None:
            self._set_sheet_current_row(origin_sheet, new_idx)
            origin_sheet.focus_set()
        return "break"

    def _on_data_up_down_key(self, event):
        if self._is_text_input_widget(event.widget) or self.analysis is None:
            return None
        delta = -1 if event.keysym == "Up" else 1
        selected = self.data_sheet.get_currently_selected()
        if selected and selected.row is not None:
            row = selected.row
        elif self.selected_x is not None:
            row = int(np.argmin(np.abs(self.analysis["x"] - float(self.selected_x))))
        else:
            auto_peak = self.analysis.get("auto_peak")
            if auto_peak is None:
                return "break"
            row = int(np.argmin(np.abs(self.analysis["x"] - float(auto_peak))))
        row = max(0, min(len(self.analysis["x"]) - 1, row + delta))
        self._set_selected_x(self.analysis["x"][row], sync_table=True)
        self.data_sheet.focus_set()
        return "break"

    def _route_arrow_key(self, event):
        widget = self.focus_get() or event.widget
        if widget is None or self._is_text_input_widget(widget):
            return None

        if hasattr(self, "target_sheet") and self._widget_belongs_to_sheet(widget, self.target_sheet):
            if event.keysym in ("Left", "Right"):
                return self._on_list_prev_next_key(event)
            if event.keysym in ("Up", "Down"):
                return self._on_target_summary_up_down_key(event)

        if hasattr(self, "summary_sheet") and self._widget_belongs_to_sheet(widget, self.summary_sheet):
            if event.keysym in ("Left", "Right"):
                return self._on_list_prev_next_key(event)
            if event.keysym in ("Up", "Down"):
                return self._on_target_summary_up_down_key(event)

        if hasattr(self, "data_sheet") and self._widget_belongs_to_sheet(widget, self.data_sheet):
            if event.keysym in ("Left", "Right"):
                return self._on_data_prev_next_key(event)
            if event.keysym in ("Up", "Down"):
                return self._on_data_up_down_key(event)

        if hasattr(self, "canvas_widget") and widget is self.canvas_widget:
            if event.keysym in ("Left", "Right"):
                return self._on_data_prev_next_key(event)
            if event.keysym in ("Up", "Down"):
                return self._on_data_up_down_key(event)

        return None

    def _clipboard_tsv_rows(self):
        try:
            text = self.clipboard_get()
        except tk.TclError:
            return []
        if not text:
            return []
        return list(csv.reader(text.splitlines(), delimiter="\t"))

    def _extract_first_float(self, row):
        for cell in row:
            value = str(cell).strip()
            if not value:
                continue
            try:
                return float(value)
            except ValueError:
                continue
        return None

    def _find_isotope_by_label(self, label):
        needle = str(label).strip()
        if not needle:
            return None
        for isotopes in self.isotope_db.get("by_element", {}).values():
            for isotope in isotopes:
                if isotope["label"] == needle:
                    return isotope
        return None

    def _paste_targets_from_clipboard(self):
        rows = self._clipboard_tsv_rows()
        added = 0
        for row in rows:
            value = self._extract_first_float(row)
            if value is None or value in self.targets:
                continue
            self._add_target(value)
            added += 1
        if added:
            if self.targets:
                self._navigate_to(min(self.current_idx, len(self.targets) - 1))
            self._status(f"Added {added} target(s) from clipboard TSV.")
        else:
            self._status("No valid target masses found in clipboard TSV.")
        return "break"

    def _paste_selected_isotopes_from_clipboard(self):
        rows = self._clipboard_tsv_rows()
        added = 0
        for row in rows:
            if not row:
                continue
            isotope = self._find_isotope_by_label(row[0])
            if isotope is None:
                continue
            key = isotope["label"]
            if key in self._selected_isotope_keys:
                continue
            self.selected_isotopes.append(isotope)
            self._selected_isotope_keys.add(key)
            added += 1
        if added:
            self.selected_isotopes.sort(key=lambda x: (x["symbol"], x["mass_number"]))
            self._refresh_selected_isotope_table()
            self._refresh_label_matches()
            self._refresh_plot()
            self._status(f"Added {added} isotope(s) from clipboard TSV.")
        else:
            self._status("No isotope labels from the clipboard matched the loaded isotope database.")
        return "break"

    def _on_global_paste(self, _event=None):
        widget = self.focus_get()
        if widget is None or self._is_text_input_widget(widget):
            return None
        if hasattr(self, "target_sheet") and self._widget_belongs_to_sheet(widget, self.target_sheet):
            return self._paste_targets_from_clipboard()
        if self.selector_window is not None and self.selector_window.winfo_exists():
            if hasattr(self, "selected_iso_sheet") and self._widget_belongs_to_sheet(widget, self.selected_iso_sheet):
                return self._paste_selected_isotopes_from_clipboard()
        return None

    # ── UI layout ────────────────────────────────────────────────────────────

    def _build_ui(self):
        tb = ttk.Frame(self, padding=(4, 3))
        tb.pack(side=tk.TOP, fill=tk.X)
        self._build_toolbar(tb)

        main = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        main.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(main, width=200)
        main.add(left, weight=0)
        self._build_left(left)

        center = ttk.PanedWindow(main, orient=tk.VERTICAL)
        main.add(center, weight=1)
        self._build_plot_frame(center)

        lower = ttk.PanedWindow(center, orient=tk.HORIZONTAL)
        center.add(lower, weight=2)
        self._build_table_frame(lower)

        summary = ttk.Frame(lower, width=460)
        lower.add(summary, weight=0)
        self._build_summary_frame(summary)

        ttk.Label(self, textvariable=self.status_var,
                  relief=tk.SUNKEN, anchor=tk.W, padding=(4, 1)
                  ).pack(side=tk.BOTTOM, fill=tk.X)

    def _build_toolbar(self, p):
        ttk.Label(p, text="Spectrum:").pack(side=tk.LEFT)
        ttk.Entry(p, textvariable=self.spectrum_path_var, width=28).pack(side=tk.LEFT, padx=2)
        ttk.Button(p, text="Browse…", command=self._open_spectrum).pack(side=tk.LEFT)
        ttk.Separator(p, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=8)
        ttk.Button(p, text="Reanalyze", command=self._reanalyze).pack(side=tk.LEFT, padx=6)
        ttk.Button(
            p,
            text="Analysis Preferences…",
            command=lambda: self._open_analysis_preferences("General"),
        ).pack(side=tk.LEFT)
        ttk.Checkbutton(p, text="Show derivs", variable=self.show_deriv_var,
                        command=self._refresh_plot).pack(side=tk.LEFT, padx=4)
        ttk.Separator(p, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=8)
        ttk.Button(p, text="Isotope Selector...", command=self._open_isotope_selector).pack(side=tk.LEFT)

    def _build_left(self, p):
        ttk.Label(p, text="Targets", font=("", 9, "bold")).pack(pady=(6, 2))

        lf = ttk.Frame(p)
        lf.pack(fill=tk.BOTH, expand=True, padx=4)
        self.target_sheet = self._build_sheet(
            lf,
            ("Assigned Peak", "Mass Target"),
            show_row_index=True,
            default_column_width=112,
        )
        self._bind_sheet_arrow_keys(self.target_sheet, self._on_list_prev_next_key, self._on_target_summary_up_down_key)
        self.target_sheet.bind("<<SheetSelect>>", self._on_target_select, add="+")

        af = ttk.Frame(p)
        af.pack(fill=tk.X, padx=4, pady=2)
        ttk.Entry(af, textvariable=self.new_target_var, width=10).pack(side=tk.LEFT)
        ttk.Button(af, text="Add",    command=self._add_target_from_entry).pack(side=tk.LEFT, padx=1)
        ttk.Button(af, text="Remove", command=self._remove_target).pack(side=tk.LEFT)

        ttk.Separator(p, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=6)

        for text, cmd in [
            ("Set Left  [L]",    lambda: self._assign("left")),
            ("Set Peak  [P]",    lambda: self._assign("peak")),
            ("Set Right [R]",    lambda: self._assign("right")),
            ("Clear Manual [C]", self._clear_manual),
        ]:
            ttk.Button(p, text=text, command=cmd, width=18).pack(pady=1)

        ttk.Separator(p, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=6)

        nf = ttk.Frame(p)
        nf.pack()
        ttk.Button(nf, text="◀ Prev", command=self._prev, width=7).pack(side=tk.LEFT, padx=1)
        ttk.Button(nf, text="Next ▶", command=self._next, width=7).pack(side=tk.LEFT)
        ttk.Label(p, textvariable=self.sel_label_var,
                  foreground="royalblue", wraplength=170).pack(pady=6)
        ttk.Label(p, textvariable=self.selected_isotope_count_var,
              foreground="darkgreen", wraplength=170).pack(pady=(2, 6))

    def _build_plot_frame(self, paned):
        pf = ttk.Frame(paned)
        paned.add(pf, weight=3)
        self.fig = Figure(figsize=(9, 4))
        self.fig.subplots_adjust(left=0.08, right=0.82, top=0.91, bottom=0.12)
        self.ax1 = self.fig.add_subplot(111)
        self.ax2 = self.ax1.twinx()
        self.ax3 = self.ax1.twinx()
        self.ax3.spines["right"].set_position(("axes", 1.12))
        self.ax3.set_frame_on(True)
        self.ax3.patch.set_visible(False)
        self.canvas = SafeFigureCanvasTkAgg(self.fig, master=pf)
        canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget = canvas_widget
        canvas_widget.pack(fill=tk.BOTH, expand=True)
        canvas_widget.bind("<Destroy>", self._on_canvas_destroy, add="+")
        canvas_widget.bind("<Left>", self._on_data_prev_next_key, add="+")
        canvas_widget.bind("<Right>", self._on_data_prev_next_key, add="+")
        canvas_widget.bind("<Up>", self._on_data_up_down_key, add="+")
        canvas_widget.bind("<Down>", self._on_data_up_down_key, add="+")
        self.canvas.mpl_connect("button_press_event", self._on_plot_click)
        nav = NavigationToolbar2Tk(self.canvas, pf)
        nav.update()

    def _build_table_frame(self, paned):
        tf = ttk.Frame(paned)
        paned.add(tf, weight=2)

        top = ttk.Frame(tf)
        top.pack(fill=tk.BOTH, expand=True)
        self.data_sheet = self._build_sheet(
            top,
            ("m/z", "Intensity", "Smoothed", "dy/dx", "d²y/dx²"),
            show_row_index=True,
            default_column_width=118,
        )
        self._bind_sheet_arrow_keys(self.data_sheet, self._on_data_prev_next_key, self._on_data_up_down_key)
        self.data_sheet.bind("<<SheetSelect>>", self._on_table_select, add="+")

        ttk.Separator(tf, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=3)
        bottom = ttk.Frame(tf)
        bottom.pack(fill=tk.BOTH, expand=True)

        ttk.Label(bottom, text="Candidate isotope labels (closest first)",
                  font=("", 9, "bold")).pack(anchor=tk.W, padx=2)
        self.match_sheet = self._build_sheet(
            bottom,
            ("Label", "Calc m/z", "Error (mmu)", "Error (ppm)"),
            show_row_index=True,
            default_column_width=118,
            height=180,
        )

    def _build_summary_frame(self, parent):
        container = ttk.PanedWindow(parent, orient=tk.VERTICAL)
        container.pack(fill=tk.BOTH, expand=True)

        summary_panel = ttk.Frame(container)
        container.add(summary_panel, weight=3)

        ttk.Label(
            summary_panel,
            text="Sum of Intensities",
            font=("", 10, "bold"),
        ).pack(anchor=tk.W, padx=6, pady=(6, 2))
        ttk.Label(
            summary_panel,
            text="Rows show assigned peak and mass target; remaining columns compare manual and automated bounds.",
            wraplength=430,
            foreground="#505050",
        ).pack(anchor=tk.W, padx=6, pady=(0, 6))

        frame = ttk.Frame(summary_panel)
        frame.pack(fill=tk.BOTH, expand=True, padx=4, pady=(0, 4))

        headers = ["Assigned Peak", "Mass Target"] + [label for _column_id, label, _prefix in SUMMARY_METHODS]
        self.summary_sheet = self._build_sheet(
            frame,
            headers,
            show_row_index=True,
            default_column_width=92,
        )
        self._bind_sheet_arrow_keys(self.summary_sheet, self._on_list_prev_next_key, self._on_target_summary_up_down_key)
        self.summary_sheet.bind("<<SheetSelect>>", self._on_summary_select, add="+")

        ratio_panel = ttk.Frame(container)
        container.add(ratio_panel, weight=2)
        self._build_isotope_ratio_frame(ratio_panel)

    def _build_isotope_ratio_frame(self, parent):
        ttk.Label(
            parent,
            text="Isotope Ratios",
            font=("", 10, "bold"),
        ).pack(anchor=tk.W, padx=6, pady=(6, 2))
        ttk.Label(
            parent,
            text=(
                "Select one or more principal peaks, then assign secondary peaks for each principal. "
                "Ratios are calculated from the selected sum-of-intensities column."
            ),
            wraplength=430,
            foreground="#505050",
        ).pack(anchor=tk.W, padx=6, pady=(0, 6))

        controls = ttk.Frame(parent)
        controls.pack(fill=tk.X, padx=6, pady=(0, 4))
        ttk.Label(controls, text="Intensity source").pack(side=tk.LEFT)
        self.ratio_source_combo = ttk.Combobox(
            controls,
            textvariable=self.ratio_source_var,
            values=[prefix for _column_id, _label, prefix in SUMMARY_METHODS],
            state="readonly",
            width=10,
        )
        self.ratio_source_combo.pack(side=tk.LEFT, padx=(6, 10))
        self.ratio_source_combo.bind("<<ComboboxSelected>>", lambda _event: self._refresh_isotope_ratio_panel())
        ttk.Button(controls, text="Clear groups", command=self._clear_ratio_groups).pack(side=tk.LEFT)
        ttk.Label(
            controls,
            textvariable=self.ratio_status_var,
            foreground="#2F5597",
            wraplength=190,
        ).pack(side=tk.LEFT, padx=(10, 0))

        editor = ttk.Frame(parent)
        editor.pack(fill=tk.BOTH, padx=6, pady=(0, 4))

        principal_frame = ttk.Frame(editor)
        principal_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        ttk.Label(principal_frame, text="Principal peaks").pack(anchor=tk.W)
        principal_scroll = ttk.Scrollbar(principal_frame, orient=tk.VERTICAL)
        principal_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.ratio_primary_lb = tk.Listbox(
            principal_frame,
            selectmode=tk.EXTENDED,
            exportselection=False,
            height=6,
            yscrollcommand=principal_scroll.set,
        )
        self.ratio_primary_lb.pack(fill=tk.BOTH, expand=True)
        principal_scroll.config(command=self.ratio_primary_lb.yview)
        self.ratio_primary_lb.bind("<<ListboxSelect>>", self._on_ratio_primary_select)

        secondary_frame = ttk.Frame(editor)
        secondary_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(10, 0))
        ttk.Label(secondary_frame, text="Secondary peaks for principal").pack(anchor=tk.W)
        self.ratio_active_principal_combo = ttk.Combobox(
            secondary_frame,
            textvariable=self.ratio_active_principal_var,
            values=[],
            state="readonly",
            width=12,
        )
        self.ratio_active_principal_combo.pack(fill=tk.X, pady=(0, 4))
        self.ratio_active_principal_combo.bind("<<ComboboxSelected>>", lambda _event: self._refresh_ratio_secondary_list())
        secondary_scroll = ttk.Scrollbar(secondary_frame, orient=tk.VERTICAL)
        secondary_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.ratio_secondary_lb = tk.Listbox(
            secondary_frame,
            selectmode=tk.EXTENDED,
            exportselection=False,
            height=6,
            yscrollcommand=secondary_scroll.set,
        )
        self.ratio_secondary_lb.pack(fill=tk.BOTH, expand=True)
        secondary_scroll.config(command=self.ratio_secondary_lb.yview)
        self.ratio_secondary_lb.bind("<<ListboxSelect>>", self._on_ratio_secondary_select)

        table_frame = ttk.Frame(parent)
        table_frame.pack(fill=tk.BOTH, expand=True, padx=4, pady=(0, 4))
        self.ratio_sheet = self._build_sheet(
            table_frame,
            ("Principal", "Peak", "Role", "Intensity Sum", "Ratio"),
            show_row_index=True,
            default_column_width=112,
            height=200,
        )
        self._refresh_isotope_ratio_panel()

    # ── Key bindings ─────────────────────────────────────────────────────────

    def _bind_keys(self):
        self.bind("<Right>", lambda _: self._next())
        self.bind("<Left>",  lambda _: self._prev())
        self.bind("<l>",     lambda _: self._assign("left"))
        self.bind("<p>",     lambda _: self._assign("peak"))
        self.bind("<r>",     lambda _: self._assign("right"))
        self.bind("<c>",     lambda _: self._clear_manual())
        self.bind_all("<Left>", self._route_arrow_key, add="+")
        self.bind_all("<Right>", self._route_arrow_key, add="+")
        self.bind_all("<Up>", self._route_arrow_key, add="+")
        self.bind_all("<Down>", self._route_arrow_key, add="+")
        self.bind_all("<Control-v>", self._on_global_paste, add="+")
        self.bind_all("<Control-V>", self._on_global_paste, add="+")
        self.bind_all("<Shift-Insert>", self._on_global_paste, add="+")

    # ── File I/O ─────────────────────────────────────────────────────────────

    def _open_spectrum(self):
        path = filedialog.askopenfilename(
            title="Select Spectrum CSV",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if path:
            self._load_spectrum_file(path)

    def _load_spectrum_file(self, path):
        try:
            start_time = time.perf_counter()
            LOGGER.info("Loading spectrum file: %s", path)
            self.spectrum_df = load_spectrum(path)
            self.spectrum_path_var.set(path)
            self._status(f"Loaded {len(self.spectrum_df):,} points — {os.path.basename(path)}")
            self._reanalyze()
            self._refresh_summary_table()
            LOGGER.info("Loaded spectrum file in %s", _format_duration_ms(start_time))
        except Exception as e:
            LOGGER.exception("Failed to load spectrum file: %s", path)
            messagebox.showerror("Load Error", str(e))

    def _open_targets(self):
        path = filedialog.askopenfilename(
            title="Select Targets File (.txt or .csv, one m/z per line)",
            filetypes=[("Text/CSV", "*.txt *.csv"), ("All", "*.*")])
        if not path:
            return
        try:
            LOGGER.info("Loading targets file: %s", path)
            with open(path) as fh:
                lines = fh.read().splitlines()
            added = 0
            for ln in lines:
                ln = ln.strip()
                if not ln or ln.startswith("#"):
                    continue
                try:
                    self._add_target(float(ln.split(",")[0]))
                    added += 1
                except ValueError:
                    pass
            self._status(f"Added {added} targets from file.")
            if self.targets:
                self._navigate_to(0)
        except Exception as e:
            LOGGER.exception("Failed to load targets file: %s", path)
            messagebox.showerror("Error", str(e))

    def _export_results(self):
        path = filedialog.asksaveasfilename(
            title="Save Results",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")])
        if not path:
            return

        rows = []
        method_rows = []
        for t in self.targets:
            metrics = self._get_target_metrics(t)
            ap = metrics["auto_peak"]
            mp = metrics["manual_peak"]
            al = metrics["auto_left"]
            ar = metrics["auto_right"]
            gl = metrics["grad_left"]
            gr = metrics["grad_right"]
            ml = metrics["manual_left"]
            mr = metrics["manual_right"]
            w_auto = metrics["width_auto"]
            w_gradient = metrics["width_gradient"]
            w_manual = metrics["width_manual"]
            s_manual = metrics["sum_intensity_manual"]
            s_auto = metrics["sum_intensity_auto"]
            s_gradient = metrics["sum_intensity_gradient"]
            p_auto = metrics["percent_error_auto_vs_manual"]
            p_gradient = metrics["percent_error_gradient_vs_manual"]
            r_auto = metrics["range_diff_auto_vs_manual"]
            r_gradient = metrics["range_diff_gradient_vs_manual"]

            label_match_rows, _tested, _truncated = compute_isotope_label_matches(
                self.selected_isotopes,
                ap,
                max_combo_size=self.max_combo_var.get(),
                max_results=1,
            )
            best_label = label_match_rows[0] if label_match_rows else None

            method_rows.append({
                "target_mz": t,
                "selection_type": "manual",
                "method": "manual",
                "left": ml,
                "right": mr,
                "width": w_manual,
                "intensity_sum": s_manual,
                "range_diff_vs_manual": 0.0 if w_manual is not None else None,
                "percent_error_vs_manual": 0.0 if s_manual is not None else None,
            })
            for method_name, left_v, right_v, width_v, sum_v, range_diff_v, pct_v in [
                ("auto_active", al, ar, w_auto, s_auto, r_auto, p_auto),
                ("gradient", gl, gr, w_gradient, s_gradient, r_gradient, p_gradient),
            ]:
                method_rows.append({
                    "target_mz": t,
                    "selection_type": "automated",
                    "method": method_name,
                    "left": left_v,
                    "right": right_v,
                    "width": width_v,
                    "intensity_sum": sum_v,
                    "range_diff_vs_manual": range_diff_v,
                    "percent_error_vs_manual": pct_v,
                })

            rows.append({
                "target_mz":    t,
                "peak_center":  ap,
                "auto_peak":    ap,
                "auto_left":    al,
                "auto_right":   ar,
                "grad_left":    gl,
                "grad_right":   gr,
                "manual_peak":  mp,
                "manual_left":  ml,
                "manual_right": mr,
                "delta_peak":   range_width(ap, mp),
                "delta_left":   range_width(al, ml),
                "delta_right":  range_width(ar, mr),
                "width_auto":   w_auto,
                "width_gradient": w_gradient,
                "width_manual": w_manual,
                "delta_width":  range_width(w_auto, w_manual),
                "sum_intensity_manual": s_manual,
                "sum_intensity_auto": s_auto,
                "sum_intensity_gradient": s_gradient,
                "range_diff_auto_vs_manual": r_auto,
                "range_diff_gradient_vs_manual": r_gradient,
                "percent_error_auto_vs_manual": p_auto,
                "percent_error_gradient_vs_manual": p_gradient,
                "best_isotope_label": best_label["label"] if best_label else None,
                "best_label_calc_mz": best_label["calc_mz"] if best_label else None,
                "best_label_error_mmu": best_label["error_mmu"] if best_label else None,
                "best_label_error_ppm": best_label["error_ppm"] if best_label else None,
            })

        summary_df = pd.DataFrame(rows)
        summary_df.to_csv(path, index=False)

        methods_path = os.path.splitext(path)[0] + "_methods.csv"
        methods_df = pd.DataFrame(method_rows)
        methods_df.sort_values(["target_mz", "selection_type", "method"], inplace=True)
        methods_df.to_csv(methods_path, index=False)

        self._status(
            f"Exported summary: {os.path.basename(path)} | methods: {os.path.basename(methods_path)}"
        )

    # ── Target management ─────────────────────────────────────────────────────

    def _add_target(self, value):
        start_time = time.perf_counter()
        if not isinstance(value, float):
            try:
                value = float(value)
            except (ValueError, TypeError):
                return
        if value in self.targets:
            LOGGER.info("Skipped duplicate target %.4f", value)
            return
        LOGGER.info("Adding target %.4f", value)
        self.targets.append(value)
        self.targets.sort()
        self._refresh_target_lb()
        if self.spectrum_df is not None:
            self._run_analysis(value)
        self._refresh_summary_table()
        LOGGER.info("Added target %.4f in %s", value, _format_duration_ms(start_time))

    def _add_target_from_entry(self):
        raw = self.new_target_var.get().strip()
        try:
            self._add_target(float(raw))
            self.new_target_var.set("")
        except ValueError:
            messagebox.showwarning("Invalid", f"'{raw}' is not a valid m/z number.")

    def _remove_target(self):
        if not self.targets:
            return
        selected_rows = self._get_sheet_selected_rows(self.target_sheet)
        idx = selected_rows[0] if selected_rows else self.current_idx
        if idx < 0 or idx >= len(self.targets):
            return
        t   = self.targets[idx]
        self.targets.pop(idx)
        self.results.pop(t, None)
        self.manual.pop(t, None)
        self._refresh_target_lb()
        if self.targets:
            self._navigate_to(min(idx, len(self.targets) - 1))
        else:
            self.analysis = None
            self._refresh_table()
            self._refresh_summary_table()
            self._refresh_label_matches()
            self._refresh_plot()

    def _refresh_target_lb(self):
        rows = [
            (self._get_assigned_peak_label(target), self._format_target_value(target))
            for target in self.targets
        ]
        self._set_sheet_data(self.target_sheet, rows)
        self._suppress_target_select = True
        try:
            if self.targets:
                self._set_sheet_current_row(
                    self.target_sheet,
                    min(self.current_idx, len(self.targets) - 1),
                )
        finally:
            self._suppress_target_select = False

    def _on_target_select(self, _event):
        if self._suppress_target_select:
            return
        idx = self._get_sheet_current_row(self.target_sheet)
        if idx is None:
            return
        if idx == self.current_idx:
            return
        self._navigate_to(idx)

    # ── Navigation ────────────────────────────────────────────────────────────

    def _next(self):
        if self.targets and self.current_idx < len(self.targets) - 1:
            self._navigate_to(self.current_idx + 1)

    def _prev(self):
        if self.targets and self.current_idx > 0:
            self._navigate_to(self.current_idx - 1)

    def _navigate_to(self, idx):
        start_time = time.perf_counter()
        if not self.targets or idx < 0 or idx >= len(self.targets):
            return
        LOGGER.info("Navigate begin: idx=%d target=%.4f", idx, self.targets[idx])
        self.current_idx = idx
        target = self.targets[idx]
        self._suppress_target_select = True
        try:
            self._set_sheet_current_row(self.target_sheet, idx)
        finally:
            self._suppress_target_select = False
        self.analysis = self.results.get(target)
        self._refresh_table()
        self._select_auto_peak(refresh_plot=False, focus_table=True)
        self._refresh_summary_table()
        self._refresh_label_matches()
        self._refresh_plot()
        self._status(f"Target {idx+1}/{len(self.targets)}: {self._format_target_tag(target)}")
        LOGGER.info("Navigate end: idx=%d in %s", idx, _format_duration_ms(start_time))

    def _set_startup_window_state(self):
        self.update_idletasks()
        self._startup_after_id = self.after_idle(self._maximize_main_window)

    def _on_canvas_destroy(self, _event=None):
        if not hasattr(self, "canvas"):
            return
        try:
            canvas_widget = self.canvas.get_tk_widget()
        except (tk.TclError, RuntimeError):
            self._closing = True
            return
        if _event is not None and _event.widget is not canvas_widget:
            return
        self._closing = True

    def _cancel_pending_canvas_draws(self):
        if not hasattr(self, "canvas"):
            return
        try:
            idle_draw_id = getattr(self.canvas, "_idle_draw_id", None)
            if idle_draw_id:
                self.canvas.get_tk_widget().after_cancel(idle_draw_id)
                self.canvas._idle_draw_id = None
        except (tk.TclError, RuntimeError, ValueError):
            pass

    def _on_close(self):
        if self._closing:
            return
        self._closing = True

        if self._default_load_after_id is not None:
            try:
                self.after_cancel(self._default_load_after_id)
            except tk.TclError:
                pass
            self._default_load_after_id = None

        if self._startup_after_id is not None:
            try:
                self.after_cancel(self._startup_after_id)
            except tk.TclError:
                pass
            self._startup_after_id = None

        self._cancel_pending_canvas_draws()

        try:
            self.destroy()
        except tk.TclError:
            pass

    def _maximize_main_window(self):
        self._startup_after_id = None
        if self._closing or not self.winfo_exists():
            return
        screen_w = self.winfo_screenwidth()
        screen_h = self.winfo_screenheight()

        try:
            self.deiconify()
        except tk.TclError:
            pass

        try:
            self.geometry(f"{screen_w}x{screen_h}+0+0")
            self.update_idletasks()
        except tk.TclError:
            pass

        if sys.platform.startswith("win"):
            try:
                self.state("normal")
                self.update_idletasks()
                self.geometry(f"{screen_w}x{screen_h}+0+0")
                self.state("zoomed")
                return
            except tk.TclError:
                pass

        try:
            self.state("zoomed")
            return
        except tk.TclError:
            pass

        try:
            self.attributes("-zoomed", True)
            return
        except tk.TclError:
            pass

        self.geometry(f"{screen_w}x{screen_h}+0+0")

    # ── Signal analysis ───────────────────────────────────────────────────────

    def _run_analysis(self, target):
        if self.spectrum_df is None:
            return
        start_time = time.perf_counter()
        LOGGER.info("Analysis begin: target=%.4f", target)
        self.results[target] = analyze_window(
            self.spectrum_df, target,
            half_window=self.half_win_var.get(),
            use_smoothing=self.use_sg_smoothing_var.get(),
            sg_window=self.sg_window_var.get(),
            sg_poly=self.sg_poly_var.get(),
            intensity_fraction=self.intensity_fraction_var.get(),
        )
        LOGGER.info("Analysis end: target=%.4f in %s", target, _format_duration_ms(start_time))

    def _reanalyze(self):
        start_time = time.perf_counter()
        LOGGER.info("Reanalyze begin: targets=%d", len(self.targets))
        if self.spectrum_df is None:
            self._refresh_summary_table()
            LOGGER.info("Reanalyze skipped: no spectrum")
            return
        for t in self.targets:
            self._run_analysis(t)
        if self.targets:
            self.analysis = self.results.get(self.targets[self.current_idx])
            self._refresh_table()
            if self.selected_x is None and self.analysis:
                self._select_auto_peak(refresh_plot=False)
            self._refresh_summary_table()
            self._refresh_label_matches()
            self._refresh_plot()
        else:
            self._refresh_summary_table()
        LOGGER.info("Reanalyze end in %s", _format_duration_ms(start_time))

    def _refresh_label_matches(self):
        start_time = time.perf_counter()
        LOGGER.info("Refresh label matches begin")
        self.label_matches = []
        if not self.targets or not self.analysis or not self.selected_isotopes:
            self._refresh_match_table()
            LOGGER.info("Refresh label matches end: empty context in %s", _format_duration_ms(start_time))
            return

        observed = self.analysis.get("auto_peak")
        if observed is None:
            self._refresh_match_table()
            LOGGER.info("Refresh label matches end: no observed peak in %s", _format_duration_ms(start_time))
            return

        matches, tested, truncated = compute_isotope_label_matches(
            self.selected_isotopes,
            observed,
            max_combo_size=self.max_combo_var.get(),
            max_results=self.max_label_rows_var.get(),
        )
        self.label_matches = matches
        self._refresh_match_table()

        if truncated:
            self._status(
                f"Label candidates truncated after {tested:,} combinations. Reduce isotopes or max combo."
            )
        LOGGER.info(
            "Refresh label matches end: rows=%d tested=%d truncated=%s in %s",
            len(self.label_matches),
            tested,
            truncated,
            _format_duration_ms(start_time),
        )

    def _refresh_match_table(self):
        rows = [
            (
                row["label"],
                f"{row['calc_mz']:.8f}",
                f"{row['error_mmu']:+.5f}",
                f"{row['error_ppm']:+.2f}",
            )
            for row in self.label_matches
        ]
        self._set_sheet_data(self.match_sheet, rows)

    # ── Isotope selector ─────────────────────────────────────────────────────

    def _close_isotope_selector(self):
        if self.selector_window is not None and self.selector_window.winfo_exists():
            self.selector_window.destroy()
        self.selector_window = None

    def _fit_selector_window_to_content(self, win):
        win.update_idletasks()

        virtual_x = win.winfo_vrootx()
        virtual_y = win.winfo_vrooty()
        virtual_w = max(win.winfo_vrootwidth(), 1)
        virtual_h = max(win.winfo_vrootheight(), 1)
        max_w = max(720, int(virtual_w * 0.94))
        max_h = max(520, int(virtual_h * 0.92))

        req_w = win.winfo_reqwidth() + 24
        req_h = win.winfo_reqheight() + 24
        width = min(max(980, req_w), max_w)
        height = min(max(600, req_h), max_h)
        min_w = min(max(760, min(req_w, width)), max_w)
        min_h = min(max(480, min(req_h, height)), max_h)

        try:
            self.update_idletasks()
            base_x = self.winfo_rootx()
            base_y = self.winfo_rooty()
            base_w = self.winfo_width()
            base_h = self.winfo_height()
        except tk.TclError:
            base_x = virtual_x
            base_y = virtual_y
            base_w = virtual_w
            base_h = virtual_h

        x = base_x + max((base_w - width) // 2, 0)
        y = base_y + max((base_h - height) // 2, 0)
        x = max(virtual_x, min(x, (virtual_x + virtual_w) - width))
        y = max(virtual_y, min(y, (virtual_y + virtual_h) - height))

        win.minsize(min_w, min_h)
        win.geometry(f"{width}x{height}+{x}+{y}")

    def _open_isotope_selector(self):
        if self.selector_window is not None and self.selector_window.winfo_exists():
            try:
                self.selector_window.transient(self)
            except tk.TclError:
                pass
            self.selector_window.lift()
            self.selector_window.focus_force()
            return

        win = tk.Toplevel(self)
        self.selector_window = win
        win.title("Periodic Table Isotope Selector")
        win.protocol("WM_DELETE_WINDOW", self._close_isotope_selector)
        try:
            win.transient(self)
        except tk.TclError:
            pass

        root = ttk.Frame(win, padding=8)
        root.pack(fill=tk.BOTH, expand=True)

        left = ttk.LabelFrame(root, text="Elements", padding=6)
        left.pack(side=tk.LEFT, fill=tk.BOTH, expand=False)
        right = ttk.Frame(root)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(8, 0))

        self._build_periodic_grid(left)
        self._build_selector_right_panel(right)

        default_symbol = self.selector_element_var.get()
        if not default_symbol:
            symbols = self.isotope_db.get("elements", [])
            default_symbol = symbols[0] if symbols else ""
        self._set_selector_element(default_symbol)
        self._refresh_selected_isotope_table()
        self._fit_selector_window_to_content(win)
        win.deiconify()
        win.lift()
        win.focus_force()

    def _build_periodic_grid(self, parent):
        available = set(self.isotope_db.get("elements", []))
        self._element_buttons = {}

        for r, row in enumerate(PERIODIC_TABLE_LAYOUT):
            for c, symbol in enumerate(row):
                if not symbol:
                    ttk.Label(parent, text="", width=4).grid(row=r, column=c, padx=1, pady=1)
                    continue

                state = tk.NORMAL if symbol in available else tk.DISABLED
                btn = tk.Button(parent, text=symbol, width=4, relief=tk.RIDGE,
                                state=state,
                                command=lambda s=symbol: self._set_selector_element(s))
                btn.grid(row=r, column=c, padx=1, pady=1, sticky="nsew")
                self._element_buttons[symbol] = btn

    def _build_selector_right_panel(self, parent):
        top = ttk.Frame(parent)
        top.pack(fill=tk.BOTH, expand=True)

        title = ttk.Frame(top)
        title.pack(fill=tk.X)
        ttk.Label(title, text="Element:").pack(side=tk.LEFT)
        ttk.Label(title, textvariable=self.selector_element_var,
                  font=("", 10, "bold")).pack(side=tk.LEFT, padx=(4, 0))

        lists = ttk.PanedWindow(top, orient=tk.HORIZONTAL)
        lists.pack(fill=tk.BOTH, expand=True, pady=(4, 6))

        left = ttk.Frame(lists)
        lists.add(left, weight=1)
        ttk.Label(left, text="Isotopes for selected element").pack(anchor=tk.W)
        lb_frame = ttk.Frame(left)
        lb_frame.pack(fill=tk.BOTH, expand=True)
        self.element_iso_lb = tk.Listbox(lb_frame, selectmode=tk.EXTENDED, exportselection=False)
        self.element_iso_lb.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        lb_scroll = ttk.Scrollbar(lb_frame, orient=tk.VERTICAL, command=self.element_iso_lb.yview)
        lb_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.element_iso_lb.config(yscrollcommand=lb_scroll.set)

        right = ttk.Frame(lists)
        lists.add(right, weight=1)
        ttk.Label(right, text="Selected isotopes for matching").pack(anchor=tk.W)
        self.selected_iso_sheet = self._build_sheet(
            right,
            ("Isotope", "Mass"),
            show_row_index=True,
            default_column_width=160,
            height=220,
        )

        controls = ttk.Frame(parent)
        controls.pack(fill=tk.X)
        ttk.Button(controls, text="Add selected isotopes ->", command=self._add_isotopes_from_element_list).pack(side=tk.LEFT)
        ttk.Button(controls, text="<- Remove highlighted", command=self._remove_selected_isotopes).pack(side=tk.LEFT, padx=4)
        ttk.Button(controls, text="Clear all", command=self._clear_all_selected_isotopes).pack(side=tk.LEFT)
        ttk.Button(controls, text="Close", command=self._close_isotope_selector).pack(side=tk.RIGHT)

        settings = ttk.Frame(parent)
        settings.pack(fill=tk.X, pady=(6, 0))
        ttk.Label(settings, text="Matching settings:").pack(side=tk.LEFT)
        ttk.Label(settings, text="Max combo").pack(side=tk.LEFT, padx=(8, 1))
        ttk.Entry(settings, textvariable=self.max_combo_var, width=3).pack(side=tk.LEFT)
        ttk.Label(settings, text="Rows").pack(side=tk.LEFT, padx=(8, 1))
        ttk.Entry(settings, textvariable=self.max_label_rows_var, width=4).pack(side=tk.LEFT)
        ttk.Button(settings, text="Refresh labels", command=self._refresh_label_matches).pack(side=tk.LEFT, padx=(8, 0))

    def _set_selector_element(self, symbol):
        self.selector_element_var.set(symbol)
        if not hasattr(self, "element_iso_lb"):
            return

        self.element_iso_lb.delete(0, tk.END)
        for iso in self.isotope_db.get("by_element", {}).get(symbol, []):
            self.element_iso_lb.insert(tk.END, f"{iso['label']} ({iso['mass']:.10f})")

        for sym, btn in getattr(self, "_element_buttons", {}).items():
            if sym == symbol:
                btn.configure(relief=tk.SUNKEN, bg="#d9ead3")
            else:
                btn.configure(relief=tk.RIDGE, bg="SystemButtonFace")

    def _add_isotopes_from_element_list(self):
        symbol = self.selector_element_var.get().strip()
        if not symbol:
            return
        indices = self.element_iso_lb.curselection() if hasattr(self, "element_iso_lb") else ()
        if not indices:
            return

        source = self.isotope_db.get("by_element", {}).get(symbol, [])
        added = 0
        for i in indices:
            if i < 0 or i >= len(source):
                continue
            iso = source[i]
            key = iso["label"]
            if key in self._selected_isotope_keys:
                continue
            self.selected_isotopes.append(iso)
            self._selected_isotope_keys.add(key)
            added += 1

        if added:
            self.selected_isotopes.sort(key=lambda x: (x["symbol"], x["mass_number"]))
            self._refresh_selected_isotope_table()
            self._refresh_label_matches()
            self._refresh_plot()

    def _remove_selected_isotopes(self):
        selected_rows = self._get_sheet_selected_rows(self.selected_iso_sheet)
        if not selected_rows:
            return
        remove_keys = {
            self.selected_isotopes[row]["label"]
            for row in selected_rows
            if 0 <= row < len(self.selected_isotopes)
        }
        self.selected_isotopes = [i for i in self.selected_isotopes if i["label"] not in remove_keys]
        self._selected_isotope_keys = {i["label"] for i in self.selected_isotopes}
        self._refresh_selected_isotope_table()
        self._refresh_label_matches()
        self._refresh_plot()

    def _clear_all_selected_isotopes(self):
        self.selected_isotopes = []
        self._selected_isotope_keys = set()
        self._refresh_selected_isotope_table()
        self._refresh_label_matches()
        self._refresh_plot()

    def _refresh_selected_isotope_table(self):
        self.selected_isotope_count_var.set(f"Selected isotopes: {len(self.selected_isotopes)}")
        if not hasattr(self, "selected_iso_sheet"):
            return
        rows = [(iso["label"], f"{iso['mass']:.10f}") for iso in self.selected_isotopes]
        self._set_sheet_data(self.selected_iso_sheet, rows)

    # ── Table ─────────────────────────────────────────────────────────────────

    def _refresh_table(self):
        start_time = time.perf_counter()
        LOGGER.info("Refresh data table begin")
        rows = []
        if not self.analysis:
            self._set_sheet_data(self.data_sheet, rows)
            LOGGER.info("Refresh data table end: no analysis in %s", _format_duration_ms(start_time))
            return
        sub   = self.analysis["sub"]
        ap_mz = self.analysis["auto_peak"]
        for _, row in sub.iterrows():
            rows.append((
                f"{row['x']:.6f}",
                f"{row['y']:.4f}",
                f"{row['smoothed']:.4f}",
                f"{row['dy']:.6f}",
                f"{row['d2y']:.6f}",
            ))
        self._set_sheet_data(self.data_sheet, rows)
        peak_i = int(np.argmin(np.abs(self.analysis["x"] - ap_mz)))
        if rows and peak_i < len(rows):
            self.data_sheet.see(peak_i, 0)
        LOGGER.info("Refresh data table end: rows=%d in %s", len(rows), _format_duration_ms(start_time))

    def _get_target_metrics(self, target):
        analysis = self.results.get(target)
        manual = self.manual.get(target, {})
        sub = analysis["sub"] if analysis else None

        metrics = {
            "target_mz": target,
            "auto_peak": analysis.get("auto_peak") if analysis else None,
            "auto_left": analysis.get("auto_left") if analysis else None,
            "auto_right": analysis.get("auto_right") if analysis else None,
            "grad_left": analysis.get("grad_left") if analysis else None,
            "grad_right": analysis.get("grad_right") if analysis else None,
            "manual_peak": manual.get("peak"),
            "manual_left": manual.get("left"),
            "manual_right": manual.get("right"),
        }

        metrics["width_manual"] = range_width(metrics["manual_left"], metrics["manual_right"])
        metrics["width_auto"] = range_width(metrics["auto_left"], metrics["auto_right"])
        metrics["width_gradient"] = range_width(metrics["grad_left"], metrics["grad_right"])

        metrics["sum_intensity_manual"] = sum_intensity_in_bounds(sub, metrics["manual_left"], metrics["manual_right"])
        metrics["sum_intensity_auto"] = sum_intensity_in_bounds(sub, metrics["auto_left"], metrics["auto_right"])
        metrics["sum_intensity_gradient"] = sum_intensity_in_bounds(sub, metrics["grad_left"], metrics["grad_right"])

        manual_sum = metrics["sum_intensity_manual"]
        metrics["range_diff_auto_vs_manual"] = range_width(metrics["width_manual"], metrics["width_auto"])
        metrics["range_diff_gradient_vs_manual"] = range_width(metrics["width_manual"], metrics["width_gradient"])
        metrics["percent_error_auto_vs_manual"] = percent_error_vs_manual(metrics["sum_intensity_auto"], manual_sum)
        metrics["percent_error_gradient_vs_manual"] = percent_error_vs_manual(metrics["sum_intensity_gradient"], manual_sum)
        return metrics

    def _format_summary_value(self, value):
        if value is None:
            return "—"
        return f"{value:.4f}"

    def _format_ratio_value(self, value):
        if value is None:
            return "—"
        return f"{value:.6f}"

    def _format_target_value(self, value):
        if value is None:
            return "—"
        return f"{float(value):.4f}"

    def _get_assigned_peak_observed_mz(self, target):
        manual_peak = self.manual.get(target, {}).get("peak")
        if manual_peak is not None:
            return float(manual_peak)

        analysis = self.results.get(target)
        if analysis is None:
            return None
        auto_peak = analysis.get("auto_peak")
        if auto_peak is None:
            return None
        return float(auto_peak)

    def _get_assigned_peak_label(self, target):
        observed_mz = self._get_assigned_peak_observed_mz(target)
        if observed_mz is None:
            return "—"
        if not self.selected_isotopes:
            return self._format_target_value(observed_mz)

        matches, _tested, _truncated = compute_isotope_label_matches(
            self.selected_isotopes,
            observed_mz,
            max_combo_size=self.max_combo_var.get(),
            max_results=1,
        )
        if matches:
            return matches[0]["label"]
        return self._format_target_value(observed_mz)

    def _target_iid(self, target):
        return f"target::{float(target):.10f}"

    def _format_target_tag(self, target):
        return (
            f"Peak {self._get_assigned_peak_label(target)} | "
            f"Target {self._format_target_value(target)}"
        )

    def _clear_ratio_groups(self):
        self.ratio_principals.clear()
        self.ratio_secondary_map.clear()
        self.ratio_active_principal_var.set("")
        self._refresh_isotope_ratio_panel()

    def _sync_ratio_groups_with_targets(self):
        valid_targets = {float(target) for target in self.targets}
        self.ratio_principals = {target for target in self.ratio_principals if target in valid_targets}

        filtered_map = {}
        for principal, secondary_targets in self.ratio_secondary_map.items():
            if principal not in self.ratio_principals:
                continue
            filtered_map[principal] = [
                target for target in secondary_targets
                if target in valid_targets and target not in self.ratio_principals and target != principal
            ]
        self.ratio_secondary_map = filtered_map

    def _get_active_ratio_principal(self):
        return self._ratio_active_principal_lookup.get(self.ratio_active_principal_var.get())

    def _refresh_ratio_primary_list(self):
        if not hasattr(self, "ratio_primary_lb"):
            return
        self._suppress_ratio_primary_select = True
        try:
            self.ratio_primary_lb.delete(0, tk.END)
            for index, target in enumerate(self.targets):
                self.ratio_primary_lb.insert(tk.END, self._get_assigned_peak_label(target))
                if target in self.ratio_principals:
                    self.ratio_primary_lb.selection_set(index)
        finally:
            self._suppress_ratio_primary_select = False

    def _refresh_ratio_active_principal_options(self):
        if not hasattr(self, "ratio_active_principal_combo"):
            return
        values = []
        lookup = {}
        for target in sorted(self.ratio_principals):
            display = self._get_assigned_peak_label(target)
            if display in lookup:
                display = f"{display} [{self._format_target_value(target)}]"
            values.append(display)
            lookup[display] = target
        self._ratio_active_principal_lookup = lookup
        self.ratio_active_principal_combo.configure(values=values)
        current = self.ratio_active_principal_var.get()
        if current not in lookup:
            self.ratio_active_principal_var.set(values[0] if values else "")

    def _refresh_ratio_secondary_list(self):
        if not hasattr(self, "ratio_secondary_lb"):
            return
        active_principal = self._get_active_ratio_principal()
        available_targets = [
            target for target in self.targets
            if target not in self.ratio_principals and target != active_principal
        ]
        selected_targets = set(self.ratio_secondary_map.get(active_principal, [])) if active_principal is not None else set()

        self._ratio_secondary_targets_view = available_targets
        self._suppress_ratio_secondary_select = True
        try:
            self.ratio_secondary_lb.delete(0, tk.END)
            for index, target in enumerate(available_targets):
                self.ratio_secondary_lb.insert(tk.END, self._get_assigned_peak_label(target))
                if target in selected_targets:
                    self.ratio_secondary_lb.selection_set(index)
        finally:
            self._suppress_ratio_secondary_select = False
        self._refresh_isotope_ratio_table()
        self._refresh_ratio_status()

    def _refresh_ratio_status(self):
        if not self.ratio_principals:
            self.ratio_status_var.set("Select one or more principal peaks.")
            return
        missing = [
            self._get_assigned_peak_label(target)
            for target in sorted(self.ratio_principals)
            if not self.ratio_secondary_map.get(target)
        ]
        if missing:
            self.ratio_status_var.set(f"Assign secondaries for: {', '.join(missing)}")
            return
        self.ratio_status_var.set(f"Configured {len(self.ratio_principals)} principal peak(s).")

    def _compute_isotope_ratio_rows(self):
        rows = []
        prefix = self.ratio_source_var.get() or "combined"
        for principal in sorted(self.ratio_principals):
            principal_metrics = self._get_target_metrics(principal)
            principal_sum = principal_metrics.get(f"sum_intensity_{prefix}")
            rows.append({
                "principal": principal,
                "principal_label": self._get_assigned_peak_label(principal),
                "peak_label": self._get_assigned_peak_label(principal),
                "role": "Principal",
                "intensity_sum": principal_sum,
                "ratio": 1.0 if principal_sum is not None else None,
            })
            for peak in self.ratio_secondary_map.get(principal, []):
                peak_metrics = self._get_target_metrics(peak)
                peak_sum = peak_metrics.get(f"sum_intensity_{prefix}")
                ratio = None
                if peak_sum is not None and principal_sum not in (None, 0):
                    ratio = peak_sum / principal_sum
                rows.append({
                    "principal": principal,
                    "principal_label": self._get_assigned_peak_label(principal),
                    "peak_label": self._get_assigned_peak_label(peak),
                    "role": "Secondary",
                    "intensity_sum": peak_sum,
                    "ratio": ratio,
                })
        return rows

    def _refresh_isotope_ratio_table(self):
        if not hasattr(self, "ratio_sheet"):
            return
        rows = [
            (
                row["principal_label"],
                row["peak_label"],
                row["role"],
                self._format_summary_value(row["intensity_sum"]),
                self._format_ratio_value(row["ratio"]),
            )
            for row in self._compute_isotope_ratio_rows()
        ]
        self._set_sheet_data(self.ratio_sheet, rows)

    def _refresh_isotope_ratio_panel(self):
        self._sync_ratio_groups_with_targets()
        self._refresh_ratio_primary_list()
        self._refresh_ratio_active_principal_options()
        self._refresh_ratio_secondary_list()

    def _on_ratio_primary_select(self, _event=None):
        if self._suppress_ratio_primary_select:
            return
        selected_indices = self.ratio_primary_lb.curselection()
        self.ratio_principals = {self.targets[index] for index in selected_indices}
        self._sync_ratio_groups_with_targets()
        self._refresh_ratio_active_principal_options()
        self._refresh_ratio_secondary_list()

    def _on_ratio_secondary_select(self, _event=None):
        if self._suppress_ratio_secondary_select:
            return
        active_principal = self._get_active_ratio_principal()
        if active_principal is None:
            return
        selected_indices = self.ratio_secondary_lb.curselection()
        self.ratio_secondary_map[active_principal] = [
            self._ratio_secondary_targets_view[index]
            for index in selected_indices
        ]
        self._refresh_isotope_ratio_table()
        self._refresh_ratio_status()

    def _refresh_summary_table(self):
        start_time = time.perf_counter()
        LOGGER.info("Refresh summary table begin: targets=%d", len(self.targets))
        if not hasattr(self, "summary_sheet"):
            self._refresh_isotope_ratio_panel()
            return
        rows = []
        for target in self.targets:
            metrics = self._get_target_metrics(target)
            values = [
                self._get_assigned_peak_label(target),
                self._format_target_value(target),
            ]
            for _column_id, _label, prefix in SUMMARY_METHODS:
                values.append(self._format_summary_value(metrics.get(f"sum_intensity_{prefix}")))
            rows.append(values)
        self._set_sheet_data(self.summary_sheet, rows)
        if self.targets:
            self._set_sheet_current_row(self.summary_sheet, self.current_idx)
        self._refresh_isotope_ratio_panel()
        LOGGER.info("Refresh summary table end in %s", _format_duration_ms(start_time))

    def _on_summary_select(self, _event):
        if not hasattr(self, "summary_sheet"):
            return
        idx = self._get_sheet_current_row(self.summary_sheet)
        if idx is None:
            return
        if idx != self.current_idx and 0 <= idx < len(self.targets):
            self._navigate_to(idx)

    def _on_table_select(self, _event):
        if self._suppress_table_select:
            return
        selected = self.data_sheet.get_currently_selected()
        if not selected or selected.row is None:
            return
        try:
            selected_x = float(self.analysis["x"][selected.row])
            if self.selected_x is not None and abs(float(self.selected_x) - selected_x) < 1e-12:
                return
            self._set_selected_x(selected_x, sync_table=False)
        except (IndexError, TypeError, ValueError):
            pass

    def _select_table_row_for_x(self, x_val, focus_table=False):
        if not hasattr(self, "data_sheet") or self.analysis is None:
            return

        row_count = self.data_sheet.get_total_rows()
        if not row_count:
            return

        idx = int(np.argmin(np.abs(self.analysis["x"] - float(x_val))))
        if idx < 0 or idx >= row_count:
            return

        current_selection = self.data_sheet.get_currently_selected()
        if current_selection and current_selection.row == idx and current_selection.column == 0:
            if focus_table:
                self.data_sheet.see(idx, 0)
                self.data_sheet.focus_set()
            return

        self._suppress_table_select = True
        try:
            self.data_sheet.select_cell(idx, 0, redraw=False, run_binding_func=False)
            self.data_sheet.see(idx, 0, redraw=True)
            if focus_table:
                self.data_sheet.focus_set()
        finally:
            self._suppress_table_select = False

    def _select_auto_peak(self, refresh_plot=True, focus_table=False):
        auto_peak = self.analysis.get("auto_peak") if self.analysis else None
        if auto_peak is None:
            self.selected_x = None
            self.sel_label_var.set("Selected x: —")
            if refresh_plot:
                self._refresh_plot()
            return

        self.selected_x = float(auto_peak)
        self.sel_label_var.set(f"Selected x:\n{self.selected_x:.6f}")
        self._select_table_row_for_x(self.selected_x, focus_table=focus_table)
        if refresh_plot:
            self._refresh_plot()

    # ── Plot ──────────────────────────────────────────────────────────────────

    def _refresh_plot(self):
        start_time = time.perf_counter()
        LOGGER.info("Refresh plot begin")
        self.ax1.cla()
        self.ax2.cla()
        self.ax3.cla()
        if not self.analysis:
            self.canvas.draw_idle()
            LOGGER.info("Refresh plot end: no analysis in %s", _format_duration_ms(start_time))
            return

        a   = self.analysis
        x   = a["x"]
        man = self.manual.get(self.targets[self.current_idx], {})
        use_smoothing = bool(self.use_sg_smoothing_var.get())
        mode = str(self.mode_var.get())

        # raw is always visible
        self.ax1.plot(x, a["y_raw"], color="#5B9BD5", alpha=0.35,
                      linewidth=1.0, label="Raw")
        if use_smoothing:
            self.ax1.plot(x, a["y_sm"],  color="#2E75B6",
                          linewidth=1.8, label="Smoothed")
        self.ax1.set_xlabel("m/z", fontsize=8)
        self.ax1.set_ylabel("Intensity", fontsize=8, color="#2E75B6")
        self.ax1.tick_params(axis="y", labelcolor="#2E75B6", labelsize=7)
        self.ax1.tick_params(axis="x", labelsize=7)
        self.ax3.spines["right"].set_position(("axes", 1.12))
        self.ax3.set_frame_on(True)
        self.ax3.patch.set_visible(False)

        # derivatives
        if self.show_deriv_var.get():
            self.ax2.plot(x, a["dy"],  color="#ED7D31", linestyle="--",
                          alpha=0.75, linewidth=1.0, label="dy/dx")
            self.ax2.axhline(0, color="#BBBBBB", linewidth=0.7, zorder=0)
            self.ax2.set_ylabel("dy/dx", fontsize=8, color="#ED7D31")
            self.ax2.tick_params(axis="y", labelsize=7, labelcolor="#ED7D31", colors="#ED7D31")

            self.ax3.plot(x, a["d2y"], color="#C00000", linestyle=":",
                          alpha=0.75, linewidth=1.0, label="d²y/dx²")
            self.ax3.axhline(0, color="#DDDDDD", linewidth=0.7, zorder=0)
            self.ax3.set_ylabel("d²y/dx²", fontsize=8, color="#C00000")
            self.ax3.tick_params(axis="y", labelsize=7, labelcolor="#C00000", colors="#C00000")
        else:
            self.ax2.set_ylabel("")
            self.ax3.set_ylabel("")
            self.ax2.set_yticks([])
            self.ax3.set_yticks([])

        def vline(ax, x_val, color, ls, lw, label):
            ax.axvline(x_val, color=color, linestyle=ls, linewidth=lw, label=label)

        vline(self.ax1, a["auto_peak"],  "#1A7A1A", "-",  1.8, "Auto peak")
        vline(self.ax1, a["auto_left"],  "#1A7A1A", "--", 1.4, f"Auto left ({mode})")
        vline(self.ax1, a["auto_right"], "#1A7A1A", "--", 1.4, f"Auto right ({mode})")
        vline(self.ax1, a["grad_left"], "#00A3A3", ":", 1.0, "Gradient left")
        vline(self.ax1, a["grad_right"], "#00A3A3", ":", 1.0, "Gradient right")

        # manual points
        m_color = {"peak": "#7030A0", "left": "#9B59B6", "right": "#9B59B6"}
        m_style = {"peak": "-",       "left": "--",       "right": "--"}
        m_lw    = {"peak": 2.0,       "left": 1.3,        "right": 1.3}
        for role in ("left", "peak", "right"):
            if man.get(role) is not None:
                vline(self.ax1, man[role],
                      m_color[role], m_style[role], m_lw[role], f"Man {role}")

        # selected x
        if self.selected_x is not None:
            self.ax1.axvline(self.selected_x, color="red",
                             linewidth=1.4, alpha=0.5, label="Selected")

        h1, l1 = self.ax1.get_legend_handles_labels()
        h2, l2 = self.ax2.get_legend_handles_labels()
        h3, l3 = self.ax3.get_legend_handles_labels()
        self.ax1.legend(h1 + h2 + h3, l1 + l2 + l3, fontsize=6.5,
                        loc="upper right", framealpha=0.7)

        # display window around selected point (or auto peak if none selected)
        view_center = self.selected_x if self.selected_x is not None else a["auto_peak"]
        x_min = float(np.min(x))
        x_max = float(np.max(x))
        try:
            view_half = abs(float(self.view_win_var.get()))
        except Exception:
            view_half = DEFAULT_VIEW_WIN
        if view_half > 0:
            x_min = max(x_min, view_center - view_half)
            x_max = min(x_max, view_center + view_half)
        if x_max > x_min:
            self.ax1.set_xlim(x_min, x_max)

            visible_mask = (x >= x_min) & (x <= x_max)
            visible_series = [a["y_raw"]]
            if use_smoothing:
                visible_series.append(a["y_sm"])
            visible_values = np.concatenate([
                np.asarray(series)[visible_mask]
                for series in visible_series
            ])
            visible_values = visible_values[np.isfinite(visible_values)]
            if visible_values.size:
                y_visible_min = float(np.min(visible_values))
                y_visible_max = float(np.max(visible_values))
                if y_visible_max <= y_visible_min:
                    pad = max(abs(y_visible_max) * 0.05, 1.0)
                else:
                    pad = max((y_visible_max - y_visible_min) * 0.08, abs(y_visible_max) * 0.02, 1e-12)

                if y_visible_min >= 0:
                    y_lower = 0.0
                else:
                    y_lower = y_visible_min - pad
                self.ax1.set_ylim(y_lower, y_visible_max + pad)

        self.ax1.set_title(
            f"Target {self.current_idx + 1}/{len(self.targets)}: "
            f"{self._format_target_tag(self.targets[self.current_idx])}", fontsize=9)

        if self.label_matches:
            best = self.label_matches[0]
            self.ax1.text(
                0.01,
                0.98,
                f"Best label: {best['label']} | error {best['error_mmu']:+.5f} mmu",
                transform=self.ax1.transAxes,
                va="top",
                ha="left",
                fontsize=8,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="#fff2cc", alpha=0.85),
            )
        self.canvas.draw_idle()
        LOGGER.info("Refresh plot end in %s", _format_duration_ms(start_time))

    def _on_plot_click(self, event):
        if self.canvas.toolbar.mode:   # pan/zoom active — ignore clicks
            return
        if event.inaxes not in (self.ax1, self.ax2, self.ax3):
            return
        if event.xdata is not None:
            self._set_selected_x(event.xdata)

    # ── Manual assignment ─────────────────────────────────────────────────────

    def _set_selected_x(self, x_val, sync_table=True):
        x_val = float(x_val)
        if self.selected_x is not None and abs(float(self.selected_x) - x_val) < 1e-12:
            return

        self.selected_x = x_val
        self.sel_label_var.set(f"Selected x:\n{self.selected_x:.6f}")
        if sync_table:
            self._select_table_row_for_x(self.selected_x)
        self._refresh_plot()

    def _assign(self, role):
        if self.selected_x is None:
            self._status("No point selected. Click plot or table row first.")
            return
        if not self.targets:
            return
        t = self.targets[self.current_idx]
        self.manual.setdefault(t, {})[role] = self.selected_x
        self._status(f"Manual {role} → {self.selected_x:.6f}")
        # Clear selection after assignment so next role requires a fresh click.
        self.selected_x = None
        self.sel_label_var.set("Selected x: —")
        self._refresh_summary_table()
        self._refresh_plot()

    def _clear_manual(self):
        if not self.targets:
            return
        self.manual.pop(self.targets[self.current_idx], None)
        self._status("Manual points cleared.")
        self._refresh_summary_table()
        self._refresh_plot()

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _status(self, msg):
        LOGGER.info("STATUS: %s", msg)
        self.status_var.set(msg)

    def report_callback_exception(self, exc, val, tb):
        _log_exception("Tk callback exception", (exc, val, tb))
        try:
            messagebox.showerror(
                "Unexpected Error",
                f"An unexpected error occurred.\n\nSee log file:\n{LOG_FILE_PATH}",
            )
        except Exception:
            pass


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    LOGGER.info("Application startup")
    app = PeakInspector()
    app.mainloop()
