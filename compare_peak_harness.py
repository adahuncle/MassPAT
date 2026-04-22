"""Minimal harness to compare analyze_window_simple vs analyze_window.

Input:
- df: pandas DataFrame with columns x and y
- targets: list of target m/z values

Output:
- pandas DataFrame with method-specific values and absolute differences
"""

from __future__ import annotations

import argparse
from typing import Any, Iterable
from pathlib import Path
import tkinter as tk
from tkinter import ttk

import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from auto_selection import analyze_window
from terminal_peak_picker import load_spectrum
from streamline_peak_picker import (
    DEFAULT_HALF_WIN,
    DEFAULT_TARGET_TOLERANCE,
    analyze_window_simple,
)

_FIELDS = ("auto_peak", "auto_left", "auto_right", "auto_peak_int")


def _validate_input_df(df: pd.DataFrame) -> None:
    missing = {"x", "y"} - set(df.columns)
    if missing:
        raise ValueError(f"Input DataFrame is missing required columns: {sorted(missing)}")


def _extract_fields(result: dict[str, Any] | None) -> dict[str, float | None]:
    if result is None:
        return {field: None for field in _FIELDS}
    return {field: result.get(field) for field in _FIELDS}


def _abs_diff(a: float | None, b: float | None) -> float | None:
    if a is None or b is None:
        return None
    return abs(float(a) - float(b))


def _percent_diff(reference: float | None, value: float | None) -> float | None:
    """Percent difference relative to reference: ((value-reference)/reference)*100."""
    if reference is None or value is None:
        return None
    ref = float(reference)
    if ref == 0.0:
        return None
    return ((float(value) - ref) / ref) * 100.0


def _sum_selected_intensity(result: dict[str, Any] | None) -> float | None:
    """Sum raw intensities for points inside the selected [left, right] bounds."""
    if result is None:
        return None

    sub = result.get("sub")
    left = result.get("auto_left")
    right = result.get("auto_right")
    if sub is None or left is None or right is None:
        return None

    lo = min(float(left), float(right))
    hi = max(float(left), float(right))
    in_range = sub[(sub["x"] >= lo) & (sub["x"] <= hi)]
    if in_range.empty:
        return 0.0
    return float(in_range["y"].sum())


def compare_analyze_windows(
    df: pd.DataFrame,
    targets: Iterable[float],
    half_window: float = DEFAULT_HALF_WIN,
    target_tolerance: float = DEFAULT_TARGET_TOLERANCE,
    analyze_window_kwargs: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """Run both analyzers on each target and return a comparison DataFrame.

    Notes:
    - Both methods receive the same shared parameters: half_window, target_tolerance.
    - analyze_window is always run with use_smoothing=False.
    """
    _validate_input_df(df)

    base_kwargs = {
        "half_window": half_window,
        "target_tolerance": target_tolerance,
        "use_smoothing": False,
    }
    if analyze_window_kwargs:
        base_kwargs.update(analyze_window_kwargs)
    base_kwargs["use_smoothing"] = False

    rows: list[dict[str, Any]] = []
    for target in targets:
        target_value = float(target)
        simple_result = analyze_window_simple(
            df,
            target_value,
            half_window=half_window,
            target_tolerance=target_tolerance,
        )
        full_result = analyze_window(df, target_value, **base_kwargs)

        simple_sum_intensity = _sum_selected_intensity(simple_result)
        full_sum_intensity = _sum_selected_intensity(full_result)
        sum_delta = None
        if simple_sum_intensity is not None and full_sum_intensity is not None:
            sum_delta = float(full_sum_intensity) - float(simple_sum_intensity)

        simple_vals = _extract_fields(simple_result)
        full_vals = _extract_fields(full_result)

        row: dict[str, Any] = {
            "target": target_value,
            "simple_has_result": simple_result is not None,
            "full_has_result": full_result is not None,
        }

        for field in _FIELDS:
            s_val = simple_vals[field]
            f_val = full_vals[field]
            row[f"simple_{field}"] = s_val
            row[f"full_{field}"] = f_val
            row[f"abs_diff_{field}"] = _abs_diff(s_val, f_val)

        row["simple_sum_intensity"] = simple_sum_intensity
        row["full_sum_intensity"] = full_sum_intensity
        row["delta_sum_intensity"] = sum_delta
        row["abs_diff_sum_intensity"] = _abs_diff(simple_sum_intensity, full_sum_intensity)
        row["percent_diff_sum_intensity"] = _percent_diff(simple_sum_intensity, full_sum_intensity)

        rows.append(row)

    return pd.DataFrame(rows)


def filter_mismatches(
    comparison_df: pd.DataFrame,
    tolerance: float = 1e-6,
    include_missing: bool = True,
) -> pd.DataFrame:
    """Return rows where methods DISAGREE, not where both are False."""
    
    # First, rows where both have results but values differ
    diff_cols = [f"abs_diff_{field}" for field in _FIELDS] + ["abs_diff_sum_intensity"]
    value_mismatch = comparison_df[diff_cols].gt(float(tolerance)).any(axis=1)
    
    # Second, rows where one has a result and the other doesn't (disagreement)
    has_result_mismatch = comparison_df["simple_has_result"] != comparison_df["full_has_result"]
    
    # Only count as mismatch if they DISAGREE
    return comparison_df[value_mismatch | has_result_mismatch].copy()


def inspect_target(
    df: pd.DataFrame,
    target: float,
    half_window: float = DEFAULT_HALF_WIN,
    target_tolerance: float = DEFAULT_TARGET_TOLERANCE,
    analyze_window_kwargs: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """Print side-by-side values for a single target and return them as a DataFrame."""
    _validate_input_df(df)

    base_kwargs = {
        "half_window": half_window,
        "target_tolerance": target_tolerance,
        "use_smoothing": False,
    }
    if analyze_window_kwargs:
        base_kwargs.update(analyze_window_kwargs)
    base_kwargs["use_smoothing"] = False

    target_value = float(target)
    simple_result = analyze_window_simple(
        df,
        target_value,
        half_window=half_window,
        target_tolerance=target_tolerance,
    )
    full_result = analyze_window(df, target_value, **base_kwargs)

    simple_sum_intensity = _sum_selected_intensity(simple_result)
    full_sum_intensity = _sum_selected_intensity(full_result)
    sum_delta = None
    if simple_sum_intensity is not None and full_sum_intensity is not None:
        sum_delta = float(full_sum_intensity) - float(simple_sum_intensity)

    simple_vals = _extract_fields(simple_result)
    full_vals = _extract_fields(full_result)

    inspect_rows = []
    for field in _FIELDS:
        s_val = simple_vals[field]
        f_val = full_vals[field]
        inspect_rows.append(
            {
                "metric": field,
                "simple": s_val,
                "full": f_val,
                "abs_diff": _abs_diff(s_val, f_val),
            }
        )

    inspect_rows.append(
        {
            "metric": "sum_intensity_selected_peak",
            "simple": simple_sum_intensity,
            "full": full_sum_intensity,
            "abs_diff": _abs_diff(simple_sum_intensity, full_sum_intensity),
        }
    )
    inspect_rows.append(
        {
            "metric": "delta_sum_intensity(full-simple)",
            "simple": None,
            "full": sum_delta,
            "abs_diff": abs(float(sum_delta)) if sum_delta is not None else None,
        }
    )
    inspect_rows.append(
        {
            "metric": "percent_diff_sum_intensity(%)",
            "simple": None,
            "full": _percent_diff(simple_sum_intensity, full_sum_intensity),
            "abs_diff": None,
        }
    )

    out = pd.DataFrame(inspect_rows)

    print(f"Target: {target_value}")
    print(f"simple_has_result={simple_result is not None}, full_has_result={full_result is not None}")
    print(out.to_string(index=False))

    return out


def _parse_targets(raw_targets: str) -> list[float]:
    return [float(part.strip()) for part in raw_targets.split(",") if part.strip()]


def _build_base_kwargs(
    half_window: float,
    target_tolerance: float,
    analyze_window_kwargs: dict[str, Any] | None = None,
) -> dict[str, Any]:
    base_kwargs = {
        "half_window": half_window,
        "target_tolerance": target_tolerance,
        "use_smoothing": False,
    }
    if analyze_window_kwargs:
        base_kwargs.update(analyze_window_kwargs)
    base_kwargs["use_smoothing"] = False
    return base_kwargs


def _run_both_methods(
    df: pd.DataFrame,
    target: float,
    half_window: float,
    target_tolerance: float,
) -> dict[str, Any]:
    base_kwargs = _build_base_kwargs(half_window, target_tolerance)
    simple_result = analyze_window_simple(
        df,
        target,
        half_window=half_window,
        target_tolerance=target_tolerance,
    )
    full_result = analyze_window(df, target, **base_kwargs)

    simple_sum_intensity = _sum_selected_intensity(simple_result)
    full_sum_intensity = _sum_selected_intensity(full_result)
    sum_delta = None
    if simple_sum_intensity is not None and full_sum_intensity is not None:
        sum_delta = float(full_sum_intensity) - float(simple_sum_intensity)
    percent_area_diff = _percent_diff(simple_sum_intensity, full_sum_intensity)

    simple_vals = _extract_fields(simple_result)
    full_vals = _extract_fields(full_result)
    diffs = {field: _abs_diff(simple_vals[field], full_vals[field]) for field in _FIELDS}
    diffs["sum_intensity"] = _abs_diff(simple_sum_intensity, full_sum_intensity)

    return {
        "simple": simple_result,
        "full": full_result,
        "simple_vals": simple_vals,
        "full_vals": full_vals,
        "simple_sum_intensity": simple_sum_intensity,
        "full_sum_intensity": full_sum_intensity,
        "sum_delta": sum_delta,
        "percent_area_diff": percent_area_diff,
        "diffs": diffs,
    }


def _build_cases(
    csv_files: list[Path],
    targets: list[float],
    half_window: float,
    target_tolerance: float,
) -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for csv_path in csv_files:
        df = load_spectrum(csv_path)
        for target in targets:
            run = _run_both_methods(df, target, half_window, target_tolerance)
            cases.append(
                {
                    "file": csv_path,
                    "target": target,
                    "simple": run["simple"],
                    "full": run["full"],
                    "simple_vals": run["simple_vals"],
                    "full_vals": run["full_vals"],
                    "simple_sum_intensity": run["simple_sum_intensity"],
                    "full_sum_intensity": run["full_sum_intensity"],
                    "sum_delta": run["sum_delta"],
                    "percent_area_diff": run["percent_area_diff"],
                    "diffs": run["diffs"],
                }
            )
    return cases


def _draw_method(ax, result: dict[str, Any] | None, target: float, title: str) -> None:
    ax.clear()
    ax.set_title(title)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Intensity")

    if result is None:
        ax.text(0.5, 0.5, "No result", ha="center", va="center", transform=ax.transAxes)
        ax.grid(alpha=0.2)
        return

    sub = result["sub"]
    x = sub["x"].to_numpy()
    y = sub["y"].to_numpy()
    y_sm = result.get("y_sm", y)

    ax.plot(x, y, label="raw", linewidth=1.0)
    if len(y_sm) == len(x):
        ax.plot(x, y_sm, label="method", linewidth=1.0, alpha=0.85)

    peak = result.get("auto_peak")
    left = result.get("auto_left")
    right = result.get("auto_right")

    ax.axvline(target, linestyle="--", linewidth=1.0, label="target")
    if peak is not None:
        ax.axvline(float(peak), linestyle="-", linewidth=1.1, label="peak")
    if left is not None:
        ax.axvline(float(left), linestyle=":", linewidth=1.0, label="left")
    if right is not None:
        ax.axvline(float(right), linestyle=":", linewidth=1.0, label="right")

    # Show a domain that is 3x the selected peak width:
    # one peak-width to the left and right of [left, right].
    if left is not None and right is not None:
        left_f = float(left)
        right_f = float(right)
        if right_f > left_f:
            width = right_f - left_f
            ax.set_xlim(left_f - width, right_f + width)

    ax.grid(alpha=0.2)
    ax.legend(fontsize=8, loc="best")


def _fmt_num(value: float | None) -> str:
    if value is None:
        return "None"
    return f"{float(value):.6f}"


def launch_gui(cases: list[dict[str, Any]]) -> None:
    root = tk.Tk()
    root.title("Peak Comparison Viewer")
    root.geometry("1200x700")

    info_var = tk.StringVar(value="")
    diff_var = tk.StringVar(value="")

    top = ttk.Frame(root, padding=8)
    top.pack(fill="x")
    ttk.Label(top, textvariable=info_var, justify="left").pack(anchor="w")
    ttk.Label(top, textvariable=diff_var, justify="left").pack(anchor="w", pady=(4, 0))

    fig = Figure(figsize=(12, 5), dpi=100)
    ax_left = fig.add_subplot(1, 2, 1)
    ax_right = fig.add_subplot(1, 2, 2)
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack(fill="both", expand=True)

    nav = ttk.Frame(root, padding=8)
    nav.pack(fill="x")

    state = {"index": 0}

    def render_current() -> None:
        case = cases[state["index"]]
        file_name = case["file"].name
        target = case["target"]
        i = state["index"] + 1
        n = len(cases)

        info_var.set(f"[{i}/{n}] File: {file_name} | Target: {target:.6f}")
        diff_var.set(
            " | ".join(
                [
                    f"peak diff={_fmt_num(case['diffs']['auto_peak'])}",
                    f"left diff={_fmt_num(case['diffs']['auto_left'])}",
                    f"right diff={_fmt_num(case['diffs']['auto_right'])}",
                    f"int diff={_fmt_num(case['diffs']['auto_peak_int'])}",
                    f"simple area={_fmt_num(case['simple_sum_intensity'])}",
                    f"full area={_fmt_num(case['full_sum_intensity'])}",
                    f"delta area={_fmt_num(case['sum_delta'])}",
                    f"abs area diff={_fmt_num(case['diffs']['sum_intensity'])}",
                    f"pct area diff={_fmt_num(case['percent_area_diff'])}%",
                ]
            )
        )

        _draw_method(ax_left, case["simple"], target, "analyze_window_simple")
        _draw_method(ax_right, case["full"], target, "analyze_window (no smoothing)")
        fig.tight_layout()
        canvas.draw_idle()

    def prev_case(*_args: Any) -> None:
        state["index"] = (state["index"] - 1) % len(cases)
        render_current()

    def next_case(*_args: Any) -> None:
        state["index"] = (state["index"] + 1) % len(cases)
        render_current()

    ttk.Button(nav, text="Previous", command=prev_case).pack(side="left")
    ttk.Button(nav, text="Next", command=next_case).pack(side="left", padx=(8, 0))
    ttk.Label(nav, text="Use Left/Right arrow keys to cycle").pack(side="left", padx=(16, 0))

    root.bind("<Left>", prev_case)
    root.bind("<Right>", next_case)

    render_current()
    root.mainloop()


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compare analyze_window_simple vs analyze_window")
    parser.add_argument(
        "--samples-dir",
        default="data/samples",
        help="Directory containing CSV files to analyze.",
    )
    parser.add_argument(
        "--targets",
        default="157.901, 165.916",
        help="Comma-separated target m/z values.",
    )
    parser.add_argument(
        "--half-window",
        type=float,
        default=DEFAULT_HALF_WIN,
        help="Half-width of analysis window around each target.",
    )
    parser.add_argument(
        "--target-tolerance",
        type=float,
        default=DEFAULT_TARGET_TOLERANCE,
        help="Target tolerance used by both methods.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-6,
        help="Absolute difference threshold for mismatch filtering.",
    )
    parser.add_argument(
        "--inspect-target",
        type=float,
        default=None,
        help="If provided, also print side-by-side details for one target.",
    )
    parser.add_argument(
        "--gui",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Open simple side-by-side viewer window.",
    )
    return parser


def _find_csv_files(samples_dir: str) -> list[Path]:
    base = Path(samples_dir)
    if not base.exists() or not base.is_dir():
        raise ValueError(f"Samples directory does not exist or is not a directory: {samples_dir}")
    return sorted(base.glob("*.csv"))


def main() -> None:
    args = _build_arg_parser().parse_args()
    targets = _parse_targets(args.targets)
    if not targets:
        raise ValueError("No targets provided. Set --targets to one or more comma-separated values.")
    csv_files = _find_csv_files(args.samples_dir)
    if not csv_files:
        print(f"No CSV files found in {args.samples_dir}")
        return

    cases = _build_cases(
        csv_files,
        targets,
        half_window=args.half_window,
        target_tolerance=args.target_tolerance,
    )
    if not cases:
        print("No file/target cases to display.")
        return

    total_files = 0
    total_mismatches = 0

    for csv_path in csv_files:
        df = load_spectrum(csv_path)

        comparison = compare_analyze_windows(
            df,
            targets,
            half_window=args.half_window,
            target_tolerance=args.target_tolerance,
        )
        mismatches = filter_mismatches(comparison, tolerance=args.tolerance)

        print(f"\n=== {csv_path} ===")
        print("Comparison:")
        print(comparison.to_string(index=False))
        print()
        print(f"Mismatches (tolerance={args.tolerance}): {len(mismatches)}")
        if len(mismatches):
            print(mismatches.to_string(index=False))

        if args.inspect_target is not None:
            print()
            print("Inspect target:")
            inspect_target(
                df,
                args.inspect_target,
                half_window=args.half_window,
                target_tolerance=args.target_tolerance,
            )

        total_files += 1
        total_mismatches += len(mismatches)

    print(f"\nProcessed {total_files} file(s). Total mismatches: {total_mismatches}")

    if args.gui:
        launch_gui(cases)


if __name__ == "__main__":
    main()
