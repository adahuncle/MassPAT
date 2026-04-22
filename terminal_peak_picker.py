#!/usr/bin/env python3
"""
Minimal terminal peak picker.

Runs automatic peak selection for one or more target m/z values against one or
more spectrum files using a JSON config. Outputs:
1) A summary CSV with selected peak bounds and summed intensity per peak.
2) A per-peak CSV with points inside/outside selected peak bounds.
3) A plot image around each target peak window.

Usage:
    python terminal_peak_picker.py --config peak_picker_config.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from auto_selection import analyze_window


APP_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG_PATH = APP_DIR / "peak_picker_config.json"


def load_spectrum(path: Path) -> pd.DataFrame:
    """Load a two-column CSV (m/z, intensity), with or without header."""
    for skip in (0, 1):
        try:
            df = pd.read_csv(
                path,
                skiprows=skip,
                header=None,
                names=["x", "y"],
                usecols=[0, 1],
            )
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
        f"Could not parse CSV at {path}. Expected two numeric columns (m/z, intensity)."
    )


def load_config(config_path: Path) -> dict[str, Any]:
    with config_path.open("r", encoding="utf-8") as fh:
        config = json.load(fh)

    if not isinstance(config, dict):
        raise ValueError("Config root must be a JSON object.")

    targets = config.get("targets")
    if not isinstance(targets, list) or not targets:
        raise ValueError("Config field 'targets' must be a non-empty list of m/z values.")

    if "data_files" in config:
        data_files = config.get("data_files")
    else:
        single = config.get("data_file")
        data_files = [single] if single else []

    if not isinstance(data_files, list) or not data_files:
        raise ValueError(
            "Config must include 'data_files' (list) or 'data_file' (single path)."
        )

    config["targets"] = [float(v) for v in targets]
    config["data_files"] = [str(v) for v in data_files]
    config.setdefault("output_dir", "output/peak_picker")
    config.setdefault("analysis", {})

    if not isinstance(config["analysis"], dict):
        raise ValueError("Config field 'analysis' must be an object.")

    return config


def plot_peak_window(
    sub: pd.DataFrame,
    target: float,
    auto_peak: float,
    left: float,
    right: float,
    out_png: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.8), dpi=130)

    ax.plot(sub["x"], sub["y"], label="Raw intensity", linewidth=1.1)
    ax.plot(sub["x"], sub["smoothed"], label="Smoothed", linewidth=1.0, alpha=0.9)

    lo = min(left, right)
    hi = max(left, right)
    in_peak = sub[(sub["x"] >= lo) & (sub["x"] <= hi)]
    if not in_peak.empty:
        ax.fill_between(in_peak["x"], in_peak["y"], alpha=0.25, label="Selected peak area")

    ax.axvline(target, linestyle="--", linewidth=1.0, label=f"Target {target:.6f}")
    ax.axvline(auto_peak, linestyle="-", linewidth=1.1, label=f"Auto peak {auto_peak:.6f}")
    ax.axvline(left, linestyle=":", linewidth=1.0, label=f"Left {left:.6f}")
    ax.axvline(right, linestyle=":", linewidth=1.0, label=f"Right {right:.6f}")

    ax.set_title(f"Peak window around target m/z {target:.6f}")
    ax.set_xlabel("m/z")
    ax.set_ylabel("Intensity")
    ax.legend(loc="best", fontsize=8)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def run(config_path: Path) -> None:
    config = load_config(config_path)
    output_dir = (APP_DIR / config["output_dir"]).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    analysis_cfg = config.get("analysis", {})
    allowed_analysis_keys = {
        "half_window",
        "sg_window",
        "sg_poly",
        "use_smoothing",
        "intensity_fraction",
        "target_tolerance",
        "min_height_fraction",
        "use_prominence",
        "gradient_min_distance_points",
        "gradient_baseline_run_length",
        "gradient_positive_only",
    }
    analyze_kwargs = {
        k: analysis_cfg[k] for k in allowed_analysis_keys if k in analysis_cfg
    }

    summary_rows: list[dict[str, Any]] = []

    for file_entry in config["data_files"]:
        data_path = Path(file_entry)
        if not data_path.is_absolute():
            data_path = (APP_DIR / data_path).resolve()

        if not data_path.exists():
            print(f"[skip] Missing file: {data_path}")
            continue

        print(f"[load] {data_path}")
        df = load_spectrum(data_path)

        for target in config["targets"]:
            result = analyze_window(df, float(target), **analyze_kwargs)
            if result is None:
                print(f"  [no-data] target={target:.6f}")
                summary_rows.append(
                    {
                        "file": str(data_path),
                        "target_mz": float(target),
                        "status": "insufficient_window_data",
                    }
                )
                continue

            sub = result["sub"].copy()
            left = float(result["auto_left"])
            right = float(result["auto_right"])
            peak_mz = float(result["auto_peak"])
            lo = min(left, right)
            hi = max(left, right)

            sub["in_selected_peak"] = (sub["x"] >= lo) & (sub["x"] <= hi)
            selected = sub[sub["in_selected_peak"]]
            sum_intensity = float(selected["y"].sum()) if not selected.empty else 0.0

            base = f"{data_path.stem}_target_{target:.6f}".replace(" ", "_")
            points_csv = output_dir / f"{base}_points.csv"
            plot_png = output_dir / f"{base}.png"

            sub.to_csv(points_csv, index=False)
            plot_peak_window(sub, float(target), peak_mz, left, right, plot_png)

            summary_rows.append(
                {
                    "file": str(data_path),
                    "target_mz": float(target),
                    "status": "ok",
                    "auto_peak_mz": peak_mz,
                    "left_bound_mz": left,
                    "right_bound_mz": right,
                    "points_in_selected_peak": int(len(selected)),
                    "sum_intensity_selected_peak": sum_intensity,
                    "points_csv": str(points_csv),
                    "plot_png": str(plot_png),
                }
            )

            print(
                "  [ok] "
                f"target={target:.6f} peak={peak_mz:.6f} "
                f"sum={sum_intensity:.4f} points={len(selected)}"
            )

    summary_path = output_dir / "peak_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    print(f"[done] Wrote summary: {summary_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Minimal terminal peak picker")
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help="Path to config JSON file.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.config)
