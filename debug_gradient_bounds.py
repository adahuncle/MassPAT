import numpy as np

import auto_selection as sel


# Hardcoded fixed dataset (20 rows)
X = np.array([
    165.9104230,
    165.9109062,
    165.9113894,
    165.9118727,
    165.9123559,
    165.9128391,
    165.9133223,
    165.9138056,
    165.9142888,
    165.9147720,
    165.9152553,
    165.9157385,
    165.9162217,
    165.9167050,
    165.9171882,
    165.9176715,
    165.9181547,
    165.9186380,
    165.9191212,
    165.9196045,
], dtype=float)

Y = np.array([
    -0.291569535,
    0.119162865,
    1.046596322,
    0.612238989,
    -0.243709642,
    0.945130113,
    2.245365940,
    0.646761579,
    6.434522114,
    44.857117430,
    128.376757800,
    226.440072000,
    278.024188600,
    245.636306900,
    152.681839700,
    60.109361000,
    9.951115871,
    -1.138804480,
    0.906038766,
    0.650151950,
], dtype=float)

DY = np.array([
    850.025662267836,
    1384.691491127352,
    510.458285593214,
    -1335.119743888753,
    344.465153017666,
    2575.616289301252,
    -308.089681284275,
    4335.938686107076,
    45747.470873821308,
    126159.240829842282,
    187879.946210629161,
    154850.404390364827,
    19879.577680137183,
    -129699.829037397518,
    -191957.608812402817,
    -147668.856710128253,
    -63379.463308906197,
    -9355.777851638015,
    1851.456351099823,
    -529.457512948506,
], dtype=float)


def fmt_row(idx):
    return f"row={idx + 1:02d} idx={idx:02d}"


def find_index_for_x(x_val):
    return int(np.argmin(np.abs(X - float(x_val))))


def compute_auto_eps(peak_idx):
    n = len(DY)
    abs_dy = np.abs(DY)
    max_abs_dy = float(np.max(abs_dy)) if n else 0.0

    flank_radius = max(
        sel.DEFAULT_GRADIENT_MIN_DISTANCE_POINTS,
        max(1, int(n * 0.05)),
    )
    flank_lo = max(0, peak_idx - flank_radius)
    flank_hi = min(n, peak_idx + flank_radius + 1)
    peak_region_scale = (
        float(np.max(abs_dy[flank_lo:flank_hi])) if flank_hi > flank_lo else max_abs_dy
    )
    return sel.DEFAULT_GRADIENT_EPS_FRACTION * max(peak_region_scale, max_abs_dy * 0.001)


def format_inclusive_rows(left_idx, right_idx):
    return f"rows {left_idx + 1}-{right_idx + 1} inclusive"


def format_exclusive_rows(left_idx, right_idx):
    left_exclusive = left_idx
    right_exclusive = right_idx + 2
    return f"rows {left_exclusive}-{right_exclusive} exclusive"


def matches_intended_region(left_idx, right_idx):
    return left_idx == 7 and right_idx == 16


def run_case(name, peak_idx, eps, intensity_fraction):
    left_x, right_x, debug = sel.get_gradient_bounds(
        X,
        Y,
        DY,
        peak_idx,
        eps=eps,
        intensity_fraction=intensity_fraction,
        return_debug=True,
    )

    left_idx = find_index_for_x(left_x)
    right_idx = find_index_for_x(right_x)

    eps_label = "auto" if eps is None else f"{float(eps):.6f}"

    print(f"\n=== {name} ===")
    print(f"eps used: {eps_label}")
    print(f"intensity_fraction: {float(intensity_fraction):.6f}")
    print(
        "left transition:",
        f"idx={debug['left']['transition_idx']} (row {debug['left']['transition_idx'] + 1}),",
        f"reason={debug['left']['transition_reason']}",
    )
    print(
        "right transition:",
        f"idx={debug['right']['transition_idx']} (row {debug['right']['transition_idx'] + 1}),",
        f"reason={debug['right']['transition_reason']}",
    )
    print(
        "left stop:",
        f"idx={debug['left']['stop_idx']} (row {debug['left']['stop_idx'] + 1}),",
        f"reason={debug['left']['stop_reason']}",
    )
    print(
        "right stop:",
        f"idx={debug['right']['stop_idx']} (row {debug['right']['stop_idx'] + 1}),",
        f"reason={debug['right']['stop_reason']}",
    )
    print(
        "final left bound:",
        f"idx={left_idx} (row {left_idx + 1}), x={left_x:.7f}",
    )
    print(
        "final right bound:",
        f"idx={right_idx} (row {right_idx + 1}), x={right_x:.7f}",
    )
    print(f"kept rows inclusive: {format_inclusive_rows(left_idx, right_idx)}")
    print(f"equivalent exclusive boundaries: {format_exclusive_rows(left_idx, right_idx)}")
    print(
        "matches intended rows 8-17 inclusive:",
        "YES" if matches_intended_region(left_idx, right_idx) else "NO",
    )


def main():
    # Known peak around 165.9162 region
    peak_idx = 12
    auto_eps = compute_auto_eps(peak_idx)
    intensity_fraction = sel.DEFAULT_INTENSITY_FRACTION

    print("=== Peak selection ===")
    print(f"peak_idx={peak_idx} (row {peak_idx + 1}), x={X[peak_idx]:.7f}, y={Y[peak_idx]:.9f}")
    print(f"auto eps estimate={auto_eps:.6f}")
    print(f"intensity_fraction={float(intensity_fraction):.6f}")
    print("intended kept region = rows 8-17 inclusive")
    print("equivalent expected boundaries = rows 7-18 exclusive")

    run_case("1) eps = auto", peak_idx, eps=None, intensity_fraction=intensity_fraction)
    run_case("2) eps = auto * 0.5", peak_idx, eps=auto_eps * 0.5, intensity_fraction=intensity_fraction)
    run_case("3) eps = auto * 0.25", peak_idx, eps=auto_eps * 0.25, intensity_fraction=intensity_fraction)
    run_case("4) eps = auto * 0.1", peak_idx, eps=auto_eps * 0.1, intensity_fraction=intensity_fraction)
    run_case("5) eps = auto * 0.05", peak_idx, eps=auto_eps * 0.05, intensity_fraction=intensity_fraction)


if __name__ == "__main__":
    main()
