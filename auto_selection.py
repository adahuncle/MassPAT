"""
Auto peak-selection logic extracted from peak_inspector.py.
Edit this module to tweak smoothing, peak finding, and bound selection behavior.
"""

import numpy as np
from scipy.signal import find_peaks, savgol_filter


DEFAULT_HALF_WIN = 0.15
DEFAULT_SG_WINDOW = 0
DEFAULT_SG_POLY = 4
DEFAULT_INTENSITY_FRACTION = 0.01
DEFAULT_TARGET_TOLERANCE = 0.005
DEFAULT_MIN_HEIGHT_FRACTION = 0.0
DEFAULT_GRADIENT_EPS_FRACTION = 0.01
DEFAULT_GRADIENT_CONSECUTIVE_VIOLATIONS = 1
DEFAULT_GRADIENT_MIN_DISTANCE_POINTS = 3
DEFAULT_GRADIENT_BASELINE_RUN_LENGTH = 3
DEFAULT_GRADIENT_POSITIVE_ONLY = True


def _nearest_zero_crossing(arr, start_idx, direction):
    """Walk from start_idx in direction (+1 or -1); return index of first sign change."""
    n = len(arr)
    i = start_idx
    while 0 < i < n - 1:
        j = i + direction
        if 0 <= j < n and arr[i] * arr[j] < 0:
            return i if direction < 0 else j
        i += direction
    if n == 0:
        return 0
    if direction < 0:
        if start_idx <= 0:
            return 0
        return int(np.argmax(np.abs(arr[:start_idx])))
    if start_idx >= n:
        return n - 1
    return int(start_idx + np.argmax(np.abs(arr[start_idx:])))


def _is_local_minimum(y_arr, idx):
    return 0 < idx < len(y_arr) - 1 and y_arr[idx - 1] > y_arr[idx] < y_arr[idx + 1]


def _find_gradient_transition(dy, peak_idx, direction, eps, k, min_distance_points):
    n = len(dy)
    if n == 0:
        return None, "window edge"

    min_distance_points = max(1, int(min_distance_points))
    if direction < 0:
        start_idx = max(0, peak_idx - min_distance_points)
        stop = -1
        sign_reason = "negative run"
    else:
        start_idx = min(n - 1, peak_idx + min_distance_points)
        stop = n
        sign_reason = "positive run"

    violation_run = 0
    violation_run_start = start_idx
    flat_run = 0
    flat_run_start = start_idx
    last_idx = start_idx
    prev_slope = None

    for i in range(start_idx, stop, direction):
        last_idx = i
        slope = float(dy[i])

        if prev_slope is not None:
            if direction > 0 and prev_slope < -eps and slope > eps:
                return i, "sign_reversal"
            if direction < 0 and prev_slope > eps and slope < -eps:
                return i, "sign_reversal"
        prev_slope = slope

        if abs(slope) <= eps:
            if flat_run == 0:
                flat_run_start = i
            flat_run += 1
            if flat_run >= k:
                return flat_run_start, "flat slope"
        else:
            flat_run = 0

        sign_violation = slope < -eps if direction < 0 else slope > eps
        if sign_violation:
            if violation_run == 0:
                violation_run_start = i
            violation_run += 1
            if violation_run >= k:
                return violation_run_start, sign_reason
        else:
            violation_run = 0

    return last_idx, "window edge"


def _extend_gradient_bound(y_sm, dy, peak_idx, start_idx, direction, eps,
                           intensity_fraction, flat_run_length,
                           positive_only=False,
                           exclude_at_idx=None):
    n = len(y_sm)
    if n == 0:
        return None, "window edge"

    flat_run_length = max(1, int(flat_run_length))
    peak_height = float(y_sm[peak_idx])
    threshold = peak_height * max(float(intensity_fraction), 0.0)
    stop = -1 if direction < 0 else n
    flat_run = 0
    flat_run_start = start_idx
    last_idx = start_idx

    for i in range(start_idx, stop, direction):
        last_idx = i
        if exclude_at_idx is not None and i == exclude_at_idx:
            boundary_idx = i + 1 if direction < 0 else i - 1
            return boundary_idx, "exclusive transition boundary"
        if positive_only and float(y_sm[i]) <= 0.0:
            if direction > 0:
                boundary_idx = max(0, i - 1)
                return boundary_idx, "non-positive intensity (exclusive)"
            boundary_idx = min(len(y_sm) - 1, i + 1)
            return boundary_idx, "non-positive intensity (exclusive)"
        if _is_local_minimum(y_sm, i):
            return i, "minimum"
        if threshold > 0 and float(y_sm[i]) <= threshold:
            return i, "threshold"

        slope = float(dy[i])
        if abs(slope) <= eps:
            if flat_run == 0:
                flat_run_start = i
            flat_run += 1
            if flat_run >= flat_run_length:
                return flat_run_start, "flat slope"
        else:
            flat_run = 0

    return last_idx, "window edge"


def _sanitize_sg_params(sg_window, sg_poly, n):
    if n < 3:
        return None, None

    max_window = n if n % 2 == 1 else n - 1
    if max_window < 3:
        return None, None

    sg_w = int(sg_window)
    sg_p = int(sg_poly)

    sg_w = min(sg_w, max_window)
    sg_w = max(sg_w, 3)
    if sg_w % 2 == 0:
        sg_w -= 1

    min_window = sg_p + 2 if sg_p % 2 == 0 else sg_p + 1
    sg_w = max(sg_w, min_window)
    if sg_w % 2 == 0:
        sg_w += 1
    if sg_w > max_window:
        sg_w = max_window
        if sg_w % 2 == 0:
            sg_w -= 1

    if sg_w < 3:
        return None, None

    sg_p = min(sg_p, sg_w - 1)
    return sg_w, sg_p


def get_gradient_bounds(x_arr, y_sm, dy, peak_idx, eps=None,
                        k=DEFAULT_GRADIENT_CONSECUTIVE_VIOLATIONS,
                        intensity_fraction=DEFAULT_INTENSITY_FRACTION,
                        min_distance_points=DEFAULT_GRADIENT_MIN_DISTANCE_POINTS,
                        flat_run_length=DEFAULT_GRADIENT_BASELINE_RUN_LENGTH,
                        positive_only=DEFAULT_GRADIENT_POSITIVE_ONLY,
                        return_debug=False):
    # find the end of the monotonic core, then walk out to the true edge
    n = len(dy)
    if n == 0:
        if return_debug:
            return None, None, {
                "peak_idx": int(peak_idx),
                "left": {
                    "transition_idx": None,
                    "transition_reason": "window edge",
                    "stop_idx": None,
                    "stop_reason": "window edge",
                },
                "right": {
                    "transition_idx": None,
                    "transition_reason": "window edge",
                    "stop_idx": None,
                    "stop_reason": "window edge",
                },
            }
        return None, None

    abs_dy = np.abs(dy)
    max_abs_dy = float(np.max(abs_dy)) if n else 0.0
    if eps is None:
        # Use the median of |dy| across the window as the gradient scale
        # reference.  The max is dominated by the steep peak core and inflates
        # eps so aggressively that real edge transitions get classified as flat.
        # The median is robust to those outlier peak slopes and stays close to
        # the baseline noise / shoulder gradient level.
        # Floor at 0.1% of the window-wide max so eps is never literally zero.
        median_abs_dy = float(np.median(abs_dy)) if n else 0.0
        eps = DEFAULT_GRADIENT_EPS_FRACTION * max(median_abs_dy, max_abs_dy * 0.001)

    eps = max(float(eps), 0.0)
    k = max(1, int(k))
    flat_run_length = max(1, int(flat_run_length))

    left_transition_idx, left_transition_reason = _find_gradient_transition(
        dy, peak_idx, -1, eps, k, min_distance_points
    )
    right_transition_idx, right_transition_reason = _find_gradient_transition(
        dy, peak_idx, +1, eps, k, min_distance_points
    )

    # sign_reversal, negative run (left), and positive run (right) all fire one
    # index *past* the local minimum at the peak base — the detection condition
    # triggers on the first point that has already moved away from that minimum.
    # Step back toward the peak by one so _extend_gradient_bound lands on the
    # minimum itself instead of overshooting to the next one on noisy data.
    _LEFT_OVERSHOOT  = ("sign_reversal", "negative run")
    _RIGHT_OVERSHOOT = ("sign_reversal", "positive run")
    left_exclude_idx = left_transition_idx if left_transition_reason in _LEFT_OVERSHOOT else None
    right_exclude_idx = right_transition_idx if right_transition_reason in _RIGHT_OVERSHOOT else None
    left_extend_start = (
        left_transition_idx + 1
        if left_transition_reason in _LEFT_OVERSHOOT and left_transition_idx < n - 1
        else left_transition_idx
    )
    right_extend_start = (
        right_transition_idx - 1
        if right_transition_reason in _RIGHT_OVERSHOOT and right_transition_idx > 0
        else right_transition_idx
    )

    left_idx, left_stop_reason = _extend_gradient_bound(
        y_sm, dy, peak_idx, left_extend_start, -1, eps,
        intensity_fraction, flat_run_length,
        positive_only=positive_only,
        exclude_at_idx=left_exclude_idx,
    )
    right_idx, right_stop_reason = _extend_gradient_bound(
        y_sm, dy, peak_idx, right_extend_start, +1, eps,
        intensity_fraction, flat_run_length,
        positive_only=positive_only,
        exclude_at_idx=right_exclude_idx,
    )

    if return_debug:
        return float(x_arr[left_idx]), float(x_arr[right_idx]), {
            "peak_idx": int(peak_idx),
            "left": {
                "transition_idx": int(left_transition_idx),
                "transition_reason": left_transition_reason,
                "stop_idx": int(left_idx),
                "stop_reason": left_stop_reason,
            },
            "right": {
                "transition_idx": int(right_transition_idx),
                "transition_reason": right_transition_reason,
                "stop_idx": int(right_idx),
                "stop_reason": right_stop_reason,
            },
        }

    return float(x_arr[left_idx]), float(x_arr[right_idx])


def analyze_window(df, target, half_window=DEFAULT_HALF_WIN,
                   sg_window=DEFAULT_SG_WINDOW, sg_poly=DEFAULT_SG_POLY,
                   use_smoothing=True,
                   intensity_fraction=DEFAULT_INTENSITY_FRACTION,
                   target_tolerance=DEFAULT_TARGET_TOLERANCE,
                   min_height_fraction=DEFAULT_MIN_HEIGHT_FRACTION,
                   use_prominence=False,
                   gradient_min_distance_points=DEFAULT_GRADIENT_MIN_DISTANCE_POINTS,
                   gradient_baseline_run_length=DEFAULT_GRADIENT_BASELINE_RUN_LENGTH,
                   gradient_positive_only=DEFAULT_GRADIENT_POSITIVE_ONLY):
    """
    Extract local window, smooth, compute derivatives, find auto peak + bounds.
    Returns result dict or None if insufficient data.
    """
    lo = target - half_window
    hi = target + half_window
    sub = df[(df["x"] >= lo) & (df["x"] <= hi)].copy().reset_index(drop=True)
    n = len(sub)
    if n < 4:
        return None

    x_arr = sub["x"].values
    y_raw = sub["y"].values

    if use_smoothing:
        sg_w, sg_p = _sanitize_sg_params(sg_window, sg_poly, n)

        try:
            if sg_w is None or sg_p is None:
                raise ValueError("Insufficient data for Savitzky-Golay filter")
            y_sm = savgol_filter(y_raw, window_length=sg_w, polyorder=sg_p)
        except Exception:
            y_sm = y_raw.copy()
    else:
        y_sm = y_raw.copy()

    dy = np.gradient(y_sm, x_arr)
    d2y = np.gradient(dy, x_arr)

    sub["smoothed"] = y_sm
    sub["dy"] = dy
    sub["d2y"] = d2y

    peak_kwargs = {}
    intensity_scale = float(np.max(y_sm)) if len(y_sm) else 0.0
    if min_height_fraction > 0 and intensity_scale > 0:
        threshold_value = float(min_height_fraction) * intensity_scale
        if use_prominence:
            peak_kwargs["prominence"] = threshold_value
        else:
            peak_kwargs["height"] = threshold_value

    candidate_idx, _peak_props = find_peaks(y_sm, **peak_kwargs)
    tolerance_mask = np.abs(x_arr[candidate_idx] - float(target)) <= float(target_tolerance)
    candidate_idx = candidate_idx[tolerance_mask]

    if len(candidate_idx):
        peak_idx = int(candidate_idx[np.argmin(np.abs(x_arr[candidate_idx] - float(target)))])
    else:
        # No candidate passed the tolerance filter.  argmax on noisy data will
        # land on a noise spike (highest point in window, not necessarily near
        # the target).  Instead, find all local maxima without constraints and
        # pick the one closest to the target m/z.  Fall back to argmax only if
        # there are genuinely no local maxima (monotone data).
        fallback_idx, _ = find_peaks(y_sm)
        if len(fallback_idx):
            peak_idx = int(fallback_idx[np.argmin(np.abs(x_arr[fallback_idx] - float(target)))])
        else:
            peak_idx = int(np.argmax(y_sm))

    auto_peak_mz = float(x_arr[peak_idx])
    auto_peak_int = float(y_sm[peak_idx])

    grad_left, grad_right, gradient_debug = get_gradient_bounds(
        x_arr,
        y_sm,
        dy,
        peak_idx,
        intensity_fraction=intensity_fraction,
        min_distance_points=gradient_min_distance_points,
        flat_run_length=gradient_baseline_run_length,
        positive_only=gradient_positive_only,
        return_debug=True,
    )
    bounds_fallback = grad_left >= grad_right
    if bounds_fallback:
        auto_left = float(x_arr[0])
        auto_right = float(x_arr[-1])
    else:
        auto_left = grad_left
        auto_right = grad_right

    return {
        "sub": sub,
        "x": x_arr,
        "y_raw": y_raw,
        "y_sm": y_sm,
        "dy": dy,
        "d2y": d2y,
        "auto_peak": auto_peak_mz,
        "auto_left": auto_left,
        "auto_right": auto_right,
        "grad_left": grad_left,
        "grad_right": grad_right,
        "gradient_peak_idx": gradient_debug["peak_idx"],
        "grad_left_transition_idx": gradient_debug["left"]["transition_idx"],
        "grad_left_transition_reason": gradient_debug["left"]["transition_reason"],
        "grad_left_stop_idx": gradient_debug["left"]["stop_idx"],
        "grad_left_stop_reason": gradient_debug["left"]["stop_reason"],
        "grad_right_transition_idx": gradient_debug["right"]["transition_idx"],
        "grad_right_transition_reason": gradient_debug["right"]["transition_reason"],
        "grad_right_stop_idx": gradient_debug["right"]["stop_idx"],
        "grad_right_stop_reason": gradient_debug["right"]["stop_reason"],
        "bounds_fallback": bounds_fallback,
        "auto_peak_int": auto_peak_int,
        "gradient_positive_only": bool(gradient_positive_only),
        "use_smoothing": bool(use_smoothing),
    }
