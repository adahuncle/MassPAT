import numpy as np

DEFAULT_HALF_WIN = 0.15
DEFAULT_TARGET_TOLERANCE = 0.005

def analyze_window_simple(df, target, half_window=DEFAULT_HALF_WIN, target_tolerance=DEFAULT_TARGET_TOLERANCE):
    """Cut window, find peak near target, walk to 1% height boundaries."""
    
    lo, hi = target - half_window, target + half_window
    sub = df[(df["x"] >= lo) & (df["x"] <= hi)].copy().reset_index(drop=True)
    
    if len(sub) < 5:
        return None
    
    x, y = sub["x"].to_numpy(), sub["y"].to_numpy()
    
    # Find local maxima
    peaks = [i for i in range(1, len(y)-1) if y[i] >= y[i-1] and y[i] >= y[i+1]]
    
    # Select peak closest to target
    close = [i for i in peaks if abs(x[i] - target) <= target_tolerance]
    if close:
        peak_idx = min(close, key=lambda i: abs(x[i] - target))
    elif peaks:
        peak_idx = min(peaks, key=lambda i: abs(x[i] - target))
    else:
        peak_idx = int(np.argmax(y))
    
    # Walk to 1% height boundaries
    min_h = y[peak_idx] * 0.01
    
    left = peak_idx
    while left > 0 and y[left-1] <= y[left] and y[left-1] > min_h:
        left -= 1
    
    right = peak_idx
    while right < len(y)-1 and y[right+1] <= y[right] and y[right] > min_h:
        right += 1
    
    return {
        "sub": sub, "x": x, "y_raw": y, "y_sm": y,
        "auto_peak": float(x[peak_idx]), "auto_left": float(x[left]), "auto_right": float(x[right]),
        "auto_peak_int": float(y[peak_idx]), "peak_idx": int(peak_idx),
        "left_idx": int(left), "right_idx": int(right),
    }