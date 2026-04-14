import pandas as pd 
import os
import numpy as np 
import matplotlib.pyplot as plt

from pathlib import Path
_root = Path(__file__).resolve().parent.parent
df = pd.read_csv(_root / "data" / "samples" / "neodynium.csv")

print(df.head())

df.columns = ['x', 'y']

x = df['x'].values
y = df['y'].values

# compute deriv
dy_dx = np.gradient(y, x)
df['dy_dx'] = dy_dx

#2nd deriv
d2y_dx2 = np.gradient(dy_dx, x)
df['d2y_dx2'] = d2y_dx2

# filter spectral range
xmin, xmax = 157.8, 158.0
df_range = df[(df['x'] >= xmin) & (df['x'] <= xmax)]

# plot
fig, ax1 = plt.subplots()

# left axis (original signal)
ax1.plot(df_range['x'], df_range['y'], label='Original')
ax1.set_xlabel("m/z")
ax1.set_ylabel("Intensity")

# right axis (derivative)
ax2 = ax1.twinx()
ax2.plot(df_range['x'], df_range['dy_dx'], linestyle='--', label='Derivative')
ax2.set_ylabel("d(Intensity)/dx")

# right axis (derivative)
ax3 = ax1.twinx()
ax3.plot(df_range['x'], df_range['d2y_dx2'], linestyle='--', label='2nd Derivative')
ax3.set_ylabel("d2(Intensity)/dx2")

# combine legends
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
lines3, labels3 = ax3.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2)

plt.show()