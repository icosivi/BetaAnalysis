#!/usr/bin/env python3
"""
dut_rotation.py
Estimate the rotation angle of a square DUT from the 2D tracker hit map.

Method
------
1. Build a 2D hit-density histogram (same binning as the ROOT Draw command).
2. Smooth + threshold → binary mask of the DUT footprint.
3. Compute the image gradient (Sobel); the gradient is perpendicular to each
   edge, so its angle + 90° gives the edge tangent direction.
4. Build a weighted histogram of tangent angles (0–180°), then fold it into
   0–90° (the two pairs of sides of a square are 90° apart).
5. The dominant peak = rotation angle of the DUT sides w.r.t. the x-axis.
6. Show:  raw map | angle histogram | rotation-corrected map.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.ndimage import gaussian_filter, sobel, rotate as ndimage_rotate

try:
    import uproot
    import awkward as ak
except ImportError:
    sys.exit("Need uproot and awkward:  pip install uproot awkward")

# ── Configuration ──────────────────────────────────────────────────────────────
ROOT_FILE   = "/media/SSD_4TB/TB11/stats/stats_Run3_tracks.root"
TREE_NAME   = "Analysis"

PMAX_CUT    = 20        # mV – minimum signal to tag a DUT hit
SIGNAL_CH   = range(1, 17)   # DUT channels (1-indexed, same as ROOT command)

# 2D histogram parameters (match the ROOT Draw command)
XBINS, XMIN, XMAX = 1800, -20.0, 10.0   # x_pos1
YBINS, YMIN, YMAX = 2400, -20.0, 20.0   # y_pos1

# Mask building: gaussian smoothing sigma [bins] and fraction of peak for threshold
SMOOTH_SIGMA = 5
MASK_THRESH  = 0.25

# Edge strength cut: only use gradient pixels above this fraction of the maximum
EDGE_THRESH  = 0.10

# Angle scan resolution
N_ANGLE_BINS = 1800    # → 0.1° resolution over 0–180°
# ──────────────────────────────────────────────────────────────────────────────


# ── 1. Load data ───────────────────────────────────────────────────────────────
print("Loading ROOT file …")
with uproot.open(ROOT_FILE) as f:
    tree = f[TREE_NAME]
    x    = tree["x_pos1"].array(library="np")
    y    = tree["y_pos1"].array(library="np")
    pmax = tree["pmax"].array(library="ak")   # variable-length vector per event

# Apply DUT signal cut: at least one channel above threshold
n_ch = ak.num(pmax, axis=1)
print(f"  Channels per event (min/max): {ak.min(n_ch)} / {ak.max(n_ch)}")

cut_ak = ak.any(pmax[:, list(SIGNAL_CH)] > PMAX_CUT, axis=1)
cut    = ak.to_numpy(cut_ak)

x_dut, y_dut = x[cut], y[cut]
print(f"  Events passing cut: {cut.sum():,} / {len(x):,}\n")


# ── 2. 2D histogram ────────────────────────────────────────────────────────────
H, xe, ye = np.histogram2d(x_dut, y_dut,
                            bins=[XBINS, YBINS],
                            range=[[XMIN, XMAX], [YMIN, YMAX]])
H = H.T                          # shape (YBINS, XBINS) so H[iy, ix]
xc = 0.5 * (xe[:-1] + xe[1:])   # x bin centres
yc = 0.5 * (ye[:-1] + ye[1:])   # y bin centres

# Physical pixel size [mm/bin]
dx = (XMAX - XMIN) / XBINS
dy = (YMAX - YMIN) / YBINS
assert abs(dx - dy) < 1e-6, "Non-square pixels – angle scale needs correction"


# ── 3. Binary mask of the DUT footprint ───────────────────────────────────────
H_sm = gaussian_filter(H.astype(float), sigma=SMOOTH_SIGMA)
mask = (H_sm > MASK_THRESH * H_sm.max()).astype(float)


# ── 4. Edge-orientation histogram ─────────────────────────────────────────────
# Sobel: axis=1 → d/dx (along columns), axis=0 → d/dy (along rows)
gx  = sobel(mask, axis=1)
gy  = sobel(mask, axis=0)
mag = np.hypot(gx, gy)

strong = mag > EDGE_THRESH * mag.max()
orient = np.arctan2(gy[strong], gx[strong])   # gradient normal direction, −π..π

# Tangent = perpendicular to gradient; fold to [0°, 180°)
tangent_deg = (np.degrees(orient) + 90.0) % 180.0

hist, bins = np.histogram(tangent_deg, bins=N_ANGLE_BINS, range=(0.0, 180.0),
                           weights=mag[strong])
bc = 0.5 * (bins[:-1] + bins[1:])    # bin centres

# Fold [0°,90°) ← combine angle θ and θ+90° (the two pairs of square sides)
half = N_ANGLE_BINS // 2
hist_fold = hist[:half] + hist[half:]
bc_fold   = bc[:half]

peak_idx   = np.argmax(hist_fold)
rot_angle  = bc_fold[peak_idx]

print(f"DUT rotation angle: {rot_angle:.2f}° (sides tilted from x-axis)")
print(f"To align DUT with axes: rotate coordinates by −{rot_angle:.2f}°\n")


# ── 5. Rotation-corrected map ──────────────────────────────────────────────────
H_corr = ndimage_rotate(H, rot_angle, reshape=False, order=1, mode='constant', cval=0)


# ── 6. Plots ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
fig.suptitle(f"DUT rotation analysis  —  angle = {rot_angle:.2f}°", fontsize=13)

norm_log = mcolors.LogNorm(vmin=1)

# Panel 1: raw 2D map
ax = axes[0]
im = ax.pcolormesh(xc, yc, np.where(H > 0, H, np.nan), cmap="inferno", norm=norm_log)
plt.colorbar(im, ax=ax, label="hits / bin")
ax.set_xlabel("x_pos1 (mm)")
ax.set_ylabel("y_pos1 (mm)")
ax.set_title("Raw 2D hit map")
ax.set_aspect("equal")

# Panel 2: edge-orientation histogram (folded 0–90°)
ax = axes[1]
ax.plot(bc_fold, hist_fold, lw=0.9, color="steelblue")
ax.axvline(rot_angle, color="red", lw=1.5, linestyle="--",
           label=f"peak = {rot_angle:.2f}°")
ax.set_xlabel("Edge tangent angle (°)")
ax.set_ylabel("Weighted edge count (a.u.)")
ax.set_title("Edge orientation (folded 0–90°)")
ax.legend(fontsize=10)
ax.set_xlim(0, 20)

# Panel 3: rotation-corrected map
ax = axes[2]
im2 = ax.pcolormesh(xc, yc, np.where(H_corr > 0.5, H_corr, np.nan),
                    cmap="inferno", norm=norm_log)
plt.colorbar(im2, ax=ax, label="hits / bin")
ax.set_xlabel("x (mm)")
ax.set_ylabel("y (mm)")
ax.set_title(f"Rotation-corrected (−{rot_angle:.2f}°)")
ax.set_aspect("equal")

plt.tight_layout()
out = "dut_rotation.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.show()
print(f"Plot saved → {out}")
