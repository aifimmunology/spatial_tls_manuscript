# Scratch snippet: pick x_clip / y_clip that keep ~90% of the CD4 cells.
# Paste into the notebook, then delete this file (it's not part of the pipeline).

# adata_cd4 = your CD4 subset, e.g.
#   adata_cd4 = adata[adata.obs['label_medium'] == 'CD4 act']

x_axis = 'avg_distance_to_bronchi_zone'
y_axis = 'distance_to_tls'
keep_frac = 0.90

# --- per-axis: each clip keeps `keep_frac` of cells ALONG THAT AXIS alone ---
x_clip = float(np.nanpercentile(adata_cd4.obs[x_axis], 100 * keep_frac))
y_clip = float(np.nanpercentile(adata_cd4.obs[y_axis], 100 * keep_frac))

print(f"x_clip ({x_axis}): {x_clip:.1f}  -> keeps {keep_frac:.0%} of cells along x")
print(f"y_clip ({y_axis}): {y_clip:.1f}  -> keeps {keep_frac:.0%} of cells along y")

# A cell is dropped if it's beyond EITHER clip, so the 2D window keeps fewer
# than `keep_frac`. Check the fraction that actually survives both clips:
m = (adata_cd4.obs[x_axis] <= x_clip) & (adata_cd4.obs[y_axis] <= y_clip)
print(f"{m.mean():.1%} of CD4 cells fall inside the {x_clip:.0f} x {y_clip:.0f} window")

# --- alternative: make the 2D WINDOW keep ~`keep_frac` of cells ---
# clip each axis at sqrt(keep_frac) so the joint window retains ~keep_frac.
per_axis = np.sqrt(keep_frac)            # e.g. 0.949 for 90%
x_clip_joint = float(np.nanpercentile(adata_cd4.obs[x_axis], 100 * per_axis))
y_clip_joint = float(np.nanpercentile(adata_cd4.obs[y_axis], 100 * per_axis))
m_joint = (adata_cd4.obs[x_axis] <= x_clip_joint) & (adata_cd4.obs[y_axis] <= y_clip_joint)
print(f"joint clips {x_clip_joint:.0f} x {y_clip_joint:.0f} keep {m_joint.mean():.1%} of CD4 cells")

# ---------------------------------------------------------------------------
# Shape-based upper limits (when you'd rather set the clip from the spread of
# the distribution than from a target percentile).
#
# These distances are non-negative and right-skewed, so plain mean + k*std is
# inflated by the long tail of far-from-structure cells. The two rules below
# are robust to that skew. Each prints the fraction kept, so you can tune k to
# the coverage you want -- which is really the underlying principled quantity.
# ---------------------------------------------------------------------------
k = 2.0  # number of "spreads" out to extend each axis

# (a) Robust: median + k * scaled-MAD. 1.4826 * MAD == std for a Gaussian, so k
#     keeps its usual ~normal meaning while ignoring the far-cell tail.
print("\nrobust (median + k*MAD):")
for axis in [x_axis, y_axis]:
    v = adata_cd4.obs[axis].to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    med = np.median(v)
    mad = np.median(np.abs(v - med)) * 1.4826
    hi = med + k * mad
    print(f"  {axis}: clip = {hi:.1f}  (median {med:.1f} + {k}*MAD {mad:.1f}), "
          f"keeps {(v <= hi).mean():.1%}")

# (b) Log-normal: fit mean/std on log1p(distance), then back-transform. Use if
#     the distances look roughly log-normal (multiplicative spread).
print("log-normal (back-transformed mean + k*std on log1p):")
for axis in [x_axis, y_axis]:
    v = adata_cd4.obs[axis].to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    lv = np.log1p(v)
    hi = np.expm1(lv.mean() + k * lv.std())
    print(f"  {axis}: clip = {hi:.1f}, keeps {(v <= hi).mean():.1%}")
