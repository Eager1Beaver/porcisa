# scripts/prmtrs_dr_evo.py
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import make_interp_spline

# ----------------------------
# Config
# ----------------------------
SAVE = True
OUT_DIR = "figures"
FILENAME = "dr_prmtrs_evo"
FORMATS = ("pdf", "svg", "png")
DPI = 300

USE_TEX = True         # will fall back automatically if LaTeX isn't found
FONT_SIZE = 11
LINE_W = 2.5
MARKER_SZ = 5
SMOOTH_POINTS = 200

# ----------------------------
# Data (experimental points)
# ----------------------------
x_k    = np.array([5.0, 6.8, 8.1, 8.5])
y_k    = np.array([105.0, 106.0, 110.0, 113.0])

x_inal = np.array([1.0, 2.14, 2.8, 3.0])
y_inal = np.array([115.0, 114.0, 110.0, 107.0])

x_ical = np.array([0.5, 0.65, 0.83, 1.0])
y_ical = np.array([115.0, 110.0, 102.0, 90.0])

x_ikatp = np.array([0.0, 0.08, 0.1])
y_ikatp = np.array([85.0, 110.0, 116.0])

panels = [
    dict(
        key="k_o", x=x_k, y=y_k,
        xlabel=r"Extracellular $[K^+]$ (mM)",
        title=r"(A) $[K^+]_{o}$",
        xlim=(4.0, 9.0), xticks=np.arange(4.0, 9.0 + 1.0, 1.0),
    ),
    dict(
        key="g_nal", x=x_inal, y=y_inal,
        xlabel=r"$G_{\mathrm{NaL}}$ (relative to control)",
        title=r"(B) Late sodium, $G_{\mathrm{NaL}}$",
        xlim=(0.0, 3.0), xticks=np.arange(0.0, 3.0 + 1.0, 1.0),
    ),
    dict(
        key="g_cal", x=x_ical, y=y_ical,
        xlabel=r"$G_{\mathrm{CaL}}$ (relative to control)",
        title=r"(C) L-type calcium, $G_{\mathrm{CaL}}$",
        xlim=(0.0, 1.0), xticks=np.arange(0.0, 1.0 + 0.5, 0.5),
    ),
    dict(
        key="f_katp", x=x_ikatp, y=y_ikatp,
        xlabel=r"$f_{\mathrm{K(ATP)}}$ (open fraction)",
        title=r"(D) K(ATP) open fraction, $f_{\mathrm{K(ATP)}}$",
        xlim=(0.0, 0.10), xticks=np.arange(0.0, 0.10 + 0.05, 0.05),
    ),
]

# ----------------------------
# Styling
# ----------------------------
try:
    if USE_TEX:
        plt.rcParams["text.usetex"] = True
        plt.rcParams["font.family"] = "serif"
        plt.rcParams["font.serif"] = ["Computer Modern Roman"]
except Exception:
    plt.rcParams["text.usetex"] = False

plt.rcParams.update({
    "axes.titlesize": FONT_SIZE,
    "axes.labelsize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE - 1,
    "ytick.labelsize": FONT_SIZE - 1,
    "axes.linewidth": 1.0,
    "figure.dpi": 100,
})

colors = plt.rcParams["axes.prop_cycle"].by_key().get("color", ["#4C78A8", "#F58518", "#54A24B", "#E45756"])

# ----------------------------
# Helpers
# ----------------------------
def smooth_curve(x, y, n=SMOOTH_POINTS, k=2):
    """Smoothed (x_new, y_new) via spline; falls back to linear interp if few points."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 2:
        return x, y
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    # Spline degree capped by number of points - 1
    k = min(k, len(xs) - 1)
    if k < 1:
        return xs, ys
    try:
        spline = make_interp_spline(xs, ys, k=k)
        x_new = np.linspace(xs.min(), xs.max(), n)
        y_new = spline(x_new)
        return x_new, y_new
    except Exception:
        # fallback to piecewise linear if spline fails
        x_new = np.linspace(xs.min(), xs.max(), n)
        y_new = np.interp(x_new, xs, ys)
        return x_new, y_new

# ----------------------------
# Plot
# ----------------------------
fig, axes = plt.subplots(1, 4, figsize=(10.2, 3.2), sharey=True)
ylim = (60, 120)

for i, spec in enumerate(panels):
    ax = axes[i]
    x, y = spec["x"], spec["y"]

    # raw markers
    ax.plot(
        x, y,
        linestyle="none",
        marker="o",
        markersize=MARKER_SZ,
        markerfacecolor="white",
        markeredgecolor=colors[i % len(colors)],
        markeredgewidth=1.2,
        zorder=3,
        label="data",
    )

    # smooth line
    xs, ys = smooth_curve(x, y, n=SMOOTH_POINTS, k=2)
    ax.plot(xs, ys, linewidth=LINE_W, color=colors[i % len(colors)], label="spline")

    # cosmetics
    ax.set_title(spec["title"], pad=8)
    ax.set_xlabel(spec["xlabel"])
    if i == 0:
        ax.set_ylabel(r"$\mathrm{DR}$ (ms)")
    ax.set_ylim(ylim)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.set_xlim(spec["xlim"])
    ax.set_xticks(spec["xticks"])

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(True, which="major", axis="both", alpha=0.2, linewidth=0.8)

fig.tight_layout(w_pad=2.0)

# Unobtrusive, shared legend (optional)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncols=2, frameon=False, bbox_to_anchor=(0.5, -0.02))

# Save
if SAVE:
    os.makedirs(OUT_DIR, exist_ok=True)
    meta = {
        "Title": "Depolarization vs parameter evolution",
    }
    for ext in FORMATS:
        fig.savefig(
            os.path.join(OUT_DIR, f"{FILENAME}.{ext}"),
            dpi=DPI, bbox_inches="tight", transparent=True, metadata=meta
        )

plt.show()
