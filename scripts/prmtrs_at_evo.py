# scripts/prmtrs_at_evo.py
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, PercentFormatter
from scipy.interpolate import make_interp_spline

# ----------------------------
# Config
# ----------------------------
SAVE = True
OUT_DIR = "figures"
FILENAME = "at_prmtrs_evo"
FORMATS = ("pdf", "svg", "png")   # choose any subset
DPI = 300

USE_TEX = True        # fallback to False automatically if LaTeX not found
FONT_SIZE = 11        # base font size (axes labels/ticks scale from this)
LINE_W = 2.5
MARKER_SZ = 5
SMOOTH_POINTS = 200   # resolution for spline curves

# ----------------------------
# Data (experimental points)
# ----------------------------
x_ina   = np.array([0.50, 0.65, 0.82, 1.00])
y_ina   = np.array([95.0, 76.0, 71.3, 70.0])

x_inal  = np.array([1.00, 2.80, 3.00])
y_inal  = np.array([101.0, 76.0, 73.0])

x_inak  = np.array([0.33, 0.50, 0.72, 1.00])
y_inak  = np.array([70.0, 76.0, 81.7, 85.0])

x_omega = np.array([0.05, 0.18, 0.24, 0.34, 0.57, 1.00])
y_omega = np.array([110.0, 76.0, 73.0, 72.0, 71.3, 70.0])

# Labels & axis settings for each subplot
panels = [
    dict(
        key="omega", x=x_omega, y=y_omega,
        xlabel=r"$\Omega$ (relative, % of control)",
        xlim=(0.0, 1.0), xticks=np.arange(0, 1.0 + 0.5, 0.5),
        format_percent=True,
        title=r"(A) Tissue resistivity, $\Omega$"
    ),
    dict(
        key="inak", x=x_inak, y=y_inak,
        xlabel=r"$P_{\mathrm{Na/K}}$ (relative to control)",
        xlim=(0.0, 1.0), xticks=np.arange(0, 1.0 + 0.5, 0.5),
        title=r"(B) $\mathrm{Na^+/K^+}$ pump, $P_{\mathrm{Na/K}}$"
    ),
    dict(
        key="ina", x=x_ina, y=y_ina,
        xlabel=r"$G_{\mathrm{Na}}$ (relative to control)",
        xlim=(0.0, 1.0), xticks=np.arange(0, 1.0 + 0.5, 0.5),
        title=r"(C) Fast sodium, $G_{\mathrm{Na}}$"
    ),
    dict(
        key="inal", x=x_inal, y=y_inal,
        xlabel=r"$G_{\mathrm{NaL}}$ (relative to control)",
        xlim=(0.0, 3.0), xticks=np.arange(0, 3.0 + 1.0, 1.0),
        title=r"(D) Late sodium, $G_{\mathrm{NaL}}$"
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

# Color cycle
colors = plt.rcParams["axes.prop_cycle"].by_key().get("color", None)
if colors is None:
    colors = ["#4C78A8", "#F58518", "#54A24B", "#E45756"]

# ----------------------------
# Helpers
# ----------------------------
def smooth_curve(x, y, n=SMOOTH_POINTS, k=2):
    """Return smoothed (x_new, y_new) via spline. Requires len(x) > k."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if len(x) <= k:
        # Not enough points for this spline degree; just linearly interpolate
        x_new = np.linspace(x.min(), x.max(), n)
        y_new = np.interp(x_new, x, y)
        return x_new, y_new
    # Ensure x is strictly increasing for spline
    order = np.argsort(x)
    x_sorted, y_sorted = x[order], y[order]
    spline = make_interp_spline(x_sorted, y_sorted, k=min(k, len(x_sorted)-1))
    x_new = np.linspace(x_sorted.min(), x_sorted.max(), n)
    y_new = spline(x_new)
    return x_new, y_new

# ----------------------------
# Plot
# ----------------------------
fig, axes = plt.subplots(1, 4, figsize=(10.2, 3.2), sharey=True)
ylim = (60, 120)

for i, spec in enumerate(panels):
    ax = axes[i]
    x, y = spec["x"], spec["y"]

    # raw points
    ax.plot(x, y, linestyle="none", marker="o", markersize=MARKER_SZ,
            markerfacecolor="white", markeredgecolor=colors[i % len(colors)],
            markeredgewidth=1.2, zorder=3, label="data")

    # smooth fit
    xs, ys = smooth_curve(x, y, n=SMOOTH_POINTS, k=2)
    ax.plot(xs, ys, linewidth=LINE_W, color=colors[i % len(colors)], label="spline")

    # axes cosmetics
    ax.set_title(spec["title"], pad=8)
    ax.set_xlabel(spec["xlabel"])
    if i == 0:
        ax.set_ylabel(r"$\mathrm{AT}_{\max}$ (ms)")
    ax.set_ylim(ylim)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=False, prune=None))
    ax.set_xlim(spec["xlim"])
    ax.set_xticks(spec["xticks"])

    if spec.get("format_percent", False):
        # show x as percent of control (0..1 -> 0..100%)
        ax.xaxis.set_major_formatter(PercentFormatter(xmax=1.0, decimals=0))

    # Clean spines & grid
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(True, which="major", axis="both", alpha=0.2, linewidth=0.8)

# Align layout nicely
fig.tight_layout(w_pad=2.0)

# Optional single legend (small, unobtrusive)
# Place a single legend across the bottom (only once)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncols=2, frameon=False, bbox_to_anchor=(0.5, -0.02))

# Save
if SAVE:
    os.makedirs(OUT_DIR, exist_ok=True)
    meta = {
        "Title": "AT vs parameter evolution",
    }
    for ext in FORMATS:
        path = os.path.join(OUT_DIR, f"{FILENAME}.{ext}")
        fig.savefig(path, dpi=DPI, bbox_inches="tight", transparent=True, metadata=meta)

plt.show()
