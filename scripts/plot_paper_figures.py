"""
Generate publication-quality KDE plots for the Eind and CRSS posteriors.

Three figure styles:
  - 2D KDE contour:  pairs of elastic constants (C11-C12, C11-C44, C12-C44).
                     Black contours (5 levels), no fill, 2"x2", grid, no labels.
  - 1D marginal:     one per (C11, C12, C44). Black line with gray fill, KDE
                     y-values scaled x100 (so they're visible on the same axis
                     range as the constants). 2"x2", grid, no labels.
  - CRSS marginal:   single larger plot (4"x4") with larger tick labels and
                     auto-scaled axes.

Filenames:
   ar_C12 vs C11.png       (2D contours)
   ar_C11.png  ar_C12.png  ar_C44.png    (1D marginals)
   ar_CRSS.png                            (CRSS)

Usage:
    # Generate all plots from the latest multichain results
    python scripts/plot_paper_figures.py

    # Specify CSVs explicitly
    python scripts/plot_paper_figures.py \
        --eind results/E_ind_sampler_ar.csv \
        --crss results/y_ind_sampler_ar.csv

    # Custom output folder and prefix
    python scripts/plot_paper_figures.py --output figures_v2 --prefix v2

    # Different axis range (default 0..375 for elastic constants)
    python scripts/plot_paper_figures.py --xlim-max 350
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT_DIR = ROOT / "figures"


# ---------------------------------------------------------------------------
# Plotting primitives
# ---------------------------------------------------------------------------

def plot_kde_2d(df, x, y, *, fig_width=2.0,
                xlim=(0, 375), ylim=(0, 375), tick_interval=100,
                levels=4, linewidth=1.0):
    """2D KDE contour plot — black contours, no fill, no axis labels."""
    fig, ax = plt.subplots(figsize=(fig_width, fig_width))

    # Newer seaborn: use `fill=False`; older seaborn: `shade=False`. Pass both
    # via kwargs so we work across versions.
    try:
        sns.kdeplot(x=df[x], y=df[y], color="black",
                    levels=levels, linewidths=linewidth, fill=False, ax=ax)
    except TypeError:
        sns.kdeplot(x=df[x], y=df[y], color="black",
                    levels=levels, linewidths=linewidth, shade=False, ax=ax)

    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    ax.set_xticks(np.arange(xlim[0], xlim[1], tick_interval))
    ax.set_yticks(np.arange(ylim[0], ylim[1], tick_interval))
    ax.grid(True, alpha=0.3)
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    fig.tight_layout()
    return fig, ax


def plot_kde_1d(df, x, *, fig_width=2.0,
                xlim=(0, 375), tick_interval=100, y_scale=100,
                linewidth=1.0):
    """1D KDE — black line, gray fill, KDE y-values multiplied by `y_scale`.

    The density values are scaled by `y_scale` (default 100) so they're
    visible alongside the actual GPa axis on the same plot.
    """
    fig, ax = plt.subplots(figsize=(fig_width, fig_width))

    # Compute the KDE curve by drawing an invisible one first and stealing
    # its line data.
    kde_artist = sns.kdeplot(
        x=df[x], ax=ax, color="black", linewidth=0, fill=False,
        clip=(xlim[0], xlim[1]),
    )
    line = kde_artist.lines[-1]   # most recently added line
    xs, ys = line.get_data()
    ys = ys * y_scale

    ax.plot(xs, ys, color="black", linewidth=linewidth)
    ax.fill_between(xs, ys, color="gray", alpha=0.3)

    ax.set_xlim(xlim[0], xlim[1])
    ax.set_xticks(np.arange(xlim[0], xlim[1], tick_interval))
    ax.grid(True, alpha=0.3)
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    fig.tight_layout()
    return fig, ax


def plot_kde_1d_crss(df, *, fig_width=4.0, x="CRSS", bw_adjust=1.0,
                     linewidth=1.0, tick_label_size=20):
    """CRSS-style 1D KDE — larger figure, larger ticks, auto-scaled axes."""
    fig, ax = plt.subplots(figsize=(fig_width, fig_width))

    kde_artist = sns.kdeplot(
        x=df[x], ax=ax, color="black", bw_adjust=bw_adjust,
        linewidth=0, fill=False,
    )
    line = kde_artist.lines[-1]
    xs, ys = line.get_data()

    ax.plot(xs, ys, color="black", linewidth=linewidth)
    ax.fill_between(xs, ys, color="gray", alpha=0.3)

    ax.grid(True, alpha=0.3)
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    ax.tick_params(axis="both", which="major", labelsize=tick_label_size)
    fig.tight_layout()
    return fig, ax


# ---------------------------------------------------------------------------
# High-level generators
# ---------------------------------------------------------------------------

ELASTIC_PAIRS = [("C12", "C11"), ("C44", "C11"), ("C44", "C12")]
ELASTIC_VARS  = ["C11", "C12", "C44"]


def generate_eind_plots(csv_path, output_dir, *, prefix="ar",
                        xlim=(0, 375), tick_interval=100, dpi=300):
    """Make the 3 + 3 elastic-constant plots for an Eind sampler CSV."""
    df = pd.read_csv(csv_path)
    print(f"Eind plots: loaded {len(df)} samples from {csv_path}")

    for x, y in ELASTIC_PAIRS:
        fig, _ = plot_kde_2d(df, x, y, xlim=xlim, ylim=xlim,
                             tick_interval=tick_interval)
        out = output_dir / f"{prefix}_{x} vs {y}.png"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  saved {out}")

    for x in ELASTIC_VARS:
        fig, _ = plot_kde_1d(df, x, xlim=xlim, tick_interval=tick_interval)
        out = output_dir / f"{prefix}_{x}.png"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  saved {out}")


def generate_crss_plots(csv_path, output_dir, *, prefix="ar", dpi=300):
    """Make the CRSS 1D KDE plot."""
    df = pd.read_csv(csv_path)
    if "CRSS" not in df.columns:
        raise ValueError(f"CRSS column not found in {csv_path}; "
                         f"columns are {df.columns.tolist()}")
    print(f"CRSS plot: loaded {len(df)} samples from {csv_path}")
    fig, _ = plot_kde_1d_crss(df)
    out = output_dir / f"{prefix}_CRSS.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out}")


# ---------------------------------------------------------------------------
# Default CSV resolution
# ---------------------------------------------------------------------------

def _find_default_eind_csv() -> Path | None:
    """Prefer multichain pooled, then single-chain, then original published."""
    candidates = [
        ROOT / "results" / "E_ind_sampler_ar_multichain.csv",
        ROOT / "results" / "E_ind_sampler_ar_reproduced.csv",
        ROOT / "results" / "E_ind_sampler_ar.csv",            # original published
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


def _find_default_crss_csv() -> Path | None:
    candidates = [
        ROOT / "results" / "y_ind_sampler_ar_multichain.csv",
        ROOT / "results" / "y_ind_sampler_ar_reproduced.csv",
        ROOT / "results" / "y_ind_sampler_ar.csv",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--eind", type=str, default=None,
                   help="Eind sampler CSV. If omitted, auto-picks the most "
                        "refined available result in results/.")
    p.add_argument("--crss", type=str, default=None,
                   help="CRSS sampler CSV. If omitted, auto-picks.")
    p.add_argument("--output", type=str, default=str(DEFAULT_OUTPUT_DIR),
                   help="Directory to save figures.")
    p.add_argument("--prefix", type=str, default="ar",
                   help="Filename prefix for the figures (e.g. 'ar' -> ar_C11.png).")
    p.add_argument("--xlim-max", type=float, default=375,
                   help="Upper bound of the elastic-constant axis range.")
    p.add_argument("--tick-interval", type=float, default=100,
                   help="Major tick spacing on elastic-constant axes.")
    p.add_argument("--dpi", type=int, default=300, help="Output figure DPI.")
    p.add_argument("--only", choices=["eind", "crss"], default=None,
                   help="Generate only one set of plots.")
    args = p.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Writing figures to {output_dir.resolve()}")

    # ---- Eind ----
    if args.only != "crss":
        eind_csv = Path(args.eind) if args.eind else _find_default_eind_csv()
        if eind_csv is None:
            print("(no Eind sampler CSV found, skipping elastic-constant plots)")
        elif not eind_csv.exists():
            print(f"(--eind {eind_csv} does not exist, skipping)")
        else:
            generate_eind_plots(
                eind_csv, output_dir,
                prefix=args.prefix,
                xlim=(0, args.xlim_max),
                tick_interval=args.tick_interval,
                dpi=args.dpi,
            )

    # ---- CRSS ----
    if args.only != "eind":
        crss_csv = Path(args.crss) if args.crss else _find_default_crss_csv()
        if crss_csv is None:
            print("(no CRSS sampler CSV found, skipping CRSS plot)")
        elif not crss_csv.exists():
            print(f"(--crss {crss_csv} does not exist, skipping)")
        else:
            generate_crss_plots(
                crss_csv, output_dir, prefix=args.prefix, dpi=args.dpi,
            )

    print("\nDone.")


if __name__ == "__main__":
    main()
