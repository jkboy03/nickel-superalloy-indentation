"""
Final-results analysis.

Loads the pooled multichain CSVs for the Eind (C11, C12, C44) and CRSS
calibrations and prints the publication-ready summary:

  * mean and std per parameter (and 3*std uncertainty bound)
  * full 3x3 covariance matrix on (C11, C12, C44)
  * MAP estimate (row with maximum SampleLogFunc)
  * derived quantities (anisotropy A, bulk modulus K)
  * CRSS in both GPa and MPa

Usage
-----
    # Default: auto-pick most refined results
    python scripts/03_analyze_results.py

    # Specify CSVs explicitly
    python scripts/03_analyze_results.py \\
        --eind results/E_ind_sampler_ar_multichain_stage1.csv \\
        --crss results/y_ind_sampler_ar_multichain_stage1.csv

    # Apply burnin trimming (default 0 — included for flexibility)
    python scripts/03_analyze_results.py --burnin 1000

    # Save the full report to a Markdown file
    python scripts/03_analyze_results.py --output results/summary.md
"""
import argparse
import importlib.util
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = Path(__file__).resolve().parent
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"

ELASTIC_VARS = ["C11", "C12", "C44"]


# Load plot_paper_figures.py as a module (filename doesn't start with a letter,
# so we can't use `import`).
_spec = importlib.util.spec_from_file_location(
    "ppf", SCRIPTS_DIR / "plot_paper_figures.py"
)
ppf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ppf)


# ---------------------------------------------------------------------------
# CSV resolution
# ---------------------------------------------------------------------------

def _find_eind_csv():
    for f in [
        RESULTS_DIR / "E_ind_sampler_ar_multichain.csv",
        RESULTS_DIR / "E_ind_sampler_ar_reproduced.csv",
        RESULTS_DIR / "E_ind_sampler_ar.csv",
    ]:
        if f.exists():
            return f
    return None


def _find_crss_csv():
    for f in [
        RESULTS_DIR / "y_ind_sampler_ar_multichain.csv",
        RESULTS_DIR / "y_ind_sampler_ar_reproduced.csv",
        RESULTS_DIR / "y_ind_sampler_ar.csv",
    ]:
        if f.exists():
            return f
    return None


# ---------------------------------------------------------------------------
# Summary helpers
# ---------------------------------------------------------------------------

def _mean_std(arr):
    """np.std with biased estimator (ddof=0)."""
    return float(np.mean(arr)), float(np.std(arr))


def analyze_elastic(df, out, burnin: int = 0):
    """Print and return the elastic-constant summary."""
    df = df.iloc[burnin:].copy()
    n = len(df)

    out.write(f"\n{'=' * 70}\n")
    out.write(f"Elastic constants  (n = {n}, burnin = {burnin})\n")
    out.write("=" * 70 + "\n")

    header = f"{'param':>6}  {'mean':>9}  {'std':>8}  {'3*std':>8}  " \
             f"{'min':>9}  {'max':>9}\n"
    out.write(header)
    out.write("-" * len(header) + "\n")
    stats = {}
    for v in ELASTIC_VARS:
        m, s = _mean_std(df[v])
        lo, hi = float(df[v].min()), float(df[v].max())
        stats[v] = {"mean": m, "std": s}
        out.write(f"{v:>6}  {m:9.3f}  {s:8.3f}  {3*s:8.3f}  "
                  f"{lo:9.3f}  {hi:9.3f}\n")

    # Derived quantities — anisotropy ratio A and bulk modulus K
    A = 2.0 * df["C44"] / (df["C11"] - df["C12"])
    K = (df["C11"] + 2.0 * df["C12"]) / 3.0
    m_A, s_A = _mean_std(A)
    m_K, s_K = _mean_std(K)
    out.write(f"{'A':>6}  {m_A:9.3f}  {s_A:8.3f}  {3*s_A:8.3f}  "
              f"{float(A.min()):9.3f}  {float(A.max()):9.3f}\n")
    out.write(f"{'K':>6}  {m_K:9.3f}  {s_K:8.3f}  {3*s_K:8.3f}  "
              f"{float(K.min()):9.3f}  {float(K.max()):9.3f}\n")
    stats["A"] = {"mean": m_A, "std": s_A}
    stats["K"] = {"mean": m_K, "std": s_K}

    # Covariance matrix
    cov = np.cov([df["C11"].to_numpy(), df["C12"].to_numpy(),
                  df["C44"].to_numpy()])
    out.write("\nCovariance matrix (C11, C12, C44):\n")
    out.write(f"{'':>8}  {'C11':>10}  {'C12':>10}  {'C44':>10}\n")
    for i, v in enumerate(ELASTIC_VARS):
        out.write(f"{v:>8}  ")
        for j in range(3):
            out.write(f"{cov[i, j]:10.3f}  ")
        out.write("\n")

    # Correlation matrix (more interpretable than covariance)
    diag = np.sqrt(np.diag(cov))
    corr = cov / np.outer(diag, diag)
    out.write("\nCorrelation matrix (C11, C12, C44):\n")
    out.write(f"{'':>8}  {'C11':>10}  {'C12':>10}  {'C44':>10}\n")
    for i, v in enumerate(ELASTIC_VARS):
        out.write(f"{v:>8}  ")
        for j in range(3):
            out.write(f"{corr[i, j]:10.4f}  ")
        out.write("\n")

    # MAP estimate
    if "SampleLogFunc" in df.columns:
        idx = df["SampleLogFunc"].idxmax()
        out.write(f"\nMAP estimate (row of max SampleLogFunc = "
                  f"{df.loc[idx, 'SampleLogFunc']:.4f}):\n")
        for v in ELASTIC_VARS:
            out.write(f"  {v:>4} = {float(df.loc[idx, v]):.3f}\n")
        out.write(f"  {'A':>4} = "
                  f"{2*df.loc[idx, 'C44']/(df.loc[idx, 'C11']-df.loc[idx, 'C12']):.4f}\n")
        stats["MAP"] = {v: float(df.loc[idx, v]) for v in ELASTIC_VARS}

    return stats


def analyze_crss(df, out, burnin: int = 0):
    """Print and return the CRSS summary."""
    df = df.iloc[burnin:].copy()
    n = len(df)

    out.write(f"\n{'=' * 70}\n")
    out.write(f"CRSS  (n = {n}, burnin = {burnin})\n")
    out.write("=" * 70 + "\n")

    m, s = _mean_std(df["CRSS"])
    out.write(f"  mean    = {m:.6f}  ({m*1000:.2f} MPa)\n")
    out.write(f"  std     = {s:.6f}  ({s*1000:.2f} MPa)\n")
    out.write(f"  3*std   = {3*s:.6f}  ({3*s*1000:.2f} MPa)\n")
    out.write(f"  median  = {float(df['CRSS'].median()):.6f}\n")
    out.write(f"  min/max = {float(df['CRSS'].min()):.6f} / "
              f"{float(df['CRSS'].max()):.6f}\n")

    if "SampleLogFunc" in df.columns:
        idx = df["SampleLogFunc"].idxmax()
        out.write(f"\nMAP CRSS (max SampleLogFunc = "
                  f"{df.loc[idx, 'SampleLogFunc']:.4f}): "
                  f"{float(df.loc[idx, 'CRSS']):.6f}  "
                  f"({float(df.loc[idx, 'CRSS'])*1000:.2f} MPa)\n")

    return {"mean": m, "std": s, "mean_MPa": m * 1000, "std_MPa": s * 1000}


def render_paper_block(elastic_stats, crss_stats, out):
    """Print a publication-ready one-line-per-parameter block."""
    out.write(f"\n{'=' * 70}\n")
    out.write("Publication-ready summary\n")
    out.write("=" * 70 + "\n")
    if elastic_stats:
        out.write(f"  C11  = {elastic_stats['C11']['mean']:7.2f} +/- "
                  f"{elastic_stats['C11']['std']:5.2f} GPa\n")
        out.write(f"  C12  = {elastic_stats['C12']['mean']:7.2f} +/- "
                  f"{elastic_stats['C12']['std']:5.2f} GPa\n")
        out.write(f"  C44  = {elastic_stats['C44']['mean']:7.2f} +/- "
                  f"{elastic_stats['C44']['std']:5.2f} GPa\n")
        out.write(f"  A    = {elastic_stats['A']['mean']:7.3f} +/- "
                  f"{elastic_stats['A']['std']:5.3f}\n")
        out.write(f"  K    = {elastic_stats['K']['mean']:7.2f} +/- "
                  f"{elastic_stats['K']['std']:5.2f} GPa\n")
    if crss_stats:
        out.write(f"  CRSS = {crss_stats['mean_MPa']:7.2f} +/- "
                  f"{crss_stats['std_MPa']:5.2f} MPa\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--eind", type=str, default=None,
                   help="Eind sampler CSV (auto-picks most refined if omitted).")
    p.add_argument("--crss", type=str, default=None,
                   help="CRSS sampler CSV (auto-picks if omitted).")
    p.add_argument("--burnin", type=int, default=0,
                   help="Number of initial samples to discard.")
    p.add_argument("--output", type=str, default=None,
                   help="Optional file to also write the report into.")
    p.add_argument("--only", choices=["eind", "crss"], default=None,
                   help="Analyze only one calibration.")
    p.add_argument("--no-plots", action="store_true",
                   help="Skip generating the publication figures.")
    p.add_argument("--show-plots", action="store_true",
                   help="Also pop up the plots interactively (in addition "
                        "to saving them).")
    p.add_argument("--figures-dir", type=str, default=str(FIGURES_DIR),
                   help="Directory to save publication PNG figures into.")
    p.add_argument("--fig-prefix", type=str, default="ar",
                   help="Prefix for figure filenames (e.g. 'ar' -> ar_C11.png).")
    p.add_argument("--xlim-max", type=float, default=375,
                   help="Upper bound of the elastic-constant axis range.")
    p.add_argument("--tick-interval", type=float, default=100,
                   help="Major tick spacing on elastic-constant axes.")
    p.add_argument("--dpi", type=int, default=300, help="Output figure DPI.")
    args = p.parse_args()

    buf = io.StringIO()
    elastic_stats = None
    crss_stats = None
    eind_csv = None
    crss_csv = None

    if args.only != "crss":
        eind_csv = Path(args.eind) if args.eind else _find_eind_csv()
        if eind_csv is None or not eind_csv.exists():
            buf.write(f"(no Eind CSV found, skipping)\n")
            eind_csv = None
        else:
            buf.write(f"Loading Eind chain: {eind_csv}\n")
            elastic_stats = analyze_elastic(pd.read_csv(eind_csv), buf, burnin=args.burnin)

    if args.only != "eind":
        crss_csv = Path(args.crss) if args.crss else _find_crss_csv()
        if crss_csv is None or not crss_csv.exists():
            buf.write(f"(no CRSS CSV found, skipping)\n")
            crss_csv = None
        else:
            buf.write(f"Loading CRSS chain: {crss_csv}\n")
            crss_stats = analyze_crss(pd.read_csv(crss_csv), buf, burnin=args.burnin)

    render_paper_block(elastic_stats, crss_stats, buf)

    report = buf.getvalue()
    print(report)

    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(report)
        print(f"\nReport saved to {out_path}")

    # ---- Publication figures ----
    if not args.no_plots:
        fig_dir = Path(args.figures_dir)
        fig_dir.mkdir(parents=True, exist_ok=True)
        print(f"\nGenerating publication figures into {fig_dir}/ ...")
        if eind_csv is not None:
            ppf.generate_eind_plots(
                eind_csv, fig_dir,
                prefix=args.fig_prefix,
                xlim=(0, args.xlim_max),
                tick_interval=args.tick_interval,
                dpi=args.dpi,
            )
        if crss_csv is not None:
            ppf.generate_crss_plots(
                crss_csv, fig_dir, prefix=args.fig_prefix, dpi=args.dpi,
            )
        if args.show_plots:
            import matplotlib.pyplot as plt
            print("\nPopping up figures (close each window to advance)...")
            for png in sorted(fig_dir.glob(f"{args.fig_prefix}_*.png")):
                img = plt.imread(png)
                fig, ax = plt.subplots(figsize=(6, 6))
                ax.imshow(img)
                ax.set_title(png.name)
                ax.axis("off")
                plt.show()


if __name__ == "__main__":
    main()
