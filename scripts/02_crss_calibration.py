"""
Bayesian calibration of the critical resolved shear stress (CRSS) for
FCC Ni-based superalloys using a Gaussian Process orientation factor
surrogate and ParaDRAM MCMC.

Model:
    measured_Yind(orientation) = CRSS * factor(orientation)

The factor(orientation) is learned from FEM simulations performed at a
reference shear stress s_0_sim = 0.015. The MCMC then estimates CRSS by
matching the predicted yield stresses against the experimental ones.
"""

# === Standard library ===
import os
import sys
import json
import random
import pickle
import argparse
import subprocess
from pathlib import Path

# === Scientific stack ===
import numpy as np
import pandas as pd
import scipy.io

# === ML / surrogate ===
import GPy
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

# === Sampler ===
import paramonte as pm

# === Plotting ===
import matplotlib.pyplot as plt

# === Local package: GSH utilities ===
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from gshutils import gshUtils

# === Paths ===
DATA_SIM_DIR = ROOT / "data" / "simulation_crss"
DATA_EXP_DIR = ROOT / "data" / "experimental"
RESULTS_DIR  = ROOT / "results"
CACHE_DIR    = ROOT / "cache"
SURROGATE_CACHE = CACHE_DIR / "surrogate_crss.pkl"  # delete to retrain
PARADRAM_OUT_DIR = RESULTS_DIR / "paradram_out_crss"

# Reference shear stress used to normalize the FEM simulations.
# Yields are divided by this to get the orientation factor.
S_0_SIM = 0.015


# ---------------------------------------------------------------------------
# Step 1 — Load CRSS simulation dataset
# ---------------------------------------------------------------------------
# crss_simulation.mat contains the 18-row training set. Fields (as the FEM
# simulation produced them):
#   PHI, phi2   — Euler angles of the indented crystal (radians)
#   yind        — raw yield stress per orientation (before normalization)
# SSH descriptors (l4m1, l12m2) are computed at load time via gshutils.
# ---------------------------------------------------------------------------

def load_crss_simulation_data() -> pd.DataFrame:
    """Load the CRSS training set and compute SSH descriptors."""
    mat = scipy.io.loadmat(str(DATA_SIM_DIR / "crss_simulation.mat"))
    df = pd.DataFrame({
        "PHI":  mat["PHI"].flatten(),
        "phi2": mat["phi2"].flatten(),
        "Yind": mat["yind"].flatten() / S_0_SIM,    # normalize -> orientation factor
    })

    gsh = gshUtils("cubic_triclinic", 6)
    ssh = gsh.struct_to_ssh(np.array(df[["PHI", "phi2"]])).real
    df["l4m1"]  = ssh[:, 1]
    df["l12m2"] = ssh[:, -1]
    return df


def plot_ssh_space(df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(df["l4m1"], df["l12m2"], s=40)
    ax.set_title("SSH Coefficient Space (CRSS sims)", fontsize=22)
    ax.invert_yaxis()
    ax.set_ylabel("l=12, m=2", fontsize=18, labelpad=12)
    ax.set_xlabel("l=4, m=1",  fontsize=18, labelpad=12)
    ax.tick_params(axis="both", which="major", labelsize=14)
    plt.show()


# ---------------------------------------------------------------------------
# Step 2 — Train GP surrogate (orientation -> yield factor)
# ---------------------------------------------------------------------------
# Features:  l4m1, l12m2   (2D)
# Target:    Yind          (dimensionless yield factor)
# Kernel:    Matern-5/2 with ARD + White noise
# ---------------------------------------------------------------------------

FEATURE_NAMES = ["l4m1", "l12m2"]
OUTPUT_NAME   = "Yind"


def build_surrogate(
    C_df: pd.DataFrame,
    feature_names=FEATURE_NAMES,
    output_name=OUTPUT_NAME,
    test_size: float = 0.2,
    random_state: int = 12,
    cache_path: Path = None,
    force_retrain: bool = False,
    num_restarts: int = 0,
    fix_lengthscale=None,
    fix_signal_variance: float = None,
    fix_noise_variance: float = None,
):
    """Fit a GPy Matern-5/2+White surrogate on standardized features/target.

    Parameters
    ----------
    num_restarts : int
        If > 0, runs `model.optimize_restarts(num_restarts)` instead of
        single-start optimize. Helps avoid local minima on small datasets.
    fix_lengthscale : list-like of length n_features, optional
        If given, pins the Matern-5/2 lengthscales to these values and
        skips optimizing them. The published CRSS surrogate converged to
        (7.937, 87.588) — passing those exactly reproduces it.
    fix_signal_variance / fix_noise_variance : float, optional
        Pin the kernel signal variance and/or noise variance to these
        values. Often paired with fix_lengthscale for full reproducibility.
    """
    from types import SimpleNamespace

    if cache_path is not None and not force_retrain:
        cache_path = Path(cache_path)
        if cache_path.exists() and cache_path.stat().st_size > 0:
            try:
                with open(cache_path, "rb") as f:
                    surrogate = pickle.load(f)
                print(f"[surrogate] loaded from cache: {cache_path}")
                return surrogate
            except (EOFError, pickle.UnpicklingError):
                print(f"[surrogate] cache at {cache_path} is corrupt; rebuilding")

    X = C_df[feature_names]
    y = C_df[output_name]

    X_train, X_test, y_train, y_test, train_ind, test_ind = train_test_split(
        X, y, np.arange(X.shape[0]),
        test_size=test_size, shuffle=True, random_state=random_state,
    )

    y_train = np.array(y_train).reshape(-1, 1)
    y_test  = np.array(y_test).reshape(-1, 1)
    X_train = np.array(X_train)
    X_test  = np.array(X_test)

    scaler_in  = StandardScaler().fit(X_train)
    scaler_out = StandardScaler().fit(y_train)

    X_train_s = scaler_in.transform(X_train)
    y_train_s = scaler_out.transform(y_train)
    X_test_s  = scaler_in.transform(X_test)
    y_test_s  = scaler_out.transform(y_test)

    n_features = len(feature_names)
    kernel = GPy.kern.Matern52(n_features, ARD=True) + GPy.kern.White(n_features)
    model  = GPy.models.GPRegression(X_train_s, y_train_s, kernel)

    # Optionally pin specific hyperparameters before optimization
    if fix_lengthscale is not None:
        ls = np.array(fix_lengthscale, dtype=float)
        if ls.shape != (n_features,):
            raise ValueError(f"fix_lengthscale must be length {n_features}, got {ls.shape}")
        model.sum.Mat52.lengthscale = ls
        model.sum.Mat52.lengthscale.fix()
        print(f"  [surrogate] pinned lengthscales = {ls.tolist()}")
    if fix_signal_variance is not None:
        model.sum.Mat52.variance = float(fix_signal_variance)
        model.sum.Mat52.variance.fix()
        print(f"  [surrogate] pinned signal variance = {fix_signal_variance}")
    if fix_noise_variance is not None:
        model.sum.white.variance = float(fix_noise_variance)
        model.sum.white.variance.fix()
        model.Gaussian_noise.variance = float(fix_noise_variance)
        model.Gaussian_noise.variance.fix()
        print(f"  [surrogate] pinned noise variance = {fix_noise_variance}")

    if num_restarts > 0:
        print(f"  [surrogate] multi-restart optimization ({num_restarts} restarts)...")
        model.optimize_restarts(num_restarts=num_restarts, messages=False, verbose=False)
    else:
        model.optimize(messages=True)

    surrogate = SimpleNamespace(
        model=model,
        kernel=kernel,
        scaler_in=scaler_in,
        scaler_out=scaler_out,
        X_train_scaled=X_train_s,
        X_test_scaled=X_test_s,
        y_train_scaled=y_train_s,
        y_test_scaled=y_test_s,
        train_ind=train_ind,
        test_ind=test_ind,
        feature_names=feature_names,
    )

    if cache_path is not None:
        cache_path = Path(cache_path)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "wb") as f:
            pickle.dump(surrogate, f)
        print(f"[surrogate] saved to cache: {cache_path}")
    return surrogate


def plot_parity(surrogate) -> None:
    m = surrogate.model
    scaler_out = surrogate.scaler_out
    y_pred_train_mean, _ = m.predict(surrogate.X_train_scaled)
    y_pred_test_mean,  _ = m.predict(surrogate.X_test_scaled)
    plt.title("CRSS Surrogate Parity Plot")
    plt.xlabel("Predicted Yind")
    plt.ylabel("True Yind")
    plt.scatter(scaler_out.inverse_transform(y_pred_train_mean),
                scaler_out.inverse_transform(surrogate.y_train_scaled),
                c="b", label="train")
    plt.scatter(scaler_out.inverse_transform(y_pred_test_mean),
                scaler_out.inverse_transform(surrogate.y_test_scaled),
                c="r", label="test")
    plt.legend()
    plt.show()


# ---------------------------------------------------------------------------
# Step 3 — Likelihood for CRSS
# ---------------------------------------------------------------------------
# At each MCMC step ParaDRAM proposes a single scalar `crss`. The model
# prediction for each experimental grain is then:
#     pred_Yind_i = crss * GP_factor(orientation_i)
# Likelihood is Gaussian in (measured - predicted) with fixed variance.
#
# The squared-error term omits the conventional 1/2 factor, matching the
# Eind likelihood convention used in this codebase.
# ---------------------------------------------------------------------------

LOG_REJECT = -1e300


class CRSSLikelihood:
    """Callable log-likelihood for the 1D CRSS calibration problem."""

    def __init__(self, surrogate, exp_df, fixed_var: float = 0.5):
        self.surrogate = surrogate
        self.fixed_var = fixed_var

        x_exp = np.array(exp_df[["l4m1", "l12m2"]])
        self.y_exp_unscaled = np.array(exp_df["Yind"]).reshape(-1, 1)

        # GP prediction for the experimental orientations only depends on
        # the orientations (which are fixed), NOT on the sampled CRSS.
        # So we precompute it once at construction.
        X_exp_scaled = surrogate.scaler_in.transform(x_exp)
        mean_scaled, _ = surrogate.model.predict(X_exp_scaled)
        self.mean_unscaled = surrogate.scaler_out.inverse_transform(mean_scaled)

    def __call__(self, crss, plot: bool = False) -> float:
        """Log-likelihood at scalar CRSS. ParaDRAM passes a length-1 array."""
        c = float(crss[0]) if hasattr(crss, "__len__") else float(crss)
        if c <= 0:
            return LOG_REJECT

        predicted = self.mean_unscaled * c
        residual = self.y_exp_unscaled - predicted
        var = self.fixed_var

        if plot:
            plt.scatter(self.y_exp_unscaled, predicted)
            lo = float(np.min(self.y_exp_unscaled))
            hi = float(np.max(self.y_exp_unscaled))
            plt.plot([lo, hi], [lo, hi], "k--")
            plt.xlabel("Measured Yind")
            plt.ylabel(f"Predicted Yind  (CRSS = {c:.4f})")
            plt.title(f"RMSE = {float(np.sqrt(np.mean(residual**2))):.3f}")
            plt.axis("equal")
            plt.show()

        return float(np.sum(
            np.log(1) - 0.5 * np.log(2 * np.pi * var) - residual ** 2 / var
        ))


# ---------------------------------------------------------------------------
# Step 4 — Run ParaDRAM (1D)
# ---------------------------------------------------------------------------

def pick_random_seed() -> int:
    return random.SystemRandom().randint(1, 2 ** 31 - 1)


def _clean_paradram_outputs(prefix: Path) -> None:
    prefix = Path(prefix)
    for f in prefix.parent.glob(prefix.name + "_process_*"):
        try:
            f.unlink()
        except OSError as e:
            print(f"  warning: could not delete {f}: {e}")


def run_paradram(
    likelihood,
    output_prefix: Path = PARADRAM_OUT_DIR / "gaussian",
    chain_size: int = 12_000,
    target_acceptance_rate: float = 0.15,
    random_seed=None,
    start_point=(0.3,),
    domain_lower=(0.0,),
    domain_upper=(1.0,),
    overwrite: bool = True,
) -> int:
    """Configure and run ParaDRAM for the 1D CRSS calibration."""
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    if random_seed is None:
        random_seed = pick_random_seed()
        print(f"  random_seed was None -> generated {random_seed} "
              f"(pass --seed {random_seed} to reproduce)")

    pmpd = pm.ParaDRAM()
    pmpd.spec.targetAcceptanceRate = target_acceptance_rate
    pmpd.spec.chainSize            = chain_size
    pmpd.spec.randomSeed           = random_seed
    pmpd.spec.overwriteRequested   = overwrite
    pmpd.spec.outputFileName       = str(output_prefix)
    pmpd.spec.domainLowerLimitVec  = list(domain_lower)
    pmpd.spec.domainUpperLimitVec  = list(domain_upper)
    pmpd.spec.variableNameList     = ["CRSS"]
    pmpd.spec.startPointVec        = list(start_point)

    print(f"  CRSS prior: uniform on [{domain_lower[0]}, {domain_upper[0]}]")
    pmpd.runSampler(ndim=1, getLogFunc=likelihood)

    # Sidecar JSON
    seed_path = output_prefix.parent / (output_prefix.name + ".seed.json")
    try:
        seed_path.write_text(json.dumps({
            "random_seed": int(random_seed),
            "chain_size":  int(chain_size),
            "target_acceptance_rate": float(target_acceptance_rate),
            "start_point": list(start_point),
            "domain_lower": list(domain_lower),
            "domain_upper": list(domain_upper),
        }, indent=2))
    except OSError as e:
        print(f"  (could not write seed sidecar at {seed_path}: {e})")

    return int(random_seed)


# ---------------------------------------------------------------------------
# Step 5 — Post-process the chain
# ---------------------------------------------------------------------------

def load_paradram_sample(
    output_prefix: Path = PARADRAM_OUT_DIR / "gaussian",
) -> pd.DataFrame:
    pmpd = pm.ParaDRAM()
    pmpd.readSample(str(output_prefix))
    return pmpd.sampleList[0].df


def summarize_posterior(sample_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for v in ["CRSS"]:
        rows.append({
            "param": v,
            "mean": float(np.mean(sample_df[v])),
            "std":  float(np.sqrt(np.var(sample_df[v]))),
            "mean_MPa": float(np.mean(sample_df[v])) * 1000,
            "std_MPa":  float(np.sqrt(np.var(sample_df[v]))) * 1000,
            "n": len(sample_df),
        })
    return pd.DataFrame(rows)


def save_sample_csv(sample_df: pd.DataFrame, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    sample_df.to_csv(path, index=False)


# ---------------------------------------------------------------------------
# Driver: single-chain pipeline
# ---------------------------------------------------------------------------

def run_single_chain(seed, chain_size: int, prefix: Path,
                     show_plots: bool = False,
                     domain_lower=(0.0,), domain_upper=(1.0,),
                     start_point=(0.3,),
                     num_restarts: int = 0,
                     fix_lengthscale=None,
                     fix_signal_variance=None,
                     fix_noise_variance=None,
                     exp_csv: str = "ar_final.csv"):
    # Step 1 — load training set
    C_df = load_crss_simulation_data()
    print(f"[chain seed={seed}] simulation set: {len(C_df)} orientations")

    # Step 2 — surrogate (cached). Cache key includes any hyperparameter
    # pinning so a `--fix-lengthscale` run doesn't reuse a cache trained
    # without pinning (and vice versa).
    cache_for_run = SURROGATE_CACHE
    if fix_lengthscale is not None:
        cache_for_run = CACHE_DIR / "surrogate_crss_pinned.pkl"
    elif num_restarts > 0:
        cache_for_run = CACHE_DIR / "surrogate_crss_restarts.pkl"
    surrogate = build_surrogate(
        C_df, cache_path=cache_for_run,
        num_restarts=num_restarts,
        fix_lengthscale=fix_lengthscale,
        fix_signal_variance=fix_signal_variance,
        fix_noise_variance=fix_noise_variance,
    )
    if show_plots:
        plot_ssh_space(C_df)
        plot_parity(surrogate)

    # Step 3 — experimental likelihood
    exp_csv_path = DATA_EXP_DIR / exp_csv
    exp_df = pd.read_csv(exp_csv_path)
    print(f"[chain seed={seed}] experimental data: {len(exp_df)} rows from {exp_csv_path.name}")
    likelihood = CRSSLikelihood(surrogate, exp_df, fixed_var=0.5)
    log_p_pub = likelihood(0.314)   # near published mean
    print(f"[chain seed={seed}] log L at CRSS=0.314 = {log_p_pub:.4f}")

    if show_plots:
        likelihood(0.314, plot=True)

    # Step 4 — clean + sample
    _clean_paradram_outputs(prefix)
    seed_label = "auto" if seed is None else str(seed)
    print(f"[chain seed={seed_label}] sampling {chain_size} "
          f"(writing to {prefix.name}_*)...")
    seed_used = run_paradram(
        likelihood=likelihood,
        output_prefix=prefix,
        chain_size=chain_size,
        target_acceptance_rate=0.15,
        random_seed=seed,
        start_point=start_point,
        domain_lower=domain_lower,
        domain_upper=domain_upper,
    )

    sample_df = load_paradram_sample(prefix)
    print(f"[chain seed={seed_used}] refined sample: n={len(sample_df)}, "
          f"max LogFunc = {sample_df['SampleLogFunc'].max():.3f}")
    return sample_df, seed_used


# ---------------------------------------------------------------------------
# Driver: multi-chain orchestrator
# ---------------------------------------------------------------------------

def run_multichain(n_chains: int, base_seed, chain_size: int,
                   output_dir: Path = PARADRAM_OUT_DIR,
                   num_restarts: int = 0,
                   fix_lengthscale_str: str = None,
                   fix_signal_variance: float = None,
                   fix_noise_variance: float = None,
                   exp_csv: str = "ar_final.csv"):
    output_dir.mkdir(parents=True, exist_ok=True)
    if base_seed is None:
        seeds = [pick_random_seed() for _ in range(n_chains)]
    else:
        seeds = [base_seed + i for i in range(n_chains)]

    procs, prefixes, logs = [], [], []
    for i, seed in enumerate(seeds):
        prefix = output_dir / f"chain_seed{seed}"
        _clean_paradram_outputs(prefix)
        log_path = output_dir / f"chain_seed{seed}.log"

        cmd = [sys.executable, str(Path(__file__).resolve()),
               "--worker",
               "--seed", str(seed),
               "--chain-size", str(chain_size),
               "--prefix", str(prefix)]
        if num_restarts > 0:
            cmd += ["--num-restarts", str(num_restarts)]
        if fix_lengthscale_str:
            cmd += ["--fix-lengthscale", fix_lengthscale_str]
        if fix_signal_variance is not None:
            cmd += ["--fix-signal-variance", str(fix_signal_variance)]
        if fix_noise_variance is not None:
            cmd += ["--fix-noise-variance", str(fix_noise_variance)]
        if exp_csv and exp_csv != "ar_final.csv":
            cmd += ["--exp-csv", exp_csv]
        f_log = open(log_path, "w")
        logs.append(f_log)
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                             stdout=f_log, stderr=subprocess.STDOUT)
        p.stdin.write(b"n\n"); p.stdin.flush(); p.stdin.close()
        procs.append(p)
        prefixes.append(prefix)
        print(f"  spawned chain {i+1}/{n_chains}  seed={seed}  "
              f"log={log_path.name}")

    print(f"\nWaiting for {n_chains} chains to finish "
          f"(check {output_dir}/*.log for progress)...")
    for i, p in enumerate(procs):
        rc = p.wait()
        logs[i].close()
        if rc != 0:
            print(f"  chain {i+1} FAILED (return code {rc}) -- see "
                  f"{output_dir / f'chain_seed{seeds[i]}.log'}")
        else:
            print(f"  chain {i+1}/{n_chains} done  (seed={seeds[i]})")

    dfs = []
    for i, pref in enumerate(prefixes):
        try:
            df = load_paradram_sample(pref)
            df["chain"] = seeds[i]
            dfs.append(df)
        except Exception as e:
            print(f"  could not load chain seed={seeds[i]}: {e}")
    if not dfs:
        raise RuntimeError("All chains failed.")
    return pd.concat(dfs, ignore_index=True), seeds


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser():
    p = argparse.ArgumentParser(
        description="Bayesian calibration of CRSS via GP surrogate + ParaDRAM.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed (base seed if --n-chains > 1). "
                        "If omitted, a fresh seed from OS entropy is used.")
    p.add_argument("--chain-size", type=int, default=12_000,
                   help="ParaDRAM chainSize per chain.")
    p.add_argument("--n-chains", type=int, default=1,
                   help="Number of parallel chains (subprocesses).")
    p.add_argument("--prefix", type=str, default=None,
                   help="ParaDRAM output prefix.")
    p.add_argument("--show-plots", action="store_true",
                   help="Show matplotlib pop-ups for diagnostic plots.")
    p.add_argument("--num-restarts", type=int, default=0,
                   help="Run multi-restart GP hyperparameter optimization "
                        "with this many random starts. 0 = single start.")
    p.add_argument("--fix-lengthscale", type=str, default=None,
                   help="Pin GP lengthscales to these values "
                        "(comma-separated, e.g. '7.937,87.588' to reproduce "
                        "the published surrogate exactly).")
    p.add_argument("--fix-signal-variance", type=float, default=None,
                   help="Pin GP Matern signal variance (e.g. 15.253).")
    p.add_argument("--fix-noise-variance", type=float, default=None,
                   help="Pin GP white-noise variance (e.g. 0.00908).")
    p.add_argument("--exp-csv", type=str, default="ar_final.csv",
                   help="Experimental CSV filename within data/experimental/. "
                        "ar_final.csv is the locked final dataset (v1 rows "
                        "with v2 Yind values where available, 32 total rows).")
    p.add_argument("--lower", type=float, default=0.0,
                   help="CRSS domain lower bound.")
    p.add_argument("--upper", type=float, default=1.0,
                   help="CRSS domain upper bound.")
    p.add_argument("--start", type=float, default=0.3,
                   help="CRSS starting value for ParaDRAM.")
    p.add_argument("--inspect-only", action="store_true",
                   help="Build the surrogate, print train/test split + "
                        "fitted lengthscales, show parity plot, and exit. "
                        "No sampling is performed.")
    p.add_argument("--worker", action="store_true",
                   help="(internal) flag set when spawned as a child of --n-chains.")
    return p


def inspect_surrogate(num_restarts: int = 0,
                      fix_lengthscale=None,
                      fix_signal_variance=None,
                      fix_noise_variance=None):
    """Build the surrogate and print its train/test split plus parity diagnostics."""
    C_df = load_crss_simulation_data()
    print(f"\n=== Simulation orientation set ({len(C_df)} rows) ===")
    cols = [c for c in ["l4m1", "l12m2", "Yind"] if c in C_df.columns]
    print(C_df[cols].to_string())

    cache_for_run = SURROGATE_CACHE
    if fix_lengthscale is not None:
        cache_for_run = CACHE_DIR / "surrogate_crss_pinned.pkl"
    elif num_restarts > 0:
        cache_for_run = CACHE_DIR / "surrogate_crss_restarts.pkl"
    surrogate = build_surrogate(
        C_df, cache_path=cache_for_run,
        num_restarts=num_restarts,
        fix_lengthscale=fix_lengthscale,
        fix_signal_variance=fix_signal_variance,
        fix_noise_variance=fix_noise_variance,
    )

    print(f"\n=== Train/test split (random_state=12, test_size=0.2) ===")
    print(f"  Train indices into C_df: {sorted(surrogate.train_ind.tolist())}")
    print(f"  Test indices into C_df:  {sorted(surrogate.test_ind.tolist())}")

    print(f"\n=== Training points (n={len(surrogate.train_ind)}) ===")
    train_df = C_df.iloc[surrogate.train_ind][["l4m1", "l12m2", "Yind"]].copy()
    train_df["original_index"] = surrogate.train_ind
    print(train_df.to_string())

    print(f"\n=== Test points (n={len(surrogate.test_ind)}) ===")
    test_df = C_df.iloc[surrogate.test_ind][["l4m1", "l12m2", "Yind"]].copy()
    test_df["original_index"] = surrogate.test_ind
    print(test_df.to_string())

    print(f"\n=== Trained GP model ===")
    print(surrogate.model)
    print(f"\nMatern-5/2 ARD lengthscales:")
    for name, ls in zip(surrogate.feature_names,
                        surrogate.kernel.parameters[0].lengthscale):
        print(f"  {name:6s}: {float(ls):.6f}")

    # Predictions on train/test
    y_pred_train_mean, _ = surrogate.model.predict(surrogate.X_train_scaled)
    y_pred_test_mean,  _ = surrogate.model.predict(surrogate.X_test_scaled)
    train_true = surrogate.scaler_out.inverse_transform(surrogate.y_train_scaled)
    train_pred = surrogate.scaler_out.inverse_transform(y_pred_train_mean)
    test_true  = surrogate.scaler_out.inverse_transform(surrogate.y_test_scaled)
    test_pred  = surrogate.scaler_out.inverse_transform(y_pred_test_mean)

    print(f"\n=== Predictions vs truth (train) ===")
    for i, (t, p) in enumerate(zip(train_true.flatten(), train_pred.flatten())):
        print(f"  C_df[{int(surrogate.train_ind[i]):2d}]: true={t:.5f}  pred={p:.5f}  err={t-p:+.5f}")
    train_rmse = float(np.sqrt(np.mean((train_true - train_pred) ** 2)))
    print(f"  RMSE(train) = {train_rmse:.6f}")

    print(f"\n=== Predictions vs truth (test) ===")
    for i, (t, p) in enumerate(zip(test_true.flatten(), test_pred.flatten())):
        print(f"  C_df[{int(surrogate.test_ind[i]):2d}]: true={t:.5f}  pred={p:.5f}  err={t-p:+.5f}")
    test_rmse = float(np.sqrt(np.mean((test_true - test_pred) ** 2)))
    print(f"  RMSE(test)  = {test_rmse:.6f}")

    plot_ssh_space(C_df)
    plot_parity(surrogate)


def _parse_lengthscale(s):
    if s is None:
        return None
    parts = [float(x.strip()) for x in s.split(",")]
    return parts


if __name__ == "__main__":
    args = _build_parser().parse_args()
    fix_lengthscale = _parse_lengthscale(args.fix_lengthscale)

    if args.inspect_only:
        inspect_surrogate(
            num_restarts=args.num_restarts,
            fix_lengthscale=fix_lengthscale,
            fix_signal_variance=args.fix_signal_variance,
            fix_noise_variance=args.fix_noise_variance,
        )
        sys.exit(0)

    # ---- multi-chain master ----
    if args.n_chains > 1 and not args.worker:
        if args.seed is None:
            print(f"=== Multi-chain CRSS run: {args.n_chains} chains, "
                  f"seeds auto-generated, chainSize={args.chain_size} each ===")
        else:
            print(f"=== Multi-chain CRSS run: {args.n_chains} chains, "
                  f"seeds {args.seed}..{args.seed + args.n_chains - 1}, "
                  f"chainSize={args.chain_size} each ===")
        combined, seeds_used = run_multichain(
            n_chains=args.n_chains,
            base_seed=args.seed,
            chain_size=args.chain_size,
            num_restarts=args.num_restarts,
            fix_lengthscale_str=args.fix_lengthscale,
            fix_signal_variance=args.fix_signal_variance,
            fix_noise_variance=args.fix_noise_variance,
            exp_csv=args.exp_csv,
        )
        print("\nSeeds used (pass --seed N to reproduce):")
        for s in seeds_used:
            print(f"  {s}")

        print("\n=== Per-chain summaries ===")
        for seed, grp in combined.groupby("chain"):
            m = grp["CRSS"].mean(); s = np.sqrt(np.var(grp["CRSS"]))
            print(f"  seed={seed}  n={len(grp)}  "
                  f"CRSS = {m:.4f} +/- {s:.4f}  "
                  f"({m*1000:.2f} +/- {s*1000:.2f} MPa)")

        print("\n=== Pooled posterior across all chains (Stage 3) ===")
        print(summarize_posterior(combined).to_string(index=False))

        out_csv = RESULTS_DIR / "y_ind_sampler_ar_multichain.csv"
        save_sample_csv(combined, out_csv)
        print(f"\nSaved pooled sample to: {out_csv}")
        sys.exit(0)

    # ---- single chain (default, or worker mode) ----
    prefix = Path(args.prefix) if args.prefix else PARADRAM_OUT_DIR / "gaussian"
    sample_df, seed_used = run_single_chain(
        seed=args.seed,
        chain_size=args.chain_size,
        prefix=prefix,
        show_plots=args.show_plots and not args.worker,
        domain_lower=(args.lower,),
        domain_upper=(args.upper,),
        num_restarts=args.num_restarts,
        fix_lengthscale=fix_lengthscale,
        fix_signal_variance=args.fix_signal_variance,
        fix_noise_variance=args.fix_noise_variance,
        exp_csv=args.exp_csv,
        start_point=(args.start,),
    )

    if not args.worker:
        print(f"\nSeed used: {seed_used}   (rerun with --seed {seed_used} to reproduce)")
        print("\n=== Posterior summary ===")
        print(summarize_posterior(sample_df).to_string(index=False))
        out_csv = RESULTS_DIR / "y_ind_sampler_ar_reproduced.csv"
        save_sample_csv(sample_df, out_csv)
        print(f"\nSaved sample to: {out_csv}")
