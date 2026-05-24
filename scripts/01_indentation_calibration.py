"""
Bayesian calibration of FCC cubic elastic constants (C11, C12, C44)
from nanoindentation modulus measurements using a Gaussian Process
surrogate and ParaDRAM MCMC.
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
from glob import glob

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
DATA_SIM_DIR = ROOT / "data" / "simulation_eind"
DATA_EXP_DIR = ROOT / "data" / "experimental"
RESULTS_DIR  = ROOT / "results"
CACHE_DIR    = ROOT / "cache"
SURROGATE_CACHE = CACHE_DIR / "surrogate_eind.pkl"


# ---------------------------------------------------------------------------
# Step 1 — Load FCC simulation dataset
# ---------------------------------------------------------------------------
# eind_simulation.mat holds the 539-row training set used by the published
# calibration. Fields (as the FEM simulation produced them):
#   C11, C12, C44   — cubic elastic constants (GPa)
#   PHI, phi2       — Euler angles of the indented crystal (radians)
#   Eind            — indentation modulus (GPa, target of the GP surrogate)
# SSH descriptors (l4m1, l12m2) are computed at load time via gshutils.
# ---------------------------------------------------------------------------

def load_simulation_data(sim_dir: Path = DATA_SIM_DIR) -> pd.DataFrame:
    """Load the FCC simulation training set and compute SSH descriptors."""
    mat = scipy.io.loadmat(str(sim_dir / "eind_simulation.mat"))
    df = pd.DataFrame({
        "C11":  mat["C11"].flatten(),
        "C12":  mat["C12"].flatten(),
        "C44":  mat["C44"].flatten(),
        "PHI":  mat["PHI"].flatten(),
        "phi2": mat["phi2"].flatten(),
        "Eind": mat["Eind"].flatten(),
    })

    # Convert Euler angles -> SSH coefficients in the cubic-triclinic basis.
    # The first two non-trivial coefficients (l=4,m=1) and (l=12,m=2) are
    # used as orientation descriptors by the surrogate.
    gsh = gshUtils("cubic_triclinic", 6)
    ssh = gsh.struct_to_ssh(np.array(df[["PHI", "phi2"]])).real
    df["l4m1"]  = ssh[:, 1]
    df["l12m2"] = ssh[:, -1]
    return df


def plot_simulation_space(df: pd.DataFrame) -> None:
    """3D scatter of the (C11, C12, C44) sampling, colored by Eind."""
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection="3d")

    im = ax.scatter(df["C11"], df["C12"], df["C44"],
                    c=df["Eind"], marker="o", s=30, edgecolors="black")

    ax.set_xlabel("C11", fontsize=18, labelpad=12)
    ax.set_ylabel("C12", fontsize=18, labelpad=12)
    ax.set_zlabel("C44", fontsize=18, labelpad=12)
    fig.colorbar(im, label="Eind (GPa)")
    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.view_init(30, 210)
    plt.show()


def plot_ssh_space(df: pd.DataFrame) -> None:
    """2D scatter of the simulation set in SSH coefficient space."""
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(df["l4m1"], df["l12m2"], s=40)
    ax.set_title("SSH Coefficient Space", fontsize=24)
    ax.invert_yaxis()
    ax.set_ylabel("l=12, m=2", fontsize=18, labelpad=12)
    ax.set_xlabel("l=4, m=1",  fontsize=18, labelpad=12)
    ax.tick_params(axis="both", which="major", labelsize=14)
    plt.show()


# ---------------------------------------------------------------------------
# Step 3 — Train Gaussian Process surrogate (forward model)
# ---------------------------------------------------------------------------
# Features:  C11, C12, C44, l4m1, l12m2  (5D)
# Target:    Eind
# Kernel:    Matern-5/2 with ARD + White noise
# The surrogate is the forward map (elastic constants + orientation) -> Eind
# that the MCMC likelihood will call thousands of times.
# ---------------------------------------------------------------------------

FEATURE_NAMES = ["C11", "C12", "C44", "l4m1", "l12m2"]
OUTPUT_NAME   = "Eind"


def build_surrogate(
    C_df: pd.DataFrame,
    feature_names=FEATURE_NAMES,
    output_name=OUTPUT_NAME,
    test_size: float = 0.2,
    random_state: int = 1,
    cache_path: Path = None,
    force_retrain: bool = False,
):
    """Fit a GPy Matern-5/2+White surrogate on standardized features/target.

    If `cache_path` is given and the file exists, load the cached surrogate
    instead of retraining. Pass `force_retrain=True` to ignore the cache.
    To invalidate the cache, delete the .pkl file (or set force_retrain=True).

    Returns a SimpleNamespace bundling the trained model, the two scalers,
    and the train/test arrays needed for the parity plot.
    """
    from types import SimpleNamespace

    # --- Cache hit: just load and return ---
    if cache_path is not None and not force_retrain:
        cache_path = Path(cache_path)
        if cache_path.exists():
            with open(cache_path, "rb") as f:
                surrogate = pickle.load(f)
            print(f"[surrogate] loaded from cache: {cache_path}")
            return surrogate

    # --- Cache miss (or forced): train from scratch ---
    X = C_df[feature_names]
    y = C_df[output_name]

    X_train, X_test, y_train, y_test, train_ind, test_ind = train_test_split(
        X, y, np.arange(X.shape[0]),
        test_size=test_size,
        shuffle=True,
        random_state=random_state,
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

    # --- Save to cache for next run ---
    if cache_path is not None:
        cache_path = Path(cache_path)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "wb") as f:
            pickle.dump(surrogate, f)
        print(f"[surrogate] saved to cache: {cache_path}")

    return surrogate


def plot_parity(surrogate) -> None:
    """Parity plot (predicted vs. true) in the original Eind units (GPa)."""
    m = surrogate.model
    scaler_out = surrogate.scaler_out

    y_pred_train_mean, _ = m.predict(surrogate.X_train_scaled)
    y_pred_test_mean,  _ = m.predict(surrogate.X_test_scaled)

    plt.title("Parity Plot")
    plt.xlabel("Predicted")
    plt.ylabel("True")
    plt.scatter(scaler_out.inverse_transform(y_pred_train_mean),
                scaler_out.inverse_transform(surrogate.y_train_scaled), c="b", label="train")
    plt.scatter(scaler_out.inverse_transform(y_pred_test_mean),
                scaler_out.inverse_transform(surrogate.y_test_scaled), c="r", label="test")
    plt.legend()
    plt.show()


# ---------------------------------------------------------------------------
# Step 4 — Experimental data + likelihood
# ---------------------------------------------------------------------------
# The likelihood is the function ParaDRAM will call thousands of times.
# It evaluates: for a proposed (C11, C12, C44) — assuming every experimental
# grain shares the same elastic constants but has its own orientation
# (l4m1, l12m2) — what is the log-probability that the surrogate's predicted
# Eind matches the measured Eind, under a Gaussian noise model with fixed
# variance? Two unphysical regions are also rejected outright.
#-----------------------------------

# A pathological-likelihood sentinel used by ParaDRAM to reject a sample.
LOG_REJECT = -1e300


class IndentationLikelihood:
    """Callable log-likelihood for the indentation calibration problem.

    Parameters
    ----------
    surrogate : SimpleNamespace
        Output of build_surrogate(). Provides model + scalers.
    exp_df : pd.DataFrame
        Experimental data with columns 'l4m1', 'l12m2', 'Eind'.
    fixed_var : float
        Hard-coded Gaussian likelihood variance (in scaled Eind units).
        This replaces the GP's own predictive variance with a fixed value
        so the likelihood width is controlled directly.

    Usage
    -----
        like = IndentationLikelihood(surrogate, exp_df)
        log_prob = like([260.0, 160.0, 130.0])      # for ParaDRAM
        df       = like.inverse_transform([...])    # diagnostic
    """

    def __init__(self, surrogate, exp_df, fixed_var: float = 0.05,
                 prior_mean=None, prior_std=None):
        # fixed_var = 0.05 reproduces the published max SampleLogFunc of
        # 18.52: with 32 experimental points,
        #   N * (-0.5 * log(2*pi * 0.05)) = 32 * 0.579 = 18.5.
        self.surrogate = surrogate
        self.fixed_var = fixed_var

        # Optional Gaussian prior on (C11, C12, C44). If both prior_mean and
        # prior_std are given, log_prior is added to the data log-likelihood,
        # so ParaDRAM samples from the proper posterior under this prior.
        if prior_mean is not None and prior_std is not None:
            self.prior_mean = np.array(prior_mean, dtype=float)
            self.prior_std  = np.array(prior_std,  dtype=float)
            if self.prior_mean.shape != (3,) or self.prior_std.shape != (3,):
                raise ValueError("prior_mean / prior_std must each be length 3")
            if (self.prior_std <= 0).any():
                raise ValueError("prior_std must be positive in every component")
        else:
            self.prior_mean = None
            self.prior_std  = None

        x_exp = np.array(exp_df[["l4m1", "l12m2"]])
        y_exp = np.array(exp_df["Eind"]).reshape(-1, 1)

        # Standardized target the GP works in
        self.Eind = surrogate.scaler_out.transform(y_exp)

        # Pre-allocated 5-column feature matrix: cols 0-2 are filled per
        # likelihood call (proposed C11, C12, C44); cols 3-4 are the fixed
        # per-grain orientation descriptors.
        self._proposed = np.zeros((x_exp.shape[0], 5))
        self._proposed[:, -2:] = x_exp

    # -- internal helpers ---------------------------------------------------

    def _violates_bounds(self, C) -> bool:
        # Physical lower bounds for a stable cubic crystal:
        #   C11 > C12   and   C11 + 2*C12 > 0
        if C[0] - C[1] <= 0:
            return True
        if C[0] + 2 * C[1] <= 0:
            return True
        if C[0] < C[1]:
            return True
        if C[1] < C[2]:
            return True
        return False

    def _predict(self, C):
        """Update self._proposed[:, :3] with C, then call the GP."""
        self._proposed[:, 0] = C[0]
        self._proposed[:, 1] = C[1]
        self._proposed[:, 2] = C[2]
        prop_scaled = self.surrogate.scaler_in.transform(self._proposed)
        mean, _ = self.surrogate.model.predict(prop_scaled)
        return mean

    # -- public API ---------------------------------------------------------

    def __call__(self, C, plot: bool = False) -> float:
        """Log-likelihood evaluated at C = (C11, C12, C44). This is what ParaDRAM calls."""
        if self._violates_bounds(C):
            return LOG_REJECT

        mean = self._predict(C)
        var = self.fixed_var

        if plot:
            plt.scatter(self.Eind, mean)
            plt.plot([np.min(self.Eind), np.max(self.Eind)],
                     [np.min(self.Eind), np.max(self.Eind)])
            plt.show()

        # Data log-likelihood.
        log_L_data = float(np.sum(
            np.log(1) - 0.5 * np.log(2 * np.pi * var)
            - (self.Eind.reshape(-1, 1) - mean) ** 2 / var
        ))

        # Optional Gaussian prior contribution (zero if no prior given).
        if self.prior_mean is not None:
            C_arr = np.asarray(C, dtype=float)
            log_prior = float(np.sum(
                -0.5 * np.log(2 * np.pi * self.prior_std ** 2)
                - 0.5 * ((C_arr - self.prior_mean) / self.prior_std) ** 2
            ))
            return log_L_data + log_prior

        return log_L_data

    def inverse_transform(self, C, plot: bool = False):
        """Diagnostic: predictions in original GPa units, as a DataFrame.

        Mirrors the original `likelihood_inverse_transform`. Note the original
        returns -1e300 (not a DataFrame) when bounds are violated — preserved.
        """
        if self._violates_bounds(C):
            return LOG_REJECT

        mean = self._predict(C)

        scaler_out = self.surrogate.scaler_out
        mean_orig = scaler_out.inverse_transform(mean)
        Eind_orig = scaler_out.inverse_transform(self.Eind)

        df = pd.DataFrame(self._proposed,
                          columns=["C11", "C12", "C44", "l4m1", "l12m2"])
        df["Eind"]         = Eind_orig
        df["Eind_predict"] = mean_orig
        df["error"]        = Eind_orig - mean_orig

        if plot:
            plt.figure(figsize=(6, 6))
            plt.scatter(Eind_orig, mean_orig)
            lo, hi = np.min(Eind_orig), np.max(Eind_orig)
            plt.plot([lo, hi], [lo, hi], "k--")
            mse = float(np.mean((mean_orig - Eind_orig) ** 2))
            plt.title(f"Mean Sq. Error: {mse}")
            plt.xlabel("Actual")
            plt.ylabel("Predicted")
            plt.axis("equal")
            plt.show()

        return df


# ---------------------------------------------------------------------------
# Step 5 — Run the ParaDRAM sampler
# ---------------------------------------------------------------------------
# ParaDRAM = Delayed-Rejection Adaptive Metropolis. It tunes its proposal
# covariance on the fly so the chain hits the target acceptance rate.
#
# Output files written under <output_prefix>_process_1_*.txt :
#   *_chain.txt       — raw Markov chain (unique states + weight column)
#   *_sample.txt      — refined / decorrelated sample (use this for stats)
#   *_progress.txt    — running diagnostics
#   *_report.txt      — final summary
# ---------------------------------------------------------------------------

PARADRAM_OUT_DIR = RESULTS_DIR / "paradram_out"


def pick_random_seed() -> int:
    """Generate a fresh random seed using OS entropy. Always positive 32-bit
    so paramonte's Fortran integer ingester is happy."""
    return random.SystemRandom().randint(1, 2**31 - 1)


def run_paradram(
    likelihood,
    C_df: pd.DataFrame,
    output_prefix: Path = PARADRAM_OUT_DIR / "gaussian",
    chain_size: int = 12_000,
    target_acceptance_rate: float = 0.15,
    random_seed=None,
    start_point=(300.0, 180.0, 140.0),
    variable_names=("C11", "C12", "C44"),
    overwrite: bool = True,
    lower_bounds=None,
    upper_bounds=None,
) -> int:
    """Configure and run ParaDRAM. Returns the random_seed that was used.

    If `random_seed` is None, a fresh seed is generated with OS entropy
    so each run differs. The seed is always written to the chain alongside
    the output files (so you can reproduce a run later with --seed N).

    The domain bounds come from the min/max of the surrogate's training
    set (C_df), so the sampler can never propose a C outside the region
    where the GP was trained.
    """
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

    # Default bounds come from the simulation training set min/max.
    # Override via `lower_bounds=` / `upper_bounds=` (length-3 tuples) to add
    # physical prior knowledge (e.g. C12 should never be < 100 GPa for Ni).
    default_lower = [
        float(C_df["C11"].min()),
        float(C_df["C12"].min()),
        float(C_df["C44"].min()),
    ]
    default_upper = [
        float(C_df["C11"].max()),
        float(C_df["C12"].max()),
        float(C_df["C44"].max()),
    ]
    pmpd.spec.domainLowerLimitVec = list(lower_bounds) if lower_bounds is not None else default_lower
    pmpd.spec.domainUpperLimitVec = list(upper_bounds) if upper_bounds is not None else default_upper
    print(f"  domain bounds: C11=[{pmpd.spec.domainLowerLimitVec[0]:.1f}, "
          f"{pmpd.spec.domainUpperLimitVec[0]:.1f}]  "
          f"C12=[{pmpd.spec.domainLowerLimitVec[1]:.1f}, "
          f"{pmpd.spec.domainUpperLimitVec[1]:.1f}]  "
          f"C44=[{pmpd.spec.domainLowerLimitVec[2]:.1f}, "
          f"{pmpd.spec.domainUpperLimitVec[2]:.1f}]")
    pmpd.spec.variableNameList = list(variable_names)
    pmpd.spec.startPointVec    = list(start_point)

    pmpd.runSampler(ndim=len(variable_names), getLogFunc=likelihood)

    # Drop a small sidecar JSON next to the outputs recording the seed used
    seed_path = output_prefix.parent / (output_prefix.name + ".seed.json")
    try:
        seed_path.write_text(json.dumps({
            "random_seed": int(random_seed),
            "chain_size":  int(chain_size),
            "target_acceptance_rate": float(target_acceptance_rate),
            "start_point": list(start_point),
        }, indent=2))
    except OSError as e:
        print(f"  (could not write seed sidecar at {seed_path}: {e})")

    return int(random_seed)


# ---------------------------------------------------------------------------
# Step 6 — Post-process the chain
# ---------------------------------------------------------------------------
# ParaMonte's Python API is annoying because the read methods (readSample,
# readChain, ...) populate hidden list attributes (sampleList, chainList, ...)
# rather than returning anything directly. We wrap that here.
# ---------------------------------------------------------------------------

def load_paradram_sample(
    output_prefix: Path = PARADRAM_OUT_DIR / "gaussian",
) -> pd.DataFrame:
    """Load the FINAL refined sample (= Stage 3) produced by a prior run.

    Equivalent to:
        pmpd = pm.ParaDRAM()
        pmpd.readSample(<prefix>)
        return pmpd.sampleList[0].df
    """
    pmpd = pm.ParaDRAM()
    pmpd.readSample(str(output_prefix))
    return pmpd.sampleList[0].df


def _stats(df, variables):
    """Plain (unweighted) mean/std for each variable. np.var is biased (ddof=0)."""
    return {v: (float(df[v].mean()),
                float(np.sqrt(np.var(df[v]))))
            for v in variables}


def _weighted_stats(df, variables, weight_col="SampleWeight"):
    """Frequency-weighted mean/std — what you want on a compact chain."""
    w = df[weight_col].to_numpy()
    W = w.sum()
    out = {}
    for v in variables:
        x = df[v].to_numpy()
        m = (x * w).sum() / W
        s2 = ((x - m) ** 2 * w).sum() / W
        out[v] = (float(m), float(np.sqrt(s2)))
    return out


def _try_match(label, vals, target=(292.0, 193.0, 140.0), tol=5.0):
    """Tag a (C11, C12, C44) mean triple as MATCH if all within tol of target."""
    if any(v is None for v in vals):
        return ""
    diffs = [abs(v - t) for v, t in zip(vals, target)]
    return "  <-- MATCHES PUBLISHED" if all(d <= tol for d in diffs) else ""


def _dump_df(label, df, variables):
    """Print n, cols, and unweighted (mean, std) for a DataFrame."""
    if df is None:
        print(f"  {label}: <None>")
        return None
    if not isinstance(df, pd.DataFrame):
        print(f"  {label}: not a DataFrame ({type(df).__name__})")
        return None
    n = len(df)
    cols = df.columns.tolist()
    present = [v for v in variables if v in cols]
    means = [float(df[v].mean()) if v in cols else None for v in variables]
    stds  = [float(np.sqrt(np.var(df[v]))) if v in cols else None for v in variables]
    tag = _try_match(label, means)
    print(f"  {label}: n={n}  cols={cols}{tag}")
    for v, m, s in zip(variables, means, stds):
        if m is not None:
            print(f"      {v}: mean={m:8.3f}  std={s:7.3f}")
    # Weighted (if SampleWeight present)
    if "SampleWeight" in cols:
        w = df["SampleWeight"].to_numpy()
        W = w.sum()
        wmeans = []
        for v in variables:
            if v in cols:
                x = df[v].to_numpy()
                m = (x * w).sum() / W
                s2 = ((x - m) ** 2 * w).sum() / W
                wmeans.append(float(m))
                print(f"      [weighted] {v}: mean={m:8.3f}  std={np.sqrt(s2):7.3f}")
            else:
                wmeans.append(None)
        tag2 = _try_match(label + " [weighted]", wmeans)
        if tag2:
            print(f"      {tag2.strip()}")
    return df


def _explore_object(label, obj, variables, depth=0, seen=None):
    """Walk obj's public attributes, dumping any DataFrames or arrays found."""
    if seen is None:
        seen = set()
    if id(obj) in seen or depth > 2:
        return
    seen.add(id(obj))

    indent = "    " * depth
    attrs = [a for a in dir(obj) if not a.startswith("_")]
    for a in attrs:
        try:
            val = getattr(obj, a)
        except Exception:
            continue
        if callable(val):
            continue
        if isinstance(val, pd.DataFrame):
            print(f"{indent}.{a}  (DataFrame)")
            _dump_df(f"{label}.{a}", val, variables)
        elif isinstance(val, np.ndarray):
            shape = val.shape
            if val.size > 0 and val.size <= 50:
                print(f"{indent}.{a}  (ndarray shape={shape}) = {val.flatten()[:12]}")
            else:
                print(f"{indent}.{a}  (ndarray shape={shape})")
        elif hasattr(val, "__dict__") and not isinstance(val, (str, int, float, bool, list, tuple, dict)):
            # Nested ParaMonte object — recurse one level
            sub_attrs = [x for x in dir(val) if not x.startswith("_")]
            if any(isinstance(getattr(val, x, None), (pd.DataFrame, np.ndarray)) for x in sub_attrs):
                print(f"{indent}.{a}  ({type(val).__name__})")
                _explore_object(f"{label}.{a}", val, variables, depth + 1, seen)


def diagnose_paradram_outputs(
    output_prefix: Path = PARADRAM_OUT_DIR / "gaussian",
    variables=("C11", "C12", "C44"),
) -> None:
    """Exhaustively try every ParaMonte read method and inspect every
    DataFrame/stats object exposed on each loaded list element.

    We don't know which extraction method produced the published 2238-row
    `E_ind_sampler_ar.csv`. This walks every public attribute on every
    loaded object so the matching one becomes obvious.
    """
    output_prefix = Path(output_prefix)
    print("=" * 78)
    print(f"Exhaustive diagnosis of ParaDRAM outputs at: {output_prefix}")
    print(f"Target (published) means: C11=292  C12=193  C44=140")
    print("=" * 78)

    # All known reader/list pairs
    readers = [
        ("readSample",      "sampleList"),
        ("readChain",       "chainList"),
        ("readMarkovChain", "markovChainList"),
        ("readProgress",    "progressList"),
        ("readReport",      "reportList"),
        ("readRestart",     "restartList"),
    ]

    for method_name, list_name in readers:
        print(f"\n{'-' * 78}")
        print(f"### pmpd.{method_name}(...) -> pmpd.{list_name}")
        print("-" * 78)
        try:
            pmpd = pm.ParaDRAM()
            getattr(pmpd, method_name)(str(output_prefix))
            objs = getattr(pmpd, list_name, None)
            if not objs:
                print(f"  pmpd.{list_name} is empty")
                continue
            for idx, obj in enumerate(objs):
                print(f"\n  --- {list_name}[{idx}]  type={type(obj).__name__} ---")
                # First: the canonical .df
                df = getattr(obj, "df", None)
                if isinstance(df, pd.DataFrame):
                    _dump_df(f"{list_name}[{idx}].df", df, variables)
                # Then walk every public attribute
                _explore_object(f"{list_name}[{idx}]", obj, variables, depth=1)
        except Exception as e:
            print(f"  {method_name} FAILED: {e!r}")

    # Raw file reads — last resort
    print(f"\n{'-' * 78}")
    print("### Raw pd.read_csv on output files")
    print("-" * 78)
    for suffix in [
        "_process_1_chain.txt",
        "_process_1_sample.txt",
        "_process_1_markovChain.txt",
        "_process_1_progress.txt",
        "_process_1_report.txt",
    ]:
        f = Path(str(output_prefix) + suffix)
        if not f.exists():
            print(f"  (missing)  {f.name}")
            continue
        try:
            df = pd.read_csv(f)
            _dump_df(f"pd.read_csv({f.name})", df, variables)
        except Exception as e:
            print(f"  read_csv({f.name}) FAILED: {e!r}")

    print("\n" + "=" * 78)
    print("Look above for the line tagged `<-- MATCHES PUBLISHED`.")
    print("Published target: C11=292+/-24  C12=193+/-34  C44=140+/-12")
    print("=" * 78)


def summarize_posterior(sample_df: pd.DataFrame,
                        variable_names=("C11", "C12", "C44")) -> pd.DataFrame:
    """Mean and std for each elastic constant + the anisotropy A."""
    df = sample_df.copy()
    df["A"] = 2.0 * df["C44"] / (df["C11"] - df["C12"])

    rows = []
    for v in list(variable_names) + ["A"]:
        rows.append({
            "param": v,
            "mean":  float(np.mean(df[v])),
            "std":   float(np.sqrt(np.var(df[v]))),  # matches original (biased)
            "n":     len(df),
        })
    return pd.DataFrame(rows)


def save_sample_csv(sample_df: pd.DataFrame, path: Path) -> None:
    """Save the posterior sample (including A) to CSV."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out = sample_df.copy()
    if "A" not in out.columns:
        out["A"] = 2.0 * out["C44"] / (out["C11"] - out["C12"])
    out.to_csv(path, index=False)


# ---------------------------------------------------------------------------
# Driver: single-chain pipeline
# ---------------------------------------------------------------------------

def _clean_paradram_outputs(prefix: Path) -> None:
    """Delete stale ParaDRAM output files for `prefix` (paramonte 2.5.2's
    overwriteRequested flag is unreliable; we wipe them manually)."""
    prefix = Path(prefix)
    for f in prefix.parent.glob(prefix.name + "_process_*"):
        try:
            f.unlink()
        except OSError as e:
            print(f"  warning: could not delete {f}: {e}")


def run_single_chain(seed, chain_size: int, prefix: Path,
                     show_plots: bool = False,
                     lower_bounds=None, upper_bounds=None,
                     prior_mean=None, prior_std=None):
    """Run the full pipeline once for a given seed and output prefix.

    If `seed` is None, a fresh OS-entropy seed is generated. Returns
    `(sample_df, seed_used)`.
    """
    # Step 1 — load simulation data
    C_df = load_simulation_data()
    print(f"[chain seed={seed}] simulation set: {len(C_df)} rows")

    # Step 2 — surrogate (cached)
    surrogate = build_surrogate(C_df, cache_path=SURROGATE_CACHE)
    if show_plots:
        plot_simulation_space(C_df)
        plot_ssh_space(C_df)
        plot_parity(surrogate)

    # Step 3 — experimental likelihood
    exp_df = pd.read_csv(DATA_EXP_DIR / "ar_final.csv")
    likelihood = IndentationLikelihood(
        surrogate, exp_df, fixed_var=0.05,
        prior_mean=prior_mean, prior_std=prior_std,
    )
    if prior_mean is not None:
        print(f"[chain seed={seed}] Gaussian prior on (C11,C12,C44): "
              f"mean={list(prior_mean)}  std={list(prior_std)}")
    log_prob_at_published = likelihood([292.0, 193.0, 140.0])
    print(f"[chain seed={seed}] log posterior at (292, 193, 140) = "
          f"{log_prob_at_published:.4f}")

    if show_plots:
        likelihood.inverse_transform([292.0, 193.0, 140.0], plot=True)

    # Step 4 — clean + sample
    _clean_paradram_outputs(prefix)
    seed_label = "auto" if seed is None else str(seed)
    print(f"[chain seed={seed_label}] sampling {chain_size} "
          f"(writing to {prefix.name}_*)...")
    seed_used = run_paradram(
        likelihood=likelihood,
        C_df=C_df,
        output_prefix=prefix,
        chain_size=chain_size,
        target_acceptance_rate=0.15,
        random_seed=seed,
        start_point=(300.0, 180.0, 140.0),
        lower_bounds=lower_bounds,
        upper_bounds=upper_bounds,
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
                   bounds_str: str = None,
                   prior_mean_str: str = None,
                   prior_std_str: str = None):
    """Spawn `n_chains` worker processes in parallel."""
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
        if bounds_str:
            cmd += ["--bounds", bounds_str]
        if prior_mean_str and prior_std_str:
            cmd += ["--prior-mean", prior_mean_str, "--prior-std", prior_std_str]

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
        description="Bayesian calibration of C11/C12/C44 via GP surrogate + ParaDRAM.",
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
    p.add_argument("--bounds", type=str, default=None,
                   help="Physical prior bounds as 6 floats: "
                        "'C11lo,C11hi,C12lo,C12hi,C44lo,C44hi'. "
                        "If omitted, uses simulation training-set min/max.")
    p.add_argument("--prior-mean", type=str, default=None,
                   help="Gaussian prior mean as 3 floats 'C11,C12,C44'.")
    p.add_argument("--prior-std", type=str, default=None,
                   help="Gaussian prior std as 3 floats 'sC11,sC12,sC44'.")
    p.add_argument("--worker", action="store_true",
                   help="(internal) flag for child processes of --n-chains.")
    return p


def _parse_bounds(s):
    if s is None:
        return None, None
    parts = [float(x) for x in s.split(",")]
    if len(parts) != 6:
        raise SystemExit(f"--bounds needs 6 floats (got {len(parts)})")
    lower = [parts[0], parts[2], parts[4]]
    upper = [parts[1], parts[3], parts[5]]
    for v, lo, hi in zip(("C11","C12","C44"), lower, upper):
        if lo >= hi:
            raise SystemExit(f"--bounds {v}: lower ({lo}) must be < upper ({hi})")
    return lower, upper


def _parse_prior(mean_str, std_str):
    if mean_str is None and std_str is None:
        return None, None
    if (mean_str is None) != (std_str is None):
        raise SystemExit("--prior-mean and --prior-std must be provided together")
    m = [float(x) for x in mean_str.split(",")]
    s = [float(x) for x in std_str.split(",")]
    if len(m) != 3 or len(s) != 3:
        raise SystemExit("--prior-mean / --prior-std each need 3 floats")
    if any(si <= 0 for si in s):
        raise SystemExit("--prior-std components must be positive")
    return m, s


if __name__ == "__main__":
    args = _build_parser().parse_args()
    lower_bounds, upper_bounds = _parse_bounds(args.bounds)
    prior_mean, prior_std = _parse_prior(args.prior_mean, args.prior_std)

    # ---- multi-chain master ----
    if args.n_chains > 1 and not args.worker:
        if args.seed is None:
            print(f"=== Multi-chain run: {args.n_chains} chains, "
                  f"seeds auto-generated, chainSize={args.chain_size} each ===")
        else:
            print(f"=== Multi-chain run: {args.n_chains} chains, "
                  f"seeds {args.seed}..{args.seed + args.n_chains - 1}, "
                  f"chainSize={args.chain_size} each ===")
        combined, seeds_used = run_multichain(
            n_chains=args.n_chains,
            base_seed=args.seed,
            chain_size=args.chain_size,
            bounds_str=args.bounds,
            prior_mean_str=args.prior_mean,
            prior_std_str=args.prior_std,
        )
        print("\nSeeds used (pass --seed N to reproduce):")
        for s in seeds_used:
            print(f"  {s}")

        print("\n=== Per-chain summaries ===")
        for seed, grp in combined.groupby("chain"):
            print(f"  seed={seed}  n={len(grp)}")
            for v in ["C11", "C12", "C44"]:
                print(f"    {v}: mean={grp[v].mean():7.2f}  "
                      f"std={np.sqrt(grp[v].var(ddof=0)):6.2f}")

        print("\n=== Pooled posterior ===")
        summary = summarize_posterior(combined)
        print(summary.to_string(index=False))

        out_csv = RESULTS_DIR / "E_ind_sampler_ar_multichain.csv"
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
        lower_bounds=lower_bounds,
        upper_bounds=upper_bounds,
        prior_mean=prior_mean,
        prior_std=prior_std,
    )

    if not args.worker:
        print(f"\nSeed used: {seed_used}   (rerun with --seed {seed_used} to reproduce)")
        summary = summarize_posterior(sample_df)
        print("\n=== Posterior summary ===")
        print(summary.to_string(index=False))
        out_csv = RESULTS_DIR / "E_ind_sampler_ar_reproduced.csv"
        save_sample_csv(sample_df, out_csv)
        print(f"\nSaved sample to: {out_csv}")
