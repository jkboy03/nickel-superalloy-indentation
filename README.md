# Bayesian calibration of elastic constants and CRSS

Reproduction code for the Bayesian calibration described in the paper.

## Setup

```
pip install -r requirements.txt
python -c "import paramonte"          # answer "n" to the MPI prompt first time
```

Tested with Python 3.10.

## Run the calibrations

```
python scripts/01_indentation_calibration.py --n-chains 4
python scripts/02_crss_calibration.py        --n-chains 4
python scripts/03_analyze_results.py
```

The first command runs the elastic-constants MCMC (~3 minutes, 4 cores),
the second runs the CRSS MCMC, and the third writes a Markdown summary
and seven publication PNG figures into `figures/`. All outputs land in
`results/` and `figures/`, both of which are recreated automatically.

## MATLAB IPDF plots

`scripts/plot_ipdf.m` produces the inverse-pole-figure plots. Requires
MTEX. Edit the `project_root` variable at the top of that file if you
move the project on disk.
