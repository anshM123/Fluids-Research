# Code

Run every command from the package root (`ncs-submission/`). The setup is:

```bash
python3 -m pip install -r requirements.txt      # Python 3.11; NumPy, SciPy, Numba, Matplotlib
```

No GPU is needed. All run times below are for one core of a four-core x86-64 machine.

## Layout

| Folder | Content |
|---|---|
| `lib/` | the solver library: the marching Newton–Krylov solver for IPM, its log-polar base and regridding, the independent global sparse Newton solver, the local eigenproblem and the resonance phase, and linear stability. The 2D Boussinesq modules (`bq_*`) are the base classes of the IPM modules. |
| `compute/` | the production runs: branch continuation, profile location, scans, grid shifts, origin truncation, deep continuation, the global solver, and eigenvalue counts and spectra. Also the job scripts of the seventh-profile campaign. |
| `analysis/` | scripts that turn computed states and logs into the reported quantities: phase tables, geometry, the seventh-profile analysis and the registered predictions. `make_si_tables.py` writes the Supplementary Tables. |
| `figures/` | `make_figures.py` writes every figure of the paper. |
| `tests/` | smoke tests. |

Every script in `compute/` and `analysis/` adds `lib/` to its import path. It can be run from any directory and
writes its outputs into the current directory.

## 1. Smoke tests (15 s)

```bash
python3 -m pytest code/tests -q            # or: python3 code/tests/test_smoke.py
```

The tests check four things:
- every solver module imports;
- the stored profiles λ₁ and λ₃ solve the equations, with residual below 10⁻⁹ and m − 2 below 10⁻⁸;
- the exact Jacobian of the global solver matches central differences;
- the Supplementary Tables and two figures are regenerated from `data/outputs`, with the tables identical byte
  for byte.

## 2. Figures and Supplementary Tables from the logged outputs (1 min)

```bash
python3 code/figures/make_figures.py          # figures/Fig1–6 and SupplementaryFig1–3, as PDF and PNG
python3 code/analysis/make_si_tables.py       # supplementary/tables/*.tex
```

Both read only `data/outputs`. The environment variables `NCS_DATA`, `NCS_FIGS` and `NCS_TABLES` redirect input and
output.

## 3. Analyses from the logged outputs (seconds each)

These scripts recompute the seventh-profile determination, the registered predictions and the endpoint analysis from
the logs. Their output is identical to the logged file:

```bash
cd data/outputs
python3 ../../code/analysis/ipm_lambda7_final.py     | diff - ipm_lambda7_final.out
python3 ../../code/analysis/ipm_lambda7_analysis.py  | diff - ipm_lambda7_analysis_committed_rule.out
python3 ../../code/analysis/ipm_stage3_predict.py    | diff - ipm_stage3_predictions.out
python3 ../../code/analysis/ipm_stage4_predict.py    | diff - ipm_stage4_predictions.out
python3 ../../code/analysis/ipm_endpoint_geometry.py | diff - ipm_endpoint_geometry.out
cd ../..
```

## 4. Recomputation from the stored profiles (minutes)

Run these in a scratch directory. `P` is the absolute path of `ncs-submission`.

**Resonance phase** (50 s per state). This reproduces the logged rows of `ipm_phase3_rungs.out` and
`ipm_phase3_l7.out` to all printed digits, for example Φ₀ = 11.4568 − 5.4154i at λ₃ and 23.9580 − 13.6584i at
z = 7.346:

```bash
python3 $P/code/analysis/ipm_phase3_table.py $P/data/profiles/ipm_profile_n3_lam0.2415663353_Nb32_hs0.0125.npy \
                                             $P/data/profiles/ipm_state_z7.3460_Nb32_hs0.00625.npy
```

**Independent global solver** (15 s per s_max at h_s = 0.025, about 1 min at 0.0125). This reproduces
`ipm_glob_l2.log`, giving λ = 0.3149632866 at s_max = 6 and 0.3149615196 at s_max = 8:

```bash
python3 $P/code/compute/ipm_glob_run.py $P/data/profiles/ipm_profile_n2_lam0.3149618108_Nb32_hs0.0125.npy \
        0.3149618108 l2 0.025 24 -12 6,8,10,12 0.0125
```

The production settings of Supplementary Table 4 are h_s = 0.0125, N_b = 24, s_min = −16 and s_max = 8, 10, 12. For
λ₄ they also include h_s = 0.00625.

**Profile location** (secant iteration on m(λ) = 2). At h_s = 0.025 this takes 1.5 min and recovers
λ₀ = 1.0285722975:

```bash
F=$P/data/profiles/ipm_profile_n0_lam1.0285722975_Nb32_hs0.025.npy
python3 $P/code/compute/ipm_crossing.py $F 1.0290 $F 1.0280 32 0.025 -20 100 32 0.025
```

One solve takes 3–8 min at h_s = 0.0125 and 11–23 min at 0.00625. The arguments are
`Y_a lam_a Y_b lam_b Nb hs s_start s_max Nb_in hs_in`.

## 5. The full computational campaign

The campaign that produced `data/outputs` ran in this order. Each production log records its arguments in its
first line (`# ipm_scan …`, `# ipm_deep …`, `# IPM global solver …`).

| Step | Script | Output in `data/outputs` | Cost |
|---|---|---|---|
| 1. Continue the branch from λ = 1 downwards (h_s = 0.025) | `compute/ipm_branch.py LAM_START LAM_END DZ TAG START.npy` | `ipm_branch_dn*` | about 2 core-hours |
| 2. Locate λ₀–λ₆ at h_s = 0.025, then 0.0125 | `compute/ipm_crossing.py` | `ipm_rung*.out` | about 1 core-hour at h_s = 0.0125 |
| 3. Scan past λ₆ at fixed z | `compute/ipm_scan.py START.npy LAM Z0 Z1 DZ SS TAG [hs]` | `ipm_scan_{s20,s16,p7,h7}` | about 3 core-hours |
| 4. Seventh profile: grid shifts and systematics | `compute/ipm_shift.py STATE Z DELTA [hs]` via `rich7_run.sh` and `rich7_run2.sh` (queues `rich7_queue*.txt`) | `rich7_*.out`, `sys7_*.out` | about 4.5 core-hours |
| 5. Origin truncation and noise floor | `compute/ipm_sstart.py`, `compute/ipm_noise2.py` | `ipm_sstart_*`, `ipm_noise2_*`, `ipm_nfloor_*` | minutes per solve |
| 6. Deep continuation (h_s = 0.0125) | `compute/ipm_deep.py LAM_START LAM_END DZ TAG START.npy [Nb hs NB_IN HS_IN tol]` | `ipm_deep_e5.out` | about 3 core-hours |
| 7. Independent global solver | `compute/ipm_glob_run.py` | `ipm_glob_*.log` | minutes per profile |
| 8. Resonance phase, tracker diagnostics, geometry | `analysis/ipm_phase3_table.py`, `ipm_tracker_diag.py`, `ipm_travel_time.py`, `ipm_phase_anatomy3.py`, `ipm_wall_profiles.py`, `ipm_K_checks.py`, `ipm_phase_variants.py`, `ipm_cutoff_drift.py`, `ipm_deepres_check.py` | `ipm_phase3_*`, `ipm_tracker_diag.*`, `ipm_travel_time.out`, … | 30–60 s per state |
| 9. Seventh-profile analysis and registered predictions | `analysis/ipm_lambda7_final.py`, `ipm_lambda7_analysis.py`, `ipm_stage3_predict.py`, `ipm_stage4_predict.py`, `ipm_endpoint_geometry.py` | see section 3 | seconds |
| 10. Instability index and spectra (λ₀–λ₄) | `compute/ipm_contour.py`, `compute/ipm_flow_spec.py` | `ipm_contour_U*`, `ipm_spec_U*` | minutes per contour |

Some environment variables control the runs:
- `NK_TOL` sets the Newton tolerance;
- `NB`, `SS` and `SMAX` set N_b, s_start and s_max for `ipm_shift.py` and `ipm_scan.py`;
- `IPM_SSTART` sets s_start for `ipm_deep.py`;
- `SAVE_EVERY=1` saves the state at every scan point.

Steps 1–6 save intermediate states (`*.npy`) that are not part of the package. Only the ten states in
`data/profiles` are included. The job scripts in `compute/` are records of the campaign as it was run in the
directory `programs/P03_boussinesq_ladder` of the full repository, and they refer to those intermediate states. Every
job resumes from its last saved state after an interruption.

The solvers for the 2D Boussinesq and Hou–Luo validation (Supplementary Note 1) are in the same directory of the
full repository.
