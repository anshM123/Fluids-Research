# Data

Two folders:

| Folder | Content | Size |
|---|---|---|
| `profiles/` | ten converged IPM states in double precision: the smooth profiles λ₀–λ₆ and three deep states | 31 MB |
| `outputs/` | every logged output from which a figure or table of the paper is generated (133 files) | 0.7 MB |

The data are released under CC BY 4.0 (see `LICENSE-DATA.md`).

## `profiles/`: converged states

Each file holds the unknowns of the marching solver, `X = Ψ/r²` (Ψ is the stream function), on the log-polar
grid. The array has one row per radial point `s_k = −120 + k·h_s` with `s_k ≥ −20` (the origin truncation
`s_start`) and one column per angular point `β_j = (π/4)(1 − cos(πj/32))`, `j = 0…32`. Column `j = 0` is the wall
and column `j = 32` is the symmetry axis.

| File | n | λ | z = 1/λ | h_s | rows × columns | Newton residual ǀRǀ | m − 2 |
|---|---|---|---|---|---|---|---|
| `ipm_profile_n0_lam1.0285722975_Nb32_hs0.025.npy` | 0 | 1.0285722975 | 0.9722 | 0.025 | 4800 × 33 | 3.1e-11 | +3.2e-11 |
| `ipm_profile_n1_lam0.4721297348_Nb32_hs0.0125.npy` | 1 | 0.4721297348 | 2.1181 | 0.0125 | 9600 × 33 | 6.9e-11 | −2.0e-10 |
| `ipm_profile_n2_lam0.3149618108_Nb32_hs0.0125.npy` | 2 | 0.3149618108 | 3.1750 | 0.0125 | 9600 × 33 | 1.3e-11 | −4.1e-11 |
| `ipm_profile_n3_lam0.2415663353_Nb32_hs0.0125.npy` | 3 | 0.2415663353 | 4.1396 | 0.0125 | 9600 × 33 | 8.3e-11 | +2.8e-10 |
| `ipm_profile_n4_lam0.1987224523_Nb32_hs0.0125.npy` | 4 | 0.1987224523 | 5.0321 | 0.0125 | 9600 × 33 | 1.3e-11 | +3.1e-11 |
| `ipm_profile_n5_lam0.1706180880_Nb32_hs0.0125.npy` | 5 | 0.1706180880 | 5.8610 | 0.0125 | 9600 × 33 | 4.7e-11 | +1.8e-10 |
| `ipm_profile_n6_lam0.15092_Nb32_hs0.0125.npy` | 6 | 0.15092 | 6.6260 | 0.0125 | 9600 × 33 | 4.4e-13 | +4.8e-10 |
| `ipm_state_z7.3460_Nb32_hs0.00625.npy` | (7) | 1/7.346 | 7.3460 | 0.00625 | 19200 × 33 | 6.6e-13 | +3.0e-11 |
| `ipm_state_z8.0660_Nb32_hs0.00625.npy` | — | 1/8.066 | 8.0660 | 0.00625 | 19200 × 33 | 7.3e-13 | +1.0e-9 |
| `ipm_state_z8.7860_Nb32_hs0.00625.npy` | — | 1/8.786 | 8.7860 | 0.00625 | 19200 × 33 | 9.1e-13 | +2.0e-9 |

The residual and m − 2 columns were recomputed from these files with the released solver (`code/tests`).

- For λ₀–λ₅ the residual is set by rounding λ to ten digits in the file name; the solver tolerance was
  7 × 10⁻¹³.
- λ₆ = 0.15092 is the profile located to the defect floor (main text).
- The state at z = 7.346 is the production state of the seventh-profile audit (Extended Data Table 1).
- The states at z = 8.066 and 8.786 are the resolution checks of the deep branch (Supplementary Table 12).

### Loading a state

```python
import sys, numpy as np
sys.path.insert(0, "code/lib")
from ipm_solver import IPM
from bq_newton import full, residual

lam, hs = 0.2415663353, 0.0125
Y = np.load("data/profiles/ipm_profile_n3_lam0.2415663353_Nb32_hs0.0125.npy")
B = IPM(lam, hs=hs, Nb=32, s_sw=12.0, s_start=-20.0)    # the grid of the file
R, info = residual(B, Y)                                 # max|R| ~ 1e-10, info["m"] - 2 ~ 3e-10
X = full(B, Y)                                           # X = Ψ/r² on the whole grid s ∈ [−120, 100)
r = B.march(X / B.ea2[:, None], return_all=True)         # fields: r["Th"] = R̂, r["Omega"], r["Ur"], r["m"], r["A"]
```

`B.s` and `B.beta` are the grid coordinates. The density is `R = cos^m(β)·R̂`. The wall speed is
`D = (1 + λ) + r["Ur"][:, 0]`.

## `outputs/`: logged outputs

Every figure (`code/figures/make_figures.py`) and every Supplementary Table (`code/analysis/make_si_tables.py`) is
generated from these files alone. Each `.out` file is the log written by the script named in the third column. The
`.npy` and `.npz` files are the arrays those scripts saved.

| Files | Content | Written by | Used in |
|---|---|---|---|
| `ipm_branch_dn*.{out,npy}` | coarse continuation of the branch, λ from 1 downwards (h_s = 0.025): m(λ), A(λ), wall dip | `compute/ipm_branch.py` | Fig. 2a |
| `ipm_rung{1..6}_hs0125.out` | location of the smooth profiles λ₁–λ₆ (secant iteration on m = 2), h_s = 0.0125 | `compute/ipm_crossing.py` | Figs 1, 2, 5, Supplementary Tables 4, 8 |
| `ipm_rung{0..6}.out`, `ipm_rung6_hs00625.out` | the same at h_s = 0.025 (λ₀ and the coarse values) and λ₆ at h_s = 0.00625 | `compute/ipm_crossing.py` | `ladders.json`, Supplementary Note 2 |
| `ipm_rung{4,5,6}_Nb48.out`, `ipm_rung5_hs0125_smax130.out`, `ipm_rung5_hs0125_dom12_40.out` | profile location with N_b = 48, s_max = 130 and a shortened domain | `compute/ipm_crossing.py` | Supplementary Note 2 |
| `ipm_scan_{s20,p7,h7}.{out,npy}` | fixed-z scans of m − 2 from the sixth profile onwards (`s20`, h_s = 0.0125) and around the seventh (`p7`, h_s = 0.0125; `h7`, h_s = 0.00625) | `compute/ipm_scan.py` | Figs 2, 3, Supplementary Table 2 |
| `ipm_scan_s16.{out,npy}` | the `s20` scan repeated with origin truncation s_start = −16 | `compute/ipm_scan.py` | Supplementary Note 2 |
| `rich7_*_d*_hs*.out` (32 files) | m − 2 on radial grids shifted by δh_s, δ = 0, ¼, ½, ¾, at four z and three h_s | `compute/ipm_shift.py` via `compute/rich7_run*.sh` | Fig. 3a, Supplementary Fig. 1, Supplementary Table 2 |
| `sys7_{nb48,nb64,ss22,ss24,sm130,sm160,combined}.out` | systematic variants of the production state at z = 7.346 (N_b, s_start, s_max) | `compute/ipm_shift.py` (environment variables `NB`, `SS`, `SMAX`) | Extended Data Table 1, `ipm_lambda7_final.out` |
| `ipm_lambda7_final.out` | seventh-profile determination from the converged grid with all systematics | `analysis/ipm_lambda7_final.py` | Figs 1–4, Extended Data Table 2 |
| `ipm_lambda7_analysis_committed_rule.out` | the committed λ₇ decision rule, applied as written | `analysis/ipm_lambda7_analysis.py` | Supplementary Note 3 |
| `ipm_sstart_lam{5,6}.out`, `ipm_noise2_ss{16,20}.out`, `ipm_nfloor_ss{26,30}_lam6.out` | origin truncation, Newton floor and m-noise at fixed λ | `compute/ipm_sstart.py`, `compute/ipm_noise2.py` | Fig. 3c,d, Supplementary Table 3 |
| `ipm_nbcheck*_lam{5,6}.out`, `ipm_noise_lam6_hs0125.out` | Newton–Krylov logs of single solves at λ₅, λ₆: angular resolution N_b = 32, 48, 64 and Newton tolerance | `lib/bq_newton.py` | Supplementary Note 2 |
| `ipm_glob_l{1..6}{,f,h}.log` | independent global solver: λ_n for s_max = 8, 10, 12 and extrapolation (no suffix: h_s = 0.025; `f`: 0.0125; `h`: 0.00625) | `compute/ipm_glob_run.py` | Fig. 2c, Supplementary Table 4 |
| `ipm_phase3_{rungs,l7,deep}.out` | repaired resonance phase Φ₀ on every computed state | `analysis/ipm_phase3_table.py` | Figs 1, 2, 5, 6, Supplementary Tables 8–10 |
| `ipm_wkb_e5_cut2.out`, `ipm_wkb_fine_deep_cut2.out` | phase from the original tracker on the deep states (the failure of Supplementary Note 5) | `lib/ipm_wkb.py` | Fig. 6, Supplementary Table 10 |
| `ipm_tracker_diag.{out,npz}` | K = κD̂ along the wall, original against repaired tracker | `analysis/ipm_tracker_diag.py` | Fig. 6 |
| `ipm_K_checks.{out,npz}` | exact scaling κ = K/D̂ and Re K against the wall gradient | `analysis/ipm_K_checks.py` | Fig. 5c,d |
| `ipm_phase_variants.out` | phase and defect on every audited variant at the seventh profile | `analysis/ipm_phase_variants.py` | Figs 3, 5b, Extended Data Table 1 |
| `ipm_cutoff_drift.out` | offsets δ_n for front cut-offs D/D₀ = 1.5, 2, 3, 5 | `analysis/ipm_cutoff_drift.py` | Supplementary Table 7 |
| `ipm_phase_anatomy3.out` | partition of the phase into inner layer, approach, dip core and front side | `analysis/ipm_phase_anatomy3.py` | Fig. 1b, Supplementary Fig. 2a |
| `ipm_wall_profiles.npz` | wall speed along the wall for six states and the wavenumber density at z = 7.346 | `analysis/ipm_wall_profiles.py` | Fig. 1a,b |
| `ipm_deep_e5.out` | deep continuation (h_s = 0.0125) from z = 7.466 to 10.226: m − 2, dip depth, steepest gradient, resolution indicator | `compute/ipm_deep.py` | Fig. 6, Supplementary Fig. 2, Supplementary Table 11 |
| `ipm_travel_time.out` | travel time T = ∫ds/D̂ and dip and front positions | `analysis/ipm_travel_time.py` | Supplementary Fig. 2, Supplementary Table 11 |
| `ipm_deepres_z{8.066,8.786}.out`, `deepres_z*_hs0.00625.out` | deep states re-solved on h_s = 0.00625 | `analysis/ipm_deepres_check.py` | Fig. 6, Supplementary Table 12 |
| `ipm_endpoint_geometry.out` | endpoint fits and termination indicators | `analysis/ipm_endpoint_geometry.py` | Supplementary Note 6 |
| `ipm_stage3_predictions.out`, `ipm_stage4_predictions.out` | registered predictions (stages 3 and 4), as computed at registration | `analysis/ipm_stage{3,4}_predict.py` | Fig. 4 (stage 3, values as registered), Supplementary Table 6 (stage 4) |
| `ipm_contour_U{0..4}.out`, `ipm_spec_U{0..4}.out` | right-half-plane eigenvalue counts and unstable spectra of λ₀–λ₄ | `compute/ipm_contour.py`, `compute/ipm_flow_spec.py` | Supplementary Note 1 |
| `ladders.json` | profile positions of the 2D Boussinesq, Hou–Luo and IPM hierarchies, with their sources | assembled from the logs | Supplementary Fig. 3 |

The 2D Boussinesq and Hou–Luo validation data (Supplementary Note 1) are in the full repository under
`programs/P03_boussinesq_ladder/`.
