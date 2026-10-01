# P03 — The ladder of self-similar blow-up profiles of the 2D Boussinesq equations with boundary

**Question.** DeepMind (Wang et al., arXiv:2509.14185) found with PINNs a stable profile and three unstable
self-similar blow-up profiles (plus a candidate fourth) for 2D Boussinesq with boundary, the standard proxy for
3D axisymmetric Euler with boundary (Hou–Luo scenario). Their blow-up rates follow the empirical law
λ_n ≈ 1 + 1/(1.4187 n + 1.0863). Is this ladder infinite (λ_n → 1), or finite?

**Approach.** A classical, non-neural solver computes the whole one-parameter family of *least-singular*
self-similar profiles. For every λ the family has Θ ≈ −|y₁|^m at the stagnation point, with
m(λ) = (λ−1)/(1+λ−A(λ)) and A = −∂₁U₁(0) the strain there. Smooth profiles are exactly the crossings m(λ) = 2,
which is the 2D analogue of the CCF branch analysis in P02.

## Formulation (see the docstring of `bq_logpolar.py`)
- The self-similar equations Ω + ((1+λ)y+U)·∇Ω = ∂₁Θ and (1−λ)Θ + ((1+λ)y+U)·∇Θ = 0, with −ΔΨ = Ω and
  U = ∇^⊥Ψ, are posed on the quarter plane (symmetry).
- They are solved in log-polar variables s = ln r, β ∈ [0, π/2] (β = 0 is the boundary).
- The exact stagnation-point structure is factored out: Θ = cos^m β Θ̂ and Ω = cos^{m−1} β Ω̂.
- All characteristics leave the origin, so (Θ̂, Ω̂) are marched in s from the exact local solution at s = −20.
- The unknown is X = Ψ/r² on s ≥ −20 (velocity-gradient scale). Newton–Krylov is applied to X − BS(march(X)).
- The Biot–Savart solve uses an FFT in s of e^{(2−a)s}Ω (decaying at both ends) and Chebyshev collocation in β.

## Solver versions
| file | content |
|---|---|
| `bq_logpolar.py` | original solver: FFT-based velocity, explicit RK4 march |
| `bq_solver.py` | **production solver**: 8th-order finite-difference velocity from X, and an implicit 2-stage Gauss–Legendre march for s < s_sw (RK4 beyond) |
| `bq_newton.py` | Newton–Krylov (GMRES). `fd='central'` uses scaled central-difference matvecs, which gives quadratic convergence (1e-7 → 2e-11 → 1.6e-13) |
| `bq_scan.py` / `bq_scan2.py` | λ-continuation (secant predictor, adaptive step) with the original / production solver |
| `bq_crossing.py` | secant iteration on m(λ) = 2 (a smooth profile), with Newton at each λ |
| `bq_regrid.py` | transfer of states between grids (Nb, hs) for resolution studies |

Two numerical pitfalls were found and removed:
1. **Gibbs contamination of the Jacobian.** The FFT derivative of Ψ̃ = e^{(2−a)s}X assumes periodicity in s.
   Newton directions that do not decay at large s grow like e^{17} at s = 100, and the resulting ringing
   corrupted the Jacobian in the uniform-strain direction (residual floor 3.6e-6).
2. **Roundoff-level Jacobian-vector products.** Forward differences with ε = 1e-7 on unit-2-norm Krylov vectors
   perturb each entry by only ~1e-10; the matvecs are then roundoff-limited (floor ~1e-7). Central differences
   with the largest perturbed entry set to 1e-6 restore quadratic convergence.

Stiffness: near the stagnation point the angular transport rate is w/(V_r/r) ≈ 2Aβ/ε with ε = 1+λ−A → 0 as
λ → 1 on smooth profiles. The explicit RK4 march (h = 0.025) is unstable for λ ≲ 1.3 at Nb = 32. The implicit
Gauss–Legendre march is A-stable.

## Results (Nb = 32, hs = 0.025 unless stated)
**The ladder.** Smooth profiles are the zeros of m(λ) − 2 along one branch (`bq_crossing.py`, `cross_l*_nb32.log`):

| n | λ_n | z_n = 1/(λ_n − 1) | notes |
|---|---|---|---|
| 0 | 1.9205593 | 1.08630 | stable (Chen–Hou; Wang et al. 1.9205) |
| 1 | 1.3990961 | 2.50566 | Nb 32/48/64: 1.39909599/1.39909609/1.39909609; hs 0.0125: 1.39909601 |
| 2 | 1.2523487 | 3.96277 | Nb 48 agrees to 9e-8 |
| 3 | 1.1842533 | 5.42733 | |
| 4 | 1.1449857 | 6.89723 | hs 0.0125 (hs 0.025: 1.1449864); candidate 4th unstable in Wang et al. (their line: 1.1479) |
| 5 | 1.1194738 | 8.37004 | new; hs 0.0125 (hs 0.025: 1.1194818) |
| 6 | 1.1015817 | 9.84429 | new; hs 0.0125 (hs 0.025: 1.1015235); Nb 48 changes m by −4.6e-10; hs 0.00625 running |
| 7 | 1.0883384 | 11.32010 | new; hs 0.0125 (hs 0.025: 1.0890078); amplitude of m − 2 here ≈ 1e-7 |

- Spacings z_{n+1} − z_n: 1.41937, 1.45711, 1.46455, 1.46991, 1.47281, 1.47426, 1.47581, increasing monotonically
  towards the asymptotic spacing π/C, which lies between 1.476 and ≈ 1.51 (ASYMPTOTICS.md §5).
- The hs = 0.025 values of λ₆ and λ₇ were off by 6e-5 and 7e-4 because the front sharpens as λ → 1.
  The hs = 0.0125 values are used from λ₄ on.
- Extrema of m − 2: +7.41e-2, −6.80e-3, +7.24e-4, −8.34e-5, +1.0e-5, …
- For λ > λ₀, m decreases monotonically (1.34 at λ = 4.1), so there are no further profiles.
- The branch continues at least to λ = 1.069 (scan D), with no sign of termination.

**Mechanism** (THEORY.md; status of each step in ASYMPTOTICS.md):
- A stalled (quasi-stagnant) boundary region x < x_c ≈ 0.72, where D = V₁/x = O(ε), bounded by a front that becomes
  a square-root cusp as λ → 1.
- The WKB phase Φ = ε⁻¹∫κ ds comes from the local 2D eigenproblem (`bq_local_eig.py`, `bq_wkb2d.py`). It gains
  2.98, 3.03, 3.07, 3.07, 3.115, 3.12 between consecutive profiles, tending to π.
- A fit Re Φ(λ_n) = C z_n + Φ₀ holds to ≤ 0.024 over n = 4–7. The slope depends on the front cut-off:
  π/C = 1.493 (cut D/ε = 2) … 1.4775 (cut 8).
- Extrapolation of the resonance positions depends on the correction exponent:
  - 1.478 with analytic 1/z phase corrections;
  - 1.49–1.50 with the non-analytic correction (∝ z^{0.6}) suggested by the scaling of the dip.
- 3/2 is a candidate, not derived. This is formal asymptotics plus numerics, not a proof.

**Hou–Luo** (`hl_solver.py`, `hl_scan_*.log`, `hl_wkb2.out`): 11 crossings resolved above the noise (z ≤ 13.56), spacings
1.2583, 1.2612, 1.2637, 1.2652, 1.2638 against the WKB value π/(2 Re a₀) ≈ 1.267.
- λ_n^{HL} = 1.99871, 1.44767, 1.28676, 1.21092, 1.16668, 1.13772, 1.11731, 1.10215, 1.09046, …
- Phase gained per profile 3.02–3.09, tending to π.
- Amplitude |m − 2| e^{−Im Φ} ∝ z^{−1.58}.
- The dip closes only as ε → 0 (D̂_min = 0.285 at z = 30.8).

**Figure:** `fig_bq_ladder.png`.

**Independent residual check** (`check_residual.py`: equations evaluated with 6th-order finite differences in s):
- Near the origin (r < 0.14) the relative residual is ~1e-7.
- Across the front region at λ₁ it is 1.6e-4 (Nb 32) and 1.7e-6 (Nb 48).
- The outer region (r ≳ 3) has angular structure that Nb = 32–48 under-resolves (residual 1e-3–1e-2).
- The λ_n themselves are insensitive to this (Nb 32/48 agree to ≤ 1e-7 in λ₁, λ₂ and 5e-10 in m at λ₆). The strain
  at the origin sees only low angular modes, and the boundary transport is resolved to 1e-11.

## Independent reproduction: global Newton solver (`bq_global.py`, `glob_run.py`, `glob_*.log`)
**Every numerical ingredient differs from the marching solver.**
- The unknowns are the unhatted Θ, Ω, X on the (s, β) grid, solved simultaneously (no marching).
- 6th-order upwind-biased finite differences in s; Chebyshev in β (Nb = 24).
- The Biot–Savart law is the local elliptic equation (∂_s+2)²X + X_ββ + Ω = 0, with a Robin far field. There is no
  FFT and no exponential weight.
- Smoothness is imposed through the exact local data at s_min = −12. λ is an eigenvalue fixed by the strain condition.
- Newton with the exact sparse Jacobian and SuperLU.
- Three far-field truncations s_max = 8, 10, 12, extrapolated geometrically.

| n | marching λ_n | global, hs 0.025 | global, hs 0.0125 | global − marching |
|---|---|---|---|---|
| 0 | 1.9205593 | 1.9205610 | — | +1.7e-6 |
| 1 | 1.3990961 | 1.3990960 | — | −1.1e-7 |
| 2 | 1.2523487 | 1.2523481 | — | −5.6e-7 |
| 3 | 1.1842533 | 1.1842512 | 1.1842530 | −3.2e-7 |
| 4 | 1.1449857 | 1.1449777 | 1.1449853 | −4.2e-7 |
| 5 | 1.1194738 | 1.1194516 | 1.1194739 | +1.1e-7 |
| 6 | 1.1015817 | 1.1015204 | 1.1015771; s_min = −20: 1.1015792 (s_max 12), 1.1015821 (s_max 10) | −4.6e-6 … +0.4e-6 |
| 7 | 1.0883384 | (no convergence) | 1.0883480 (s_max = 12 only); s_min = −20: 1.0883386 | +9.6e-6; s_min = −20: +2.1e-7 |

**Notes on the global solver.**
- The truncation at the origin matters: s_min = −8 → −10 → −12 moves λ₁ by −3.2e-6 → −1.7e-7 → −1.1e-7.
- For λ₇, Newton converges only at s_max = 12. The rung is weakly determined: |m − 2| ~ 1e-7 nearby, so the
  Jacobian is nearly singular in the λ direction.
- Both solvers show the same hs trend for the deep rungs (λ₆: +5.8e-5 marching, +5.7e-5 global, from hs 0.025 to
  0.0125).

## Linear stability (`bq_stability.py`, `stab_*.py`, `glob_stab.py`, `gstab_*.log`)
**Setting.** Perturbations are taken ∝ e^{μτ}, with τ = −ln(1−t), and must be regular at the stagnation point.
There are two trivial modes: μ = 1 (time translation) and μ = 0 (the scaling symmetry Θ → L⁻¹Θ(Ly), neutral).

**Method 1 (march-based).**
- μ is an eigenvalue iff the linearised march + Biot–Savart map T_μ has an eigenvalue ν = 1.
- Validation: |T₁v − v|/|v| = 3e-6 on the exact time-translation mode.
- Real-axis scans (`stab_scan.py`) count crossings. `stab_refine.py` solves ν(μ) = 1 by Brent.
- `stab_contour.py` is an argument-principle count of det(I − T_μ) on [0.06, 1.5] × [−1.5, 1.5], so it includes
  complex μ.

**Method 2 (global).**
- The generalised eigenproblem J v + μ M v = 0 on the global discretisation, solved by shift-invert Arnoldi at real
  and complex shifts.
- It reproduces the values of the eigenvalues found by method 1: all 21 of λ₁–λ₆ and the six largest of λ₇, to
  ≤ 5e-5 between the two truncations s_min = −12 and −20.
- Its spectrum also contains discretisation artifacts, none of which appears in method 1:
  - a dense family whose real part depends on the truncation (for λ₁, Re μ ≈ 0.054 at s_min = −12 and ≈ 0.15 at
    −20), with Im μ spaced by about 2π/(travel time through the stalled layer);
  - inflow-boundary modes whose eigenvectors are concentrated at s_min. An example is 0.610 ± 1.104i for λ₂, present
    for both truncations; the argument-principle count of method 1 excludes it;
  - a pseudospectral cloud (residuals ~1e-9) at shifts far from any eigenvalue.
- Method 2 therefore confirms eigenvalue *values*. The counts come from method 1: real scans plus argument-principle
  counts.

**Truncation artifacts in method 1.** The regular r² coefficient of a perturbation carries a factor 1/μ (a resonance
with the neutral scaling mode). A finite origin truncation regularises it. This produces one spurious crossing,
which moves with the truncation: 0.052 (s_start −12), 0.033 (−20), 0.022 (−30). Counts are taken for μ ≥ 0.045.

| n | unstable modes | method 1: brackets (real scan) | method 2: eigenvalues μ_k | method 1 refined |
|---|---|---|---|---|
| 0 | 0 | — | — | |
| 1 | 1 | (0.35, 0.40) | 0.37379 | 0.373789 |
| 2 | 2 | (0.55, 0.60), (0.20, 0.25) | 0.55419, 0.22048 | 0.55418, 0.22048 |
| 3 | 3 | (0.60, 0.65), (0.35, 0.40), (0.15, 0.16) | 0.63437, 0.37541, 0.15492 | |
| 4 | 4 | (0.65, 0.70), (0.45, 0.50), (0.25, 0.30), (0.11, 0.12) | 0.68004, 0.45973, 0.28629, 0.11908 | |
| 5 | 5 | (0.70, 0.75), (0.50, 0.55), (0.35, 0.40), (0.20, 0.25), (0.09, 0.10) | 0.70985, 0.51143, 0.36869, 0.23121, 0.09654 | |
| 6 | 6 | (0.70, 0.75), (0.50, 0.55), (0.40, 0.45), (0.30, 0.35), (0.19, 0.20), (0.08, 0.09) | 0.73139, 0.54481, 0.42577, 0.30749, 0.19415, 0.08154 | |
| 7 | 7 | (0.70, 0.75), (0.55, 0.60), (0.45, 0.50), (0.35, 0.40), (0.25, 0.30), (0.16, 0.17), (0.065, 0.070) | 0.74632, 0.56676, 0.46713, 0.36192, 0.26371, 0.16678; lowest inside the artifact band | 0.0699 (s_start −20 and −30) |

**Complex modes.** Argument-principle counts in [0.06, 1.5] × [−1.5, 1.5] match the real count plus the trivial
mode, so no complex unstable eigenvalues were found in the box:
- λ₁: 1.997 zeros (μ = 1 and 0.374).
- λ₂: 3.000 zeros (μ = 1, 0.554, 0.220). This also rules out the method-2 candidate 0.610 ± 1.104i, which is an
  inflow-boundary artifact.
- λ₃: 4.000 zeros. λ₄: 5.000 zeros.
- λ₅: 5.965 ≈ 6 zeros, with one unresolved phase step at the left edge near 0.06 + 0.58i, next to the
  truncation-artifact region.
- λ₆ and λ₇: not contour-checked.

**Observation (not derived).** The lower unstable eigenvalues scale with ε = (λ_n−1)/2 and are nearly equally
spaced. For n = 6, μ_k/ε = 1.61, 3.82, 6.05, 8.38, 10.73, with spacing ≈ 2.2–2.3; the largest approaches ≈ 0.75.

## Asymptotics (`ASYMPTOTICS.md`)
- An elementary lemma (proved).
- A formal WKB derivation with explicit orders of the dropped terms.
- Positivity of C: closed form for Hou–Luo, checked pointwise in 2D.
- Numerical verification of every hypothesis, and what remains open: the exact C; 3/2; the o(1) remainder near the
  front.
- §8: the Hou–Luo λ → 1 limit problem. It gives the spacing constant π/C_∞ = 1.279 ± 0.001 (`hl_limit.py`,
  `hl_limit_nk.py`, `fig_limit.py`).

## The λ → 1 limit problem (Hou–Luo)
| step | command | time |
|---|---|---|
| family for Ω_f ≤ 10 (from the exact Ω_f = 0 solution) | `python3 hl_limit.py 0.5 1 1.5 2 2.5 3 3.5 4`, then `hl_limit_large.py 5 6 8 10` | ~15 min |
| Newton–Krylov continuation to Ω_f = 20 | `python3 hl_limit_nk.py hl_limit_Omf10.000.npz 10 20 0.5` | ~4 min |
| grid check (240 modes, y = 1 + σ³) | `HL_J=240 HL_M=1800 HL_NSIG=8000 HL_PSIG=3 HL_TAG=_fine python3 hl_limit_nk.py hl_limit_Omf10.000.npz 10 14 1` | ~2 min |
| table / figure | `python3 fig_limit.py` (and the table in `hl_limit_table.txt`) | ~1 min |
| front structure of the finite-ε states | `python3 hl_limit_look.py` → `hl_limit_look.out` | ~10 min |

## IPM (third model; blind test of the one-phase theory, `PREDICTIONS_IPM.md`)
The IPM self-similar problem is the Boussinesq problem with the vorticity slaved to the density gradient,
Ω = −∂₁R. It uses the same grid, Biot–Savart and Newton machinery (`ipm_solver.IPM` subclasses `bq_solver.BQ`).

| step | command |
|---|---|
| profile at one λ | `python3 ipm_solver.py 1.0` |
| continuation in z = 1/λ | `python3 ipm_branch.py LAM_START LAM_END DZ TAG START.npy` (`ipm_branch_*.out`) |
| rung (m = 2) | `python3 ipm_crossing.py Ya lam_a Yb lam_b` (`ipm_rung*.out`, states `ipm_rung_*.npy`) |
| stability operator T_μ | `ipm_stability.py` (time-translation check: `python3 ipm_stability.py STATE LAM`) |
| unstable spectrum (parity scan) | `python3 ipm_flow_spec.py STATE LAM NUHAT_MAX STEP -30 0.1` (`ipm_spec_U*.out`) |
| right-half-plane count | `python3 ipm_contour.py STATE LAM 0.025 TAG X_LO X_HI Y_HI` |
| counts used in the paper | `python3 ipm_contour.py ipm_rung_…_lamL.npy L 0.025 Un 0.04 1.5 Y 16 -30` with Y = 1.5 (U₀–U₂) or 2.5 (U₃, U₄); results in `ipm_contour_Un.out` (0.994, 1.947, 3.082, 4.077, 4.902 for n = 0–4; ~20–60 min each) |
| stage-2 predictions | `python3 ipm_stage2_predict.py` (data in `ipm_stage2_data.py`) |

### IPM deep branch, λ₇ holdout, phase observable and validation (NCS companion, `dossier/P03_MANUSCRIPT_NCS.md`)
| step | command / output |
|---|---|
| fine-grid scans at fixed z (resumable) | `NB=… HS_IN=… python3 ipm_scan.py …` → `ipm_scan_{s16,s20,p7,h7}.out/.npy` |
| grid-shift solves (front–grid locking) | `python3 ipm_shift.py STATE Z DELTA HS` (env NB, SS, SMAX, HS_IN, NK_TOL); queues `rich7_run*.sh` → `rich7_*.out`, `sys7_*.out` |
| λ₇ (pre-committed rule; converged grid) | `ipm_lambda7_analysis.py` → `…_committed_rule.out`; `ipm_lambda7_final.py` → `ipm_lambda7_final.out` |
| deep continuation (resumable) | `IPM_SSTART=-20 python3 ipm_deep.py LAM0 LAM1 DZ TAG START [Nb hs]` → `ipm_deep_e5.out` (to z = 10.23) |
| local eigenproblem | `ipm_local_eig.py` (shooting + secant) |
| phase, repaired tracker | `ipm_wkb.wkb_phase3`; tables `ipm_phase3_table.py` → `ipm_phase3_{rungs,l7,deep}.out` |
| tracker failure diagnostic | `ipm_tracker_diag.py` → `ipm_tracker_diag.out/.npz` |
| phase under the λ₇ audit variants | `ipm_phase_variants.py` → `ipm_phase_variants.out` |
| K = κD̂ checks (scaling, saturation) | `ipm_K_checks.py` → `ipm_K_checks.out/.npz` |
| travel time, anatomy, geometry | `ipm_travel_time.py`, `ipm_phase_anatomy3.py`, `ipm_endpoint_geometry.py` → `.out` |
| deep resolution checks | `ipm_deepres_check.py Z STATE_0125 STATE_00625 [Z_0125]` → `ipm_deepres_z*.out` |
| origin truncation / Newton floor | `ipm_sstart.py`, `ipm_noise2.py`, `ipm_shift.py` with SS = −26, −30 → `ipm_nfloor_*.out` |
| offset drift vs cut-off | `ipm_cutoff_drift.py` → `ipm_cutoff_drift.out` |
| independent global solver | `ipm_global.py` (Jacobian check: `python3 ipm_global.py`), `ipm_glob_run.py SRC LAM TAG HS NB SMIN SMAX,…` → `ipm_glob_*.log` |
| stage-4 predictions (λ₈–λ₁₀) | `ipm_stage4_predict.py` → `ipm_stage4_predictions.out` (registered in `PREDICTIONS_IPM.md`, 251bf1e) |
| figures | `fig_ncs.py [1–6]` → `fig_ncs*.png`; `fig_endings2.py` |

**Two independent solvers, IPM** (h_s = 0.0125, s_min = −16, s_max → ∞; λ₁ at h_s 0.025, s_min −12):
|Δλ| = 1.4e-8, 3.7e-8, 3.4e-8, 6.5e-7, 9e-7 and 1.3e-4 for n = 1–6. That is 0.2–3.5e-9 in m at every n, the
floor of both solvers. At h_s = 0.00625, λ₅ and λ₆ did not converge in the global solver, and those runs are not
used.
