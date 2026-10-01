# Blind predictions for 2D IPM (written before the IPM ladder or spectrum was computed)

**Status when written (2026-09-30).** The only IPM computation done so far is the solver test at λ = 1.0
(`ipm_solver.py`, one Newton solve from the generic initial guess):

- A = 1.5034172850;
- m = 2.0137632055;
- min V_r/r = 0.4964.

We have not continued in λ, located any rung, or run any IPM stability computation. We have not consulted any
published IPM eigenvalues. From the literature we use only the empirical fit 1/λ_n ≈ 1.1459 n + 0.9723 of
Wang et al. (arXiv:2509.14185), and the statement there that the profiles they found are ordered by their number
of unstable directions. Every prediction below comes from the one-phase theory built on Hou–Luo (HL) and 2D
Boussinesq (`INSTABILITY_LADDER.md` §7, `ASYMPTOTICS.md` §7). No IPM data were used.

## Setting

- **Equations.** IPM in the upper half-plane with no penetration at y₂ = 0 and even density:
  ρ = (1−t)^λ R(y), y = x/(1−t)^{1+λ}, with V·∇R = λR, V = (1+λ)y + U and Ω = −∂₁R.
- **Branch.** The least-singular branch has m(λ) = λ/(1+λ−A). A profile is smooth when m = 2.
- **Transported scalar.** Its exponent is ε_s = λ, where Boussinesq and HL have ε_s = λ − 1.
- **Perturbations.** They grow like e^{μτ}.
- **Trivial modes.** μ = 1 is time translation. μ = 0 is dilation, a Jordan block with ∂_λP.

## Theory being tested

- **One phase for rungs and spectrum.** A single WKB phase Φ(λ) of the profile fixes both the smooth rungs and
  the unstable spectrum.
  - Rungs: Re Φ(λ_n) = (n + c_r)π.
  - Spectrum: Re Φ(λ) − π μ/ε_s + θ₀ = (k + ½)π.
- **Where the μ-dependence comes from.** A mode leaving the stagnation point, where V₁ ≈ (λ/m) y₁, behaves like
  |y₁|^{m(1 − μ/λ)}. With m = 2 this gives the phase shift π μ/λ. So the unstable eigenvalues form a lattice of
  step ε_s = λ.
- **Where the rung law comes from.** Φ ≈ C/ε_s as ε_s → 0, so the rungs obey 1/λ_n ≈ (π/C) n + b.
- **Index.** Two lattices with the same step interlace, which gives index(U_n) = n.

## Predictions

**P1. The ladder is infinite, as in Boussinesq and HL, not finite as in CCF.**

- The branch continues from λ = 1 down to λ → 0 without a fold, a turning point or a sonic point.
- m(λ) − 2 changes sign infinitely often, and the sign changes accumulate only at λ = 0.
- min V_r/r stays positive for λ > 0 and tends to 0 like O(λ) as λ → 0. This is the stalled region; its front
  forms only at the endpoint.

*Falsified by* a fold, or by min V_r/r = 0 at some λ* > 0, beyond which there are no further rungs. That would
be the CCF-type finite ladder: a cusp at finite parameter.

**P2. Rung law.**

- 1/λ_n = a n + b + o(1).
- The local spacings Δ_n = 1/λ_n − 1/λ_{n−1} increase monotonically for n ≥ 2.
- They converge as Δ_n = a − κ λ_n + o(λ_n) with κ > 0. For comparison, κ ≈ 0.10 for HL and 0.09 for 2D.
- The asymptotic spacing a therefore exceeds every computed local spacing.

*Falsified by* non-monotone spacings, or by spacings that decrease with n.

**P3. Index (a consistency check, not blind).**

- U_n has exactly n eigenvalues with Re μ > 0, besides the time-translation mode μ = 1.
- These n eigenvalues are real.
- The dilation mode sits at μ = 0, split to about ±0.6/|s_start| by the origin truncation, as in 2D Boussinesq.

**P4. The unstable spectrum is a lattice of step λ_n.**

- The n unstable eigenvalues of U_n satisfy μ_{n,j} ≈ s_n λ_n (j + c_n), for j = 0, …, n−1.
- The step s_n → 1, with s_n − 1 = O(λ_n) > 0.
- The lower members are equally spaced. The top member sits above the lattice, as in both Boussinesq models.
- The rung offset c_n = μ_{n,0}/(s_n λ_n) lies in (0.5, 1), decreases with n, and converges.

The step s_n → 1 is the calibration-free number: it is the transported-scalar exponent and nothing else. By
analogy with HL (s ≈ 1 + 0.67/z) and 2D (s ≈ 1 + 1.0/z), we expect the lowest gap s_n ∈ [1.10, 1.40] at
n = 3 and 4, where z = 1/λ ≈ 3–4.5.

*Falsified by* unequal lower gaps (differing by more than 15% at n = 4), or by a step that does not approach
λ_n (for example, a step tending to λ_n/2 or to 2λ_n).

**P5. Spectral flow between rungs.**

- Along the continuous branch between U_n and U_{n+1}, the whole lattice moves rigidly with the profile phase:
  μ/λ ≈ s (j + c + f).
- f ∈ [0, 1) is the fraction of the interval in 1/λ.
- So d(μ/λ)/df ≈ 1 for every member.
- Exactly one eigenvalue crosses μ = 0 per interval, at f* ≈ 1 − c_n. This ties the place where the index
  changes to the offset measured at the rung.

**P6. Localization.** The unstable eigenfunctions of U_n (n ≥ 2) peak at the front of the profile, where |∇R| is
largest away from the origin. They do not peak at the stagnation point.

## Protocol

**Stage 1 (this file).** P1–P6 are fixed before any IPM continuation.

**Stage 2 (after the first rungs and their spectra).**

- Compute the rungs and spectra up to U₃ from this code, without further tuning.
- Commit numerical predictions for U₄ before computing it:
  - λ₄, from Δ₄ = a − κλ₄ fitted on Δ₂, Δ₃;
  - index 4;
  - its four unstable eigenvalues, from P4 with the s and c trend of U₁…U₃.
- Only then compute U₄. Do the same for U₅ if resolution allows.

**Recording.** Record every outcome, including failures, in `INSTABILITY_LADDER.md` and `ledger/RLOG.md`.

---

## Stage 1 outcomes so far (recorded 2026-09-30, before stage 2a)
**Rungs.** Solver: `ipm_branch.py`, then `ipm_crossing.py`.

| rung | λ | 1/λ | grid check |
|---|---|---|---|
| λ₀ | 1.0285722975 | 0.97222 | — |
| λ₁ | 0.4721297256 | 2.11806 | h_s = 0.0125 gives 0.4721297348 |
| λ₂ | 0.3149621637 | 3.17498 | — |

- Spacings: Δ₁ = 1.14584, Δ₂ = 1.05692.
- The empirical line through λ₀ and λ₁ (Wang et al.: 1/λ = 1.1459n + 0.9723) would put λ₂ at 1/λ = 3.264. The
  computed rung is 2.7 % lower.

**P3 (index), exact so far.**
- U₀: only the time-translation mode (μ = 0.999995).
- U₁: ν̂ = μ/λ = 0.8608 (μ = 0.406428). The argument-principle count on [0.04, 1.5] × [±1.5] is 1.947 → 2
  (trivial + 1).
- U₂: ν̂ = 0.8247, 1.9444 (μ = 0.259752, 0.612423). The count is 3.082 → 3 (trivial + 2).
- All unstable modes are real.

**P1, consistent so far.**
- The branch continues without a fold to 1/λ = 3.7.
- min V_r/r = 0.496, 0.325, 0.203, 0.117 at λ = 1, 0.714, 0.490, 0.315, i.e. ≈ 0.40λ, decreasing.
- m − 2 oscillates with extrema +0.081, −0.0077, +0.00077: a factor ≈ 10 per half period.

**P2.** Not yet tested; its monotonicity clause concerns n ≥ 2.

## Stage 2a — predictions for λ₃, U₃ and the spectral flow between U₁ and U₂
These are fixed before λ₃ or any between-rung spectrum is computed. They are produced by
`ipm_stage2_predict.py` from `ipm_stage2_data.py` (λ₀–λ₂; spectra of U₁, U₂); output in
`ipm_stage2a_predictions.out`.

**Disclosure.** The branch continuation already ran past λ₂, to 1/λ = 3.7 (`ipm_branch_dn2.out`). Those m values
are not used: the rung prediction uses only the three rung positions.

**λ₃.**
- The P2 law Δ_n = a − κ/z̄_n, fitted on Δ₁ and Δ₂, gives a = 0.932 and κ = −0.330. Hence Δ₃ = 1.022 and
  **λ₃ = 0.2383 (1/λ₃ = 4.197)**.
- This fit has κ < 0 (spacings decreasing), which contradicts the monotonicity clause of P2 (κ > 0). Both
  commitments are recorded, and at most one can hold.
  - The fit predicts Δ₃ < Δ₂.
  - P2's clause predicts Δ₃ > Δ₂, i.e. λ₃ < 0.2363.
- For reference, the empirical line gives λ₃ = 0.2268.

**U₃ spectrum.**
- **index(U₃) = 3**, all real.
- **ν̂ = μ/λ₃ ≈ [0.807, 2.005, 3.169]**, i.e. μ ≈ [0.192, 0.478, 0.755] at the predicted λ₃. The ingredients:
  - lowest member 0.807, extrapolated linearly in 1/z through U₁ and U₂;
  - lower gap s = 1 + β/z with the prior β = 0.83 ± 0.17 (Hou–Luo 0.67, 2D 1.0), so s = 1.198 (range 1.157–1.238);
  - top gap 1.164, i.e. U₂'s gap grown by 4 %.
- The parameter-free lattice (s = 1, offset of U₂) would give ν̂ = [0.825, 1.825, 2.825].
- **Falsified if:** the index is not 3; or a complex pair appears; or the lower gap of U₃ differs from s ∈ [1.10, 1.30]
  by more than its range; or the lowest member is outside [0.72, 0.90].

**P5, spectral flow inside [z₁, z₂] = [2.1181, 3.1750].** The rigid lattice moves linearly in the fraction f,
with slope 1.084 per interval, and a new member enters through μ = 0 at **f* = 0.239 (z = 2.371)**. Predicted
non-trivial unstable ν̂ = μ/λ at three branch points:

| f | z | λ | predicted ν̂ |
|---|---|---|---|
| 0.25 | 2.3823 | 0.4197637 | [0.012, 1.132] (the lowest member has just entered) |
| 0.50 | 2.6465 | 0.3778542 | [0.283, 1.403] |
| 0.75 | 2.9108 | 0.3435536 | [0.554, 1.674] |

**Falsified if:** the members are not near-linear in f (deviation > 0.1 in ν̂ at f = 0.5), or the entry point is
outside f ∈ [0.1, 0.4].

## Stage 2a outcomes, part 1 (recorded before stage 2b)
**λ₃.**
- Computed: **λ₃ = 0.2415660984 (1/λ₃ = 4.13965)**.
- Predicted: 0.2383 (4.197). The error is 0.057 in 1/λ, i.e. 6 % of a spacing, or 1.4 % in λ.
- For reference, the empirical line gives 1/λ₃ = 4.410, an error of 0.27 (28 % of a spacing).

**P2 monotonicity clause: falsified for IPM.**
- The spacings decrease: Δ₁, Δ₂, Δ₃ = 1.146, 1.057, 0.965. So Δ₃ < Δ₂, whereas in Hou–Luo and 2D Boussinesq they
  increase.
- The fitted P2 law (κ < 0) anticipated the decrease, but underestimated it.
- Whether the spacing converges to a positive constant (P1 and P2's linear law) or keeps falling is open. The
  stalled-layer diagnostic favours convergence: min V_r/r ∝ λ^{1.28} over λ = 1 → 0.24 (min V_r/r / λ = 0.50 → 0.33).
  That is a slowly deepening dip, as in 2D Boussinesq, with no sign of a sonic point at finite λ.

**P5 at f = 0.25** (λ = 0.4197637).
- Computed ν̂ = [1.1331] above the artifact level ν̂ ≈ 0.05, plus the trivial mode.
- Predicted [0.012, 1.132]:
  - the continuing member is off by 0.001;
  - the entering member, predicted at 0.012, lies below the resolvable level, as expected.

## Stage 2b, part 1 — λ₄ (fixed before the branch reaches it; the branch is capped at 1/λ = 4.8)
- Rule of stage 2: Δ = a − κ/z̄ on (Δ₂, Δ₃). This gives a = 0.723 and κ = −0.883, hence Δ₄ = 0.915 and
  **λ₄ = 0.19783 (1/λ₄ = 5.0549)**.
- Other simple extrapolations span λ₄ = 0.1959 (constant spacing) to 0.1995 (linearly decreasing spacing).
- The empirical line gives 0.1800.

The U₄ spectral predictions (stage 2b, part 2) follow once the U₃ spectrum is computed, and before λ₄ is computed.

## Stage 2a outcomes, part 2 (recorded before stage 2b part 2)
**U₃ spectrum** (`ipm_spec_U3.out`). ν̂ = μ/λ₃ = [0.8054, 1.8794, 2.9371], i.e. μ = [0.19456, 0.45401, 0.70949].
The trivial mode sits at μ = 1.000011.

| quantity | predicted | computed | verdict |
|---|---|---|---|
| index(U₃) | 3 | 3, all real (contour count running) | ✓ |
| lowest member | 0.807 | 0.8054 | ✓ (0.2 %) |
| lower gap | 1.198 (prior range 1.157–1.238, from Hou–Luo and 2D) | 1.074 | ✗ |
| top gap | 1.164 | 1.058 | ✗ |
| middle / top members | 2.005 / 3.169 | 1.879 / 2.937 | ✗ (off by 0.13 / 0.23) |
| parameter-free lattice | [0.825, 1.825, 2.825] | [0.805, 1.879, 2.937] | closer: 0.02 / 0.05 / 0.11 |

The IPM lattice is *closer* to the asymptotic step λ_n than the prior taken from the other two models: s₃ − 1 =
0.074, so β ≈ 0.31 instead of 0.67–1.0. Nor is the top member displaced. The calibration transferred from Hou–Luo
and 2D failed; the prediction of the theory itself, step → λ, fared better.

**P5, spectral flow inside [z₁, z₂]** (`ipm_spec_bp025/050/075.out`).

| f | predicted ν̂ | computed ν̂ |
|---|---|---|
| 0.25 | [0.012, 1.132] | [1.1331] (entering member below the resolvable 0.05) |
| 0.50 | [0.283, 1.403] | [0.2384, 1.4038] |
| 0.75 | [0.554, 1.674] | [0.5347, 1.6741] |

- **Continuing member:** within 0.001 at all three points.
- **Entering member:** 0.045 and 0.019 behind the rigid prediction. It moves slightly faster (slope 1.17 per
  interval against 1.08) and entered at f* ≈ 0.30, against 0.24 predicted.
- **Falsification criteria:** both are met, i.e. the prediction survives: linearity is within 0.1, and f* lies in
  [0.1, 0.4].

## Stage 2b, part 2 — U₄ and the spectral flow inside [z₂, z₃] (fixed before λ₄ is computed)
From `ipm_stage2b_predictions.out`, i.e. `ipm_stage2_predict.py` on λ₀–λ₃ and U₁–U₃.

**U₄.**
- **index(U₄) = 4**, all real.
- **ν̂ = μ/λ₄ ≈ [0.794, 1.854, 2.915, 3.911]**. The ingredients:
  - lowest member linear in 1/z;
  - lower gap s = 1 + β/z with β = 0.306 from U₃'s lowest gap, giving s = 1.061;
  - top gap linear in n, giving 0.996.
- At the predicted λ₄ = 0.19783, μ ≈ [0.157, 0.367, 0.577, 0.774].
- The parameter-free lattice (s = 1, offset of U₃) gives [0.805, 1.805, 2.805, 3.805].

**P5 inside [z₂, z₃] = [3.1750, 4.1397]** (rigid lattice, slope 0.993).
- The new member enters at f* = 0.189 (z = 3.357).
- Predicted ν̂ at the branch points:

  | f | λ | predicted ν̂ |
  |---|---|---|
  | 0.25 | 0.2927270 | [0.061, 1.135, 2.193] |
  | 0.50 | 0.2734243 | [0.309, 1.383, 2.441] |
  | 0.75 | 0.2565098 | [0.557, 1.631, 2.689] |

## Stage 2b outcome, part 1, and stage 2c (λ₅), fixed before the branch reaches λ₅
**λ₄.**
- Computed: **λ₄ = 0.198730 (1/λ₄ = 5.0320)**. The refinement is converging; the final digits are in
  `ipm_rung4.out`.
- Predicted: 0.19783 (5.0549). The error is 0.023 in 1/λ (2.5 % of a spacing), or 0.45 % in λ.
- The spacings keep decreasing (1.146, 1.057, 0.965, 0.892), but the decrements shrink (−0.089, −0.092, −0.072).
  This is consistent with convergence to a positive asymptotic spacing.

**Stage 2c: λ₅** (`ipm_stage2c_predictions.out`).
- The same rule, fitted on (Δ₃, Δ₄), gives **λ₅ = 0.17010 (1/λ₅ = 5.879)**.
- Simple extrapolations bracket 0.1688–0.1707.
- The empirical line gives 0.1492.

## Stage 2c outcome and stage 2d (λ₆), recorded 2026-09-30 before the branch reached λ₆
**λ₅.**
- Computed: **λ₅ = 0.1708490449 (1/λ₅ = 5.8531)**.
- Predicted: 0.17010 (5.879). The error is 0.026 in 1/λ (3 % of a spacing), or 0.44 % in λ.

**Stage 2d: λ₆** (`ipm_stage2d_predictions.out`).
- The rule gives **λ₆ = 0.150928 (1/λ₆ = 6.626)**.
- Simple extrapolations give 0.1498–0.1513.

**Observation, recorded before λ₆ is computed.** The spacings 1.146, 1.057, 0.965, 0.892, 0.821 shrink by a nearly
*constant factor*, 0.922, 0.913, 0.925, 0.920. The fitted asymptotic spacing of the P2 law therefore keeps falling
with every new rung (0.93, 0.72, 0.61, 0.44). This is the signature of rungs accumulating at a *finite* λ, not
at λ → 0.

It is corroborated by the stalled-layer diagnostic:
- min V_r/r / λ falls roughly linearly in 1/λ, from 0.50 at λ = 1 to 0.24 at λ = 0.156;
- its local log-slope rises from 0.25 to above 1;
- extrapolated, it reaches zero near λ ≈ 0.08–0.09.

The minimum is the dip of the boundary speed just behind the front, at x ≈ 0.61.

If this persists, P1 is **falsified** for IPM: the branch would end at a finite-λ sonic cusp, as its 1D reduction
CCF does, rather than in a stalled layer at λ → 0. The branch is being continued to 1/λ ≈ 12 to test this
(`ipm_branch_dn5.out`).

**Stage 2d′ (added before λ₆ is computed): a shifted accumulation point.**
- Fitting an accumulation point λ_c so that 1/(λ_n − λ_c) has constant spacing gives **λ_c = 0.0379**. The
  spacings are then 1.294, 1.306, 1.301, 1.308, 1.304, constant to 0.2 %.
- The same procedure recovers λ_c = 0.9992 (2D) and 0.9993 (Hou–Luo), where the known value is 1.
- This law predicts **λ₆ = 0.1512, λ₇ = 0.1367, λ₈ = 0.1254**.
- The geometric model (ratio 0.92) gives nearly the same: 0.1513, 0.1369, 0.1259.
- Both imply rungs accumulating at a finite λ (0.038 or 0.065), not at λ → 0.

## Stage 2d outcome (λ₆, h_s = 0.025) and a resolution warning
**λ₆.** Computed λ₆ = 0.1517666 (1/λ₆ = 6.589). Errors in 1/λ against the predictions:

| rule | predicted λ₆ | error in 1/λ |
|---|---|---|
| stage rule | 0.150928 | −0.037 |
| shifted law (λ_c = 0.0379) | 0.1512 | −0.025 |
| geometric | 0.1513 | −0.020 |

That is 3–5 % of a spacing. Spacing Δ₆ = 0.736.

**Resolution warning.**
- Beyond λ₅ the smoothness defect is |m − 2| ≲ 1e-6.
- On the h_s = 0.025 grid, m jitters by ~1e-7 from step to step (`ipm_branch_dn5.out`). This is visible as a spurious
  "crossing" at 1/λ = 6.905, with |m − 2| growing to 1.6e-6.
- So λ₆ on this grid is uncertain by roughly 0.02–0.05 in 1/λ, and the deep branch run was stopped.
- λ₅ and λ₆ are being recomputed at h_s = 0.0125. As in 2D Boussinesq, the deep rungs need the finer grid.

## Stage 2b outcome, part 2 — U₄ spectrum (`ipm_spec_U4.out`)
Computed ν̂ = μ/λ₄ = [0.7958, 1.8463, 2.8791, 3.8684], i.e. μ = [0.15814, 0.36692, 0.57215, 0.76877]. The trivial
mode sits at μ = 1.000025.

| quantity | predicted | computed | error |
|---|---|---|---|
| index(U₄) | 4 | 4, all real; contour count 4.902 → trivial + 4 (`ipm_contour_U4.out`) | ✓ |
| ν̂ members | [0.794, 1.854, 2.915, 3.911] | [0.7958, 1.8463, 2.8791, 3.8684] | +0.002, −0.008, −0.036, −0.043 (0.2–1.2 %) |
| μ at the computed rung | [0.157, 0.367, 0.577, 0.774] | [0.158, 0.367, 0.572, 0.769] | within 0.7 % |
| parameter-free lattice | [0.805, 1.805, 2.805, 3.805] | as above | −0.01, +0.04, +0.07, +0.06 |

- Lattice gaps: 1.051, 1.033, 0.989. The step tends to λ_n (s₄ − 1 = 0.05), and the top member is not displaced.
- **Fine-grid check of λ₅** (h_s = 0.0125): 0.1706181, against 0.1708490 at h_s = 0.025. The shift is 0.008 in
  1/λ, i.e. 1 % of a spacing. The coarse-grid rungs n ≤ 5 are therefore accurate to about 1 % of a spacing.

## Stage 2d outcome at h_s = 0.0125 (λ₆)
**Result.**
- On the fine grid, λ₆ = 0.15092 ± 0.00002 (1/λ₆ = 6.626). The secant refinement wanders in 1e-9 noise of m, so the
  uncertainty comes from the bracketing points in `ipm_rung6_hs0125.out`.
- The stage-2d prediction was **0.150928**, which agrees to 1e-5.
- The coarse-grid value 0.15177 was off by 5 % of a spacing. That was the resolution error flagged above.

**The IPM ladder, n ≤ 4 from h_s = 0.025 and λ₅, λ₆ from h_s = 0.0125.**
- In 1/λ the spacings are 1.146, 1.057, 0.965, 0.892, 0.829, 0.765, with ratios 0.92–0.93.
- In 1/(λ − λ_c) with **λ_c = 0.0370** they are 1.290, 1.300, 1.291, 1.295, 1.301, 1.295: constant to 0.3 % over six
  intervals.
- The one-phase linear law therefore holds for IPM, with a *finite* accumulation point λ_c = 0.037, not λ = 0.

**Scorecard of P1–P6.**

| prediction | outcome |
|---|---|
| P1 (no fold or sonic point at λ > 0; accumulation as λ → 0; min V_r/r = O(λ)) | **failed** in its accumulation and O(λ) clauses. No fold was met down to λ = 0.14, but min V_r/r / λ falls from 0.50 to 0.21 and extrapolates to zero near λ ≈ 0.09. |
| P2 (linear law; increasing spacings) | linear law **holds in 1/(λ − λ_c)**; the increasing-spacing clause **failed** (spacings decrease) |
| P3 (index n, all real) | **held**, n = 0–4. Contour counts on [0.04, 1.5] × [±1.5 or ±2.5]: 0.994, 1.947, 3.082, 4.077, 4.902 for n = 0–4, i.e. trivial + n |
| P4 (lattice of step λ_n, s → 1) | **held**: gaps 1.12 (U₂); 1.07, 1.06 (U₃); 1.05, 1.03, 0.99 (U₄) |
| P5 (rigid spectral flow, entry at μ = 0) | **held**: continuing member within 1e-3 at 4 branch points; others within 0.05 |
| P6 (localization at the front) | **held**: U₃ modes peak at the dip x ≈ 0.61, regular at the stagnation point |

**Stage-2 numerical predictions.**

| predicted | from | error |
|---|---|---|
| λ₃ | λ₀–λ₂ | 6 % of a spacing |
| λ₄ | λ₀–λ₃ | 2.5 % |
| λ₅ | λ₀–λ₄ | 3 % (coarse) / 2 % (fine) |
| λ₆ | λ₀–λ₅ | 0.1 % (fine) |
| U₃, lowest member | — | 0.2 % |
| U₃, gaps (calibration from other models) | — | ✗ |
| U₄, all four members | — | 0.2–1.2 % |
| spectral flow | — | within 1e-3 (continuing member) |

---

## Stage 3 — post-outcome holdout (registered before λ₇ and λ₈ are computed)
*The stage-1 and stage-2 predictions above stay as they were, including the failed ones. Nothing below edits them.*

**Why a stage 3.** The accumulation at λ → 0 (P1, P2) failed. Three readings of the finite-λ behaviour remain:
- (i) infinitely many profiles accumulating at λ_c > 0;
- (ii) a finite ladder ending where the wall dip closes (a sonic point, as in CCF);
- (iii) a numerical degeneration of the branch.

Stage 3 tests them with pre-registered numbers.

**The ladder on one grid.** All rungs recomputed at h_s = 0.0125 (`ipm_rung{2,3,4}_hs0125.out`); λ₀ from h_s = 0.025.

| n | λ_n (h_s = 0.0125) | shift from h_s = 0.025 |
|---|---|---|
| 0 | 1.0285722975 (h_s = 0.025) | — |
| 1 | 0.4721297348 | +9e-9 |
| 2 | 0.3149618108 | −3.5e-7 |
| 3 | 0.2415663353 | +2.4e-7 |
| 4 | 0.1987224523 | −7.5e-6 |
| 5 | 0.1706180880 | −2.3e-4 |
| 6 | 0.15092 ± 0.00002 | −8.5e-4 |

**λ_c from successive subsets** (`ipm_stage3_predict.py` → `ipm_stage3_predictions.out`).
- 3-point windows (0:2) … (4:6): 0.0341, 0.0396, 0.0357, 0.0353, 0.0389.
- 4-point windows: 0.0368, 0.0376, 0.0355, 0.0371.
- The λ₆ uncertainty moves the (4:6) value by ±0.0005.
- So λ_c = 0.037 ± 0.003 on the present data, with no clear trend yet.

**Predictions for λ₇ and λ₈ (fixed rules, no refitting after λ₇ is known).**

| rule | λ₇ | z₇ | λ₈ | z₈ |
|---|---|---|---|---|
| H1: 1/(λ_n − λ_c) = an + b on λ₂–λ₆ (λ_c = 0.0365) | 0.13621 | 7.3415 | 0.12486 | 8.0087 |
| H1′: the same on λ₀–λ₆ (λ_c = 0.0367) | 0.13623 | 7.3404 | 0.12490 | 8.0063 |
| H2: geometric contraction in z (ratio 0.9258) | 0.13635 | 7.3343 | 0.12516 | 7.9899 |
| H3: 3-point law from (λ₄, λ₅, λ₆) (λ_c = 0.0389) | 0.13635 | 7.3342 | 0.12513 | 7.9917 |

These rules differ by 1 % (λ₇) to 3 % (λ₈) of a spacing. The rung positions therefore test the rung law at the 1 %
level; they cannot decide the endpoint.

**The endpoint test: the phase measured directly.**
- New tool: the IPM local eigenproblem (`ipm_local_eig.py`, `ipm_wkb.py`). It gives κ(s) along the wall of any
  computed profile, and Φ = (1/D₀)∫κ ds with D₀ = λ/m.
- Check on the known rungs (front cut-off D/D₀ = 3, `ipm_wkb_rungs_cut3.out`):
  - Re Φ = 5.50, 8.58, 11.69, 14.82, 17.93 (fine), 20.86 (coarse λ₆) for n = 1–6;
  - increments 3.08, 3.11, 3.13, 3.11 against π.
  - The phase therefore quantizes the IPM rungs, as in Boussinesq and Hou–Luo.
- The phase coefficient I(z) = ∫Re κ ds = D₀ Re Φ: 1.298, 1.351, 1.412, 1.473, 1.530, 1.583 at z = 1/λ = 2.12–6.59.

What each reading implies for I(z) on the deep branch (`ipm_deep.py` states, z up to 25):

| reading | I(z) as z grows | profiles |
|---|---|---|
| accumulation at λ = 0 with a linear law (original P2) | tends to a constant | 1/λ_n ≈ an + b |
| accumulation at λ = 0, faster phase | grows linearly in z | λ_n ∝ n^{−1/2} |
| accumulation at λ_c ≈ 0.037 (H1) | ∝ 1/(1 − λ_c z), diverging at z ≈ 27 | 1/(λ_n − λ_c) ≈ an + b |
| termination at the dip closure (ii) | the branch ends where min V_r/r → 0 (λ ≈ 0.08–0.09 extrapolated) | finite ladder |

Discriminating values:
- At z = 15, "linear in z" gives I ≈ 2.1 and H1 gives I ≈ 2.7.
- At z = 20, "linear in z" gives I ≈ 2.4 and H1 gives I ≈ 4.6.

**Protocol.**
- A phase-based prediction of λ₇ (H4, from Φ on the deep coarse states, which do not resolve the defect) will be
  committed separately.
- The fine-grid computation of λ₇ starts only after that commit.
- λ₈ is attempted at h_s = 0.0125 or 0.00625, depending on the measured noise floor of m (`ipm_noise.py`).
- The defect amplitude falls by about 10 per half-period, from 7e-6 at z ≈ 5.3 to about 5e-8 near λ₇ and 5e-9 near
  λ₈.

### Stage 3, H4 — phase-based prediction (registered before any fine-grid computation below λ₆)
**The phase at the rungs, cut-off D/D₀ = 2** (`ipm_wkb_rungs_fine.out`; λ₆ from the coarse phase curve interpolated
to the fine λ₆, `ipm_wkb_coarse_cut2.out`).

| n | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| δ_n = Re Φ(λ_n) − nπ | 2.038 | 2.022 | 2.014 | 2.058 | 2.057 | 2.037 |

- δ_n is constant to ±0.02; the mean is δ = 2.038.
- The mean Re Φ increment is 3.146, against π (0.1 %). With cut-off 3 or 5 the δ_n drift by −0.11 and −0.25.
- So cut-off 2 is the one that makes the quantization exact. The same held in 2D (`ASYMPTOTICS.md` §5.2).

**Rule H4.** Re Φ(λ_n) = nπ + 2.038, with Φ from the coarse-grid (h_s = 0.025) deep branch (`ipm_deep.py`, tag c1).

| | λ | z | confidence |
|---|---|---|---|
| **λ₇** | 0.13601 | 7.352 (quadratic interpolation 7.3523–7.3543) | usable |
| **λ₈** | ≈ 0.1263 | 7.91–7.92 | **low** |

- **Why λ₈ is low-confidence.** Beyond z ≈ 7.6 the coarse grid does not resolve the steepening front.
  - The identity Ω_b = −∂ₓR fails by 8 % there.
  - Re Φ jumps by 2.30 between z = 7.87 and 8.12.
  - Root-tracking failures appear.
- **What λ₈ could decide.** H4 and H1–H3 differ by 2 % of a spacing at λ₇ and by about 12 % at λ₈. So λ₈, if it
  can be resolved, discriminates between the phase-based and the fit-based rules.

**Resolution facts measured so far (not predictions).**
- Newton floor on the fine grid: |R| ≈ 4e-13. Tightening the tolerance from 1e-11 to 1e-13 changes m by 1.7e-11
  (`ipm_noise_lam6_hs0125.out`). The 1e-9 scatter of m near λ₆ is therefore discretization, not the solver.
- A shorter domain, s ∈ [−12, 40] instead of [−20, 100], changes m at λ₅ by 2e-6 and is rejected. An s_max = 130
  check is running.

### Stage 3 — measured systematics and the direct phase (recorded before λ₇ is known)
**Systematic uncertainty of the deep rungs.** All checks at λ₅ and λ₆ are on the fine grid (h_s = 0.0125)
(`ipm_sstart_lam{5,6}.out`, `ipm_nbcheck*.out`, `ipm_rung{4,5,6}_Nb48.out`, `ipm_rung5_hs0125_smax130.out`,
`ipm_rung6_hs00625.out`, `ipm_noise*.out`).

| source | effect on m at λ₆ | effect on λ₆ |
|---|---|---|
| Newton floor: 5e-13 at s_start = −20 (1e-13 at −16, 1e-11 at −30) | scatter of about 5e-10 | ±2e-5 |
| origin cut-off (s_start −19 vs −20) | 3e-10 | converged |
| angular resolution (Nb 48 vs 32 at s_start = −20) | ≤ 1e-9 (−1.0e-9, +6e-10, −3e-10 at λ₄, λ₅, λ₆) | ≤ 4e-5 |
| far field (s_max 130 vs 100) | 2e-10 | negligible |
| radial resolution (h_s 0.00625 vs 0.0125) | ≈ −7e-10 (noise-limited) | ≈ −2.5e-5 ± 4e-5 |

- **Result: λ₆ = 0.15092 ± 0.00005, i.e. ±0.25 % of a spacing.**
- **Origin cut-off in detail.** The bias is not monotonic: m changes by 2e-7 (s_start = −14), 2e-8 (−15), −5e-9 (−16)
  and −1.4e-9 (−17), all relative to −20. It is the same at λ₅ and λ₆.
- **Angular resolution in detail.** At the shallow cut-off s_start = −16, Nb 32 → 48 shifts m by 2e-8 (the same at λ₅
  and λ₆), and Nb 64 agrees with 48 to 9e-11. At s_start = −20 the angular effect disappears. The production
  settings (Nb 32, h_s 0.0125, s_start −20) are therefore converged.

**λ_c with a free exponent.** Fitting (λ_n − λ_s)^{−γ} = an + b:

| rungs | γ | λ_s |
|---|---|---|
| λ₀–λ₆ | 0.996 | 0.0372 |
| λ₁–λ₆ | 1.013 | 0.0356 |
| λ₂–λ₆ | 0.914 | 0.0437 |

- With γ fixed at ½ (the law a closing dip would give; see below), the residuals are 5–15 times larger.
- So over λ₁–λ₆ the data prefer the shifted linear law, with λ_c = 0.036–0.044.

**The direct phase, continuous cut-off.**
- The cut-off D/D₀ = 2 is now interpolated between grid points (`ipm_wkb_fine_rungs_cut2c.out`,
  `ipm_wkb_fine_l6_cut2c.out`); this removes the grid-step jitter of the phase.

| n | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| δ_n = Re Φ(λ_n) − nπ | 2.033 | 2.023 | 2.032 | 2.037 | 2.028 | 1.984 |

- For n = 1–5 the phase quantizes the rungs to ±0.007, i.e. increments of π within 0.3 %.
- At n = 6 the increment is 3.097, 1.4 % below π.

**An exact scaling of the IPM local problem.**
- Rescaling Y → Y/D̂ removes D̂ from the local eigenproblem, so κ = K(ĉ, μ, G, R_y)/D̂ exactly. Numerically D̂κ is
  constant to four digits for D̂ = 0.6 → 0.012.
- Hence **Φ = ∫K(s) ds / D(s)**: the phase is the integral of a local wavenumber K over the unnormalized
  self-similar wall speed D = V₁/x.
- This holds only for IPM. In Boussinesq the vorticity transport term breaks the scaling.
- Consequence: if the wall dip closes (D_min → 0), Φ diverges, like D_min^{−1/2} for a quadratic minimum.
- On the fine grid, D̂_min = 0.819, 0.735, 0.666, 0.597, 0.528, 0.452 at λ₁–λ₆, and 0.438, 0.426 at z = 6.75, 6.87.
  Linear extrapolation closes the dip near z ≈ 11.1–11.4 (λ ≈ 0.087–0.090). If that happens, the accumulation law
  must change before λ ≈ 0.09 (towards γ = ½). The resolved continuation `ipm_deep.py` (tag e1) tests this.

## Stage 3 outcome: λ₇ (holdout of H1–H4)
**The committed rule could not complete.**
- The decision rule `ipm_lambda7_analysis.py` (committed f0b3f06, before the crossing was known) takes σ_z from the
  front–grid scatter of the production scan p7.
- p7 was paused before its crossing to test the front–grid locking hypothesis, so the rule's σ path could not run
  (`ipm_lambda7_analysis_committed_rule.out`).
- Its other branch gives, from the h_s = 0.00625 scan alone, z₇ = 7.3457 with σ_z(noise) = 0.004.

**Converged-grid determination** (`ipm_lambda7_final.py`; h_s = 0.00625, four grid shifts per point).

| z | m − 2, mean of 4 shifts |
|---|---|
| 7.256 | +4.89e-9 ± 3e-11 |
| 7.316 | +9.1e-10 ± 4e-11 |
| 7.346 | −2.3e-11 ± 4e-11 |
| 7.376 | −7.2e-10 ± 3e-11 |

- The locking error is ±9e-11 at h_s 0.00625, against ±3e-9 to ±1.1e-8 at h_s 0.0125 and ±1.5e-6 at 0.025.
- Result: **z₇ = 7.3451 ± 0.0016 (statistical), λ₇ = 0.136145.**

**Systematics at the crossing** (single solves at z = 7.346; `sys7_*.out`):

| change | effect on m − 2 |
|---|---|
| Nb 32 → 48 → 64 | +3e-11 → −1.10e-9 → −3.39e-9 (not converging) |
| s_start −20 → −22 → −24 | +3e-11 → −4.8e-10 → −1.48e-9 (the −24 value equals the λ₅ offset: a round-off artifact of deep truncation) |
| s_max 100 → 130 → 160 | +3e-11 → +6.2e-10 → +3.6e-11 (the 130 value is solver-path noise) |
| Nb 48, s_start −22, s_max 130 together | −2.03e-9 (not additive) |

- The formulation therefore carries an erratic error floor of about 1–3 × 10⁻⁹ in m at this depth.
- The defect crosses zero with slope −2.5 × 10⁻⁸ per unit z, because the next lobe is ten times smaller.
- Hence **z₇ = 7.35 ± 0.08 from the defect**.

**Verdict under the registered rule.**

| rule | z₇ predicted | distance from measurement |
|---|---|---|
| H1 | 7.3415 | −0.04σ |
| H1′ | 7.3404 | −0.06σ |
| H2 | 7.3343 | −0.13σ |
| H3 | 7.3342 | −0.14σ |
| H4 | 7.352 | +0.09σ |

σ_z = 0.08 > 0.003, so the λ₇ holdout is **not decisive**: no hypothesis is preferred or disfavoured.

**The same profiles, read through the phase** (`ipm_wkb_h7b_cut2.out`; not a registered prediction).
- Re Φ = 23.858, 24.025, 24.166 at z = 7.316, 7.346, 7.376.
- With δ = 2.031 (the n = 1–5 mean) the phase places λ₇ at **z₇ = 7.3456 ± 0.003**, about 25 times more precisely than
  the defect.
- At the defect's central value, δ₇ = 2.03.

**The general lesson** (Fig. 4 of the manuscript).
- Rung n is located to σ_z ≈ σ_m / |dm/dz| at its crossing.
- |dm/dz| there scales with the next, smaller lobe, so it falls about twentyfold per profile (5.2 × 10⁻⁷ at λ₆,
  2.5 × 10⁻⁸ at λ₇).
- The discretization floor σ_m stays at about 10⁻⁹.
- Direct location therefore loses about a factor of twenty per profile (σ_z ≈ 0.002 at λ₆, 0.08 at λ₇, of order 1
  at λ₈), while the phase keeps σ_z ≈ 0.003.
