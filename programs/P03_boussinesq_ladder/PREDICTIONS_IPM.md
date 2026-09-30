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
