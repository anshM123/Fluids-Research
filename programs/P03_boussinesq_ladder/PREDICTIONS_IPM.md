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
