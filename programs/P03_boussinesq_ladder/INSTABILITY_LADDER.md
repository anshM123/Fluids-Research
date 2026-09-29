# The instability ladder: counts, spectra and localization of the unstable modes

*2D Boussinesq: linearization A = `bq_stability.py` (march-based map T_μ) and linearization B = `glob_stab.py`
(global shift-invert). Hou–Luo: `hl_stability.py` (map T_μ), validated by ν = 1.000000 at μ = 1 for the exact
time-translation mode. Perturbations are ∝ e^{μτ}, with τ = −ln(1−t).*

## 1. Right-half-plane counts (argument principle for det(I − T_μ))
The counts include the trivial time-translation zero μ = 1, so the expected value for rung n is n + 1.

**2D Boussinesq** (`stab_contour.py`, k = 12; `stab_contour2.py`, k = 16):
| n | variant (hs, origin truncation, rectangle) | zeros | n + 1 |
|---|---|---|---|
| 1 | 0.025, −20, [0.06,1.5]×[−1.5,1.5] | 1.997 | 2 |
| 2 | same | 3.000 | 3 |
| 3 | same | 4.000 | 4 |
| 4 | same | 5.000 | 5 |
| 5 | A: 0.025, −20, [0.06,1.5]×[−1.5,1.5] | 5.952 | 6 |
| 5 | B: 0.025, −30, [0.075,1.5]×[−1.5,1.5] | 5.954 | 6 |
| 6 | A: 0.025, −20, [0.06,1.5]×[−1.5,1.5] | 6.942 | 7 |
| 6 | B: 0.025, −30, [0.07,1.5]×[−1.5,1.5] | 6.945 | 7 |
| 6 | C: 0.0125, −20, [0.065,1.2]×[−0.8,0.8] | 7.003 | 7 |
| 7 | A: 0.025, −30, [0.05,1.5]×[−1.5,1.5] | 7.921 | 8 |
| 7 | C: 0.0125, −30, [0.055,1.2]×[−0.8,0.8] | 7.973 | 8 |

Deviations from an integer (≤ 0.08) come from eigenvalue swaps at the truncation of the top-k set. At the
unresolved steps the smallest |1 − ν| is ≥ 0.12, so no zero lies near the contour.

**Hou–Luo** (`hl_contour.py`, k = 20, [0.03,1.5]×[−1.5,1.5], N = 32768, η₀ = −30; n ≤ 2 at N = 8192):
| n | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| zeros | 1.991 | 3.000 | 4.000 | 5.000 | 6.000 | 7.000 | 8.000 | 9.000 | 9.959 | 11.000 |

Variants:
- Origin truncation η₀ = −40 with x_lo = 0.04: 6.000, 7.024, 8.000, 9.032, 9.963, 11.000 for n = 5…10.
- Grid N = 65536 (profiles re-converged, |m − 2| ≤ 6e-11; η₀ = −30): 9.000 (n = 8) and 10.989 (n = 10).

**Result.** In the rectangle, the n-th profile has exactly n unstable eigenvalues: n ≤ 7 in 2D and n ≤ 10 in
Hou–Luo. For 2D n = 5–7 the count is unchanged under a deeper origin truncation with a moved contour (B) and under
grid refinement hs 0.025 → 0.0125 (C). The strip 0 < Re μ < x_lo (0.03–0.075) is not counted. There the origin truncation produces a spurious
crossing at μ ≈ 0.6/|s_start| that moves when the truncation is moved (§4).

## 2. Real versus complex
- **2D, n ≤ 7:** all unstable eigenvalues are real. The real-axis scans show no complex ν-pair of modulus > 1 for
  μ ≥ 0.045, and each contour count equals the real count plus one.
- **Hou–Luo, n ≤ 7:** all real.
- **Hou–Luo, n ≥ 8:** weakly oscillatory complex pairs appear in the upper part of the unstable spectrum, where two
  adjacent real eigenvalues collide. The total index stays n.

| n | real unstable eigenvalues | complex pairs |
|---|---|---|
| 8 | 0.79095, 0.44373, 0.35277, 0.25583, 0.16004, 0.06372 | 0.587182 ± 0.022206i |
| 9 | 0.80084, 0.46930, 0.40407, 0.31343, 0.22823, 0.14278, 0.05681 | 0.606624 ± 0.056719i |
| 10 | 0.81290, 0.35813, 0.28329, 0.20598, 0.12891, 0.05081 | 0.623297 ± 0.083614i, 0.463219 ± 0.020897i |

**Method note.** On the real axis, det(I − T_μ) has the sign (−1)^{#real ν>1}. An eigenvalue μ therefore flips
the parity of that count. A collision of two real ν > 1 into a complex pair does not. `hl_stab_scan.py` counts
parity flips, and `hl_find_complex.py` locates the complex pairs by secant iteration on ν(μ) = 1.

## 3. The ladder μ_{n,k} ≈ (λ_n − 1)(k + c)
Sorted from the smallest (k = 0), μ_{n,k}/(λ_n − 1), with extrapolations to λ_n → 1 (`ladder_fit.py`):

| | k = 0 | k = 1 | k = 2 | spacing k=0→1 | spacing k=1→2 |
|---|---|---|---|---|---|
| 2D, n = 7 | 0.7908 | 1.8880 | 2.9852 | 1.0972 | 1.0973 |
| 2D, limit | 0.71–0.74 | 1.74–1.76 | 2.70–2.72 | 1.01–1.03 | 0.94–0.98 |
| Hou–Luo, n = 8 | 0.7046 | 1.7696 | 2.8287 | 1.0650 | 1.0591 |
| Hou–Luo, limit | 0.66 | 1.67–1.68 | 2.67–2.74 | 1.01–1.02 | 0.99–1.06 |

**Result.** In both models the lower unstable eigenvalues form an arithmetic ladder. Its spacing tends to
(λ_n − 1)(1.00 ± 0.03), and its offset c is model-dependent: ≈ 0.72 in 2D and ≈ 0.66 in Hou–Luo. The largest unstable
eigenvalue approaches ≈ 0.78–0.8 in both models (2D: 0.746 at n = 7; Hou–Luo: 0.813 at n = 10).

**A curiosity.** The single unstable eigenvalue of the first rung is 0.37379 in 2D and 0.37385 in Hou–Luo.

## 4. Localization of the unstable modes
- **Hou–Luo, n = 8** (`hl_modes.py`, `hl_modes_lam1.0904.png`). All six real modes peak sharply at the dip
  (η = −0.55) and the front (η = −0.47) of the stalled layer. Toward the stagnation point they decay like x²:
  d ln|θ'|/dη = 1.92–2.01 on η ∈ [−4.5, −2], i.e. they are regular.
- **2D, λ₄'s lowest mode** (`gstab_vec.py`, μ = 0.11908). θ' peaks at s = −0.55 (the dip), ω' at s = −0.42 (the
  front), and X' at s = 0.23.
- So the unstable modes live at the same front whose structure controls the profile ladder: the phase
  corrections, the dip scaling, and the limiting cusp.

## 5. An exact ladder identity of the transport operator
- The base temperature satisfies V·∇Θ̄ = (λ−1)Θ̄. Hence the temperature and vorticity transport operators
  L^θ_μ = (μ+1−λ) + V·∇ and L^ω_μ = (μ+1) + V·∇ obey

    L^θ_μ(Θ̄φ) = Θ̄ L^θ_{μ+λ−1} φ,     L^ω_μ(Θ̄φ) = Θ̄ L^ω_{μ+λ−1} φ,

  exactly, everywhere. Multiplication by Θ̄ intertwines the transport problem at μ and at μ + (λ − 1).
- The couplings of the linearized operator break this symmetry: buoyancy ∂₁θ', and the Biot–Savart terms
  u'·∇Θ̄ and u'·∇Ω̄.
- The observed ladder spacing λ_n − 1 is the value this symmetry would give. We have not derived why the lower
  unstable spectrum inherits it exactly in the limit, nor the offset c, nor why the count is n.

## 6. Spectral flow along the continuous branch (Hou–Luo) [N]
The branch profiles between the rungs can be followed continuously in λ (`hl_deep_spec.py`; parity-flip refinement as
in §2; origin truncation η₀ = −60/−100). This works because the states of scans E and F are smooth to |m − 2| ≲ 1e-6.
Each state is first re-converged at its exact λ (`reconverge.py`, λ from the scan logs).

**Re-convergence matters.** The branch-state file names carry λ rounded to six decimals. Used with the rounded λ, a
state is inconsistent: m shifts by Δλ/ε ≈ 1e-5 and the residual is ~1e-5. That leaves the upper ladder members
unchanged, but it displaces the lowest, front-localized member and produces spurious low roots, including a spurious
complex pair at z = 20.56. The first runs, done this way, are superseded and kept only as `hl_saw.out` and
`hl_deep_[A-D].out`.

**Definitions.**
- f is the fraction of the rung interval, f = (z − z_n)/(z_{n+1} − z_n), with z = 1/(λ−1).
- s is the local spacing of the lower ladder, in units of λ − 1.
- The offset is that of the lowest member in units of s. It is predicted as c̄ + f (mod 1), where c̄ = 0.663 is the
  mean offset at the rungs 7–10.

| z | interval | f | ν̂ = μ/(λ−1) | s | offset: observed | offset: c̄ + f |
|---|---|---|---|---|---|---|
| 10.310 | 7→8 | 0.410 | 1.148, 2.213 | 1.065 | 0.079 (+1) | 0.069 (+1) |
| 10.560 | 7→8 | 0.608 | (0.279), 1.357, 2.420, 3.489 | 1.066 | 0.273 | 0.267 |
| 10.710 | 7→8 | 0.726 | 0.409, 1.482, 2.543 | 1.067 | 0.383 | 0.385 |
| 11.110 | 8→9 | 0.042 | 0.749, 1.813, 2.872 | 1.061 | 0.706 | 0.701 |
| 11.510 | 8→9 | 0.357 | 1.085, 2.144, 3.200 | 1.058 | 0.026 (+1) | 0.016 (+1) |
| 11.910 | 8→9 | 0.673 | 0.347, 1.417, 2.472 | 1.063 | 0.326 | 0.331 |
| 12.310 | 8→9 | 0.988 | 0.687, 1.747, 2.800 | 1.057 | 0.650 | 0.647 |
| 12.710 | 9→10 | 0.303 | 1.022, 2.076, 3.128 | 1.053 | 0.970 | 0.962 |
| 13.110 | 9→10 | 0.618 | 0.285, 1.353, 2.404 | 1.060 | 0.269 | 0.277 |

"(+1)" means the lowest member present is one spacing above the predicted bottom slot, because that slot lies below the
scan floor ν̂ = 0.2. The value in parentheses at z = 10.560 comes from the first run; the re-converged scan stops at
ν̂ = 0.3.

**Result.**
- Between the rungs the whole lower ladder moves up rigidly by one spacing per rung interval, linearly in z, i.e.
  linearly in the WKB phase of the profile ladder.
- The observed offset agrees with c̄ + f within 0.011 at every point.
- **The entry is clean.** A new member appears at the bottom slot at μ = 0 when c̄ + f passes 1, i.e. at f ≈ 0.34 of
  each interval. It is then found exactly at its predicted position: 0.409/0.417, 0.347/0.363, 0.687/0.694 and
  0.285/0.286 (observed/predicted).
- So each rung has exactly one more unstable mode than the previous one.
- The local spacing decreases towards 1 (Fig. 3d): s = 1.06–1.07 at z ≈ 10.5–11, 1.053–1.060 at z ≈ 12.3–13.1. The
  deep states z = 13–38, re-converged, are running (`hl_deep_exact.out`).

## 7. Formal theory: one phase quantizes both ladders [F]
**(a) The WKB phase does not depend on the growth rate.** In the stalled layer the perturbations are
∝ exp(i∫κ ds/ε). Adding the growth rate μ to the local problem shifts the root by exactly iμ/D̂:
- In Hou–Luo this is an algebraic identity. The local relation is v(1 + v) = iΓ with v = μ + iD̂κ, so
  κ(μ) = κ(0) + iμ/D̂ for every μ.
- In 2D the same holds for the full local eigenproblem (vertical structure, shear ĉ, Θ_y). Solved by shooting, it
  gives κ(μ) − κ(0) = iμ/D̂ to 1e-10 for μ = 0.02–0.3 at every tested wall point.
- Hence Re Φ, and with it the quantization of the profile ladder, is the same for every μ. The μ-dependence of the
  WKB solutions is the real transport factor |Θ̄|^{−μ/(λ−1)}, the WKB form of the exact identity of §5.

**(b) The connection at the stagnation point supplies a phase linear in μ.**
- Near x = 0 the transport is ε[(ν − 2) + X∂_X], with X = x/ε and ν = μ/ε (m = 2, so ε = (λ−1)/2).
- In Hou–Luo the leading-order problem there is (ν − 2)θ + Xθ_X = −2X[H(θ) + strain]. Fourier transformation gives
  θ̂ ∝ |k − k₀|^{ν−3} near the WKB wavenumber k₀.
- Regularity at the stagnation point removes the branch |k| > k₀. The remaining endpoint singularity produces the
  wave Γ(ν − 2) X^{2−ν} exp(i[k₀X ± π(ν − 2)/2]).
- Its phase, ±π(ν − 2)/2 = ±π(μ/(λ−1) − 1), is the phase of the Mellin transform of the wave. It is linear in μ with
  slope π/(λ−1), whatever the details of the model.

**(c) Eigen-condition.** Matching the two gives

    Re Φ(λ) − π μ/(λ−1) + θ₀ = (k + ½)π + o(1).

Re Φ(λ) is the same phase that quantizes the profiles, Re Φ(λ_n) = nπ + δ. θ₀ collects the O(1) connection phases.

**(d) Consequences and their tests.**
| prediction | observed |
|---|---|
| lower spacing → λ − 1 exactly | 1.00 ± 0.03 (extrapolated from the rungs, both models); along the HL branch s = 1.07 → 1.035 for z = 10.5 → 18 |
| the offset follows the profile phase between rungs | §6: within 0.011 at 9 points (Hou–Luo) |
| the offset at the rungs is a model constant c = frac((δ + θ₀)/π − ½) | 0.65–0.67 (HL, n = 7–10); 0.70–0.72 (2D, n = 3–7) |
| one unstable eigenvalue per half-turn of the phase, entering at μ = 0 when the offset passes an integer; so index(n) = n + const, with const = 0 from any one rung | new members observed at their predicted positions (§6); counts n for n = 0–7 (2D) and 0–10 (HL) |
| the uniform ladder of n members ends near s·n(λ_n − 1) → C/π = 1/Δz_∞ | its computed top is 0.73–0.74 (HL, n = 8–10) and 0.65 (2D, n = 7), against C/π ≈ 0.79 and 0.66–0.68. The highest eigenvalue lies a further 0.06–0.10 above it (0.79–0.81 HL, 0.75 2D) and is not described by the leading-order condition |

**The upper end.**
- The top member sits 0.06–0.10 above its uniform-ladder position in both models, and this displacement does not
  shrink over the computed rungs.
- In Hou–Luo, adjacent upper members collide into the weakly complex pairs of §2. Each pair has the mean of the two
  ladder positions: 0.587 against 0.5865 at n = 8, and 0.623 against 0.631 at n = 10.

**(e) Not derived.**
- The constants θ₀ and c. They require the front connection and the full stagnation-point problem with the strain
  forcing.
- The 2D stagnation-point connection. The Mellin phase is model-independent, but the 2D region-I problem has not
  been solved.

## 8. Status
- **Established numerically:**
  - the counts in §1 and the real/complex structure in §2;
  - the ladder spacing in §3 and the localization in §4;
  - the phase-locked spectral flow in §6.
- **Exact:** the identity in §5, and κ(μ) = κ(0) + iμ/D̂ for the Hou–Luo local root.
- **Formal:** the eigen-condition of §7 and its consequences: unit spacing, phase locking, and index = number of
  half-turns.
- **Open:**
  - the offsets c;
  - a 2D between-rung test (running: `bq_flow_spec.py`);
  - whether 2D develops complex pairs at higher n, as Hou–Luo does.
