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
| 5 | B: 0.025, −30, [0.075,1.5]×[−1.5,1.5] | ⟨running⟩ | 6 |
| 6 | A: 0.025, −20, [0.06,1.5]×[−1.5,1.5] | 6.942 | 7 |
| 6 | B: 0.025, −30, [0.07,1.5]×[−1.5,1.5] | ⟨running⟩ | 7 |
| 6 | C: 0.0125, −20, [0.065,1.2]×[−0.8,0.8] | ⟨running⟩ | 7 |
| 7 | A: 0.025, −30, [0.05,1.5]×[−1.5,1.5] | 7.921 | 8 |

Deviations from an integer (≤ 0.08) come from eigenvalue swaps at the truncation of the top-k set. At the
unresolved steps the smallest |1 − ν| is ≥ 0.12, so no zero lies near the contour.

**Hou–Luo** (`hl_contour.py`, k = 20, [0.03,1.5]×[−1.5,1.5], N = 32768, η₀ = −30; n ≤ 2 at N = 8192):
| n | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| zeros | 1.991 | 3.000 | 4.000 | 5.000 | 6.000 | 7.000 | 8.000 | 9.000 | 9.959 | 11.000 |

Variants: η₀ = −40 with x_lo = 0.04 gives 6.000 (n = 5) and 7.024 (n = 6); ⟨rest running⟩. N = 65536: ⟨running⟩.

**Result.** In the rectangle, the n-th profile has exactly n unstable eigenvalues: n ≤ 7 in 2D and n ≤ 10 in
Hou–Luo. The strip 0 < Re μ < x_lo (0.03–0.075) is not counted. There the origin truncation produces a spurious
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

## 6. Status
- **Established numerically:** the counts in §1, the real/complex structure in §2, the ladder spacing in §3, and
  the localization in §4.
- **Exact:** the identity in §5.
- **Open:**
  - a derivation of the ladder (the limiting slow operator at the front);
  - the offsets c;
  - the mechanism fixing the index to n;
  - whether 2D, like Hou–Luo, develops complex pairs at higher n.
