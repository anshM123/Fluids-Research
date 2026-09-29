# P03 — point-by-point response to the internal review

Every item points to the file that holds the evidence. Numbers are from `programs/P03_boussinesq_ladder`.

## (A) Align the claims
| Request | Done |
|---|---|
| Do not claim proof | No document claims a proof. The one proved statement is an elementary lemma (ASYMPTOTICS.md §2). Every claim carries a tag: [N] numerical, [F] formal, [P] proved, [O] open. |
| Headline "Numerical evidence and asymptotic analysis reveal a quantized hierarchy …, with evidence for an infinite ladder accumulating at λ = 1" | Used verbatim in P03_DOSSIER.md, the top-level DOSSIER.md and the manuscript abstract. |
| Central concept: "discrete smooth blow-up profiles are regularity-selected points on a continuous singular branch" | P03_DOSSIER.md and manuscript §2. The manuscript title was changed to "Regularity-selected blow-up …". |
| The 8 core claims | P03_DOSSIER.md §1 lists them one by one, each with its status and evidence. |
| Call profiles 5–8 "next in continuation ordering" until their instability indices are verified | The indices are now verified. The march method gives 4, 5, 6, 7 for λ₄…λ₇. The independent method confirms every eigenvalue of λ₄–λ₆ and the six largest of λ₇; λ₇'s lowest is confirmed by method 1 with two truncations, and the method-2 check at s_min = −20 is ⟨pending⟩. |
| 3/2 only as a candidate unless derived | 3/2 is a candidate, not derived. Our estimate is now more careful (see D). |

## (B) Stability spectra of at least λ₄, λ₅, λ₆
Done for all eight profiles, with two independent methods (P03_DOSSIER.md Table 3, `fig_bq_stability.png`).

**Methods.**
- Method 1 is the eigen-condition ν(μ) = 1 of the linearised march + Biot–Savart map T_μ: real-axis scans, Brent
  refinement, and argument-principle counts for complex μ.
- Method 2 is a shift-invert eigen-solve of the independent global discretisation.

**Instability indices.** 0, 1, 2, 3, 4, 5, 6, 7 for λ₀…λ₇: n unstable modes for the n-th profile, as Wang et al.
reported for n ≤ 3.

**Agreement between methods.** Method 2 reproduces all 21 unstable eigenvalues of λ₁–λ₆ and the six largest of λ₇
to about 1e-5. Examples:
- λ₁: μ = 0.373789 (method 1) and 0.37379 (method 2);
- λ₂: 0.55418 / 0.22048 and 0.55419 / 0.22048.

**Checks.**
- The trivial time-translation mode μ = 1 is recovered by both methods.
- Truncation artifacts near μ ≈ 0.66/|s_start| move when the truncation is moved (0.052 → 0.033 → 0.022 for
  s_start = −12, −20, −30), and are discarded.
- Genuine eigenvalues do not move; for example λ₇'s lowest is 0.0699 for both s_start = −20 and −30.
- Complex modes: the argument-principle counts on [0.06, 1.5] × [−1.5, 1.5] give exactly real + trivial (λ₁: 2;
  λ₂: 3; ⟨λ₃–λ₅ pending⟩).

**New observation.** The lower unstable eigenvalues scale with ε = (λ_n − 1)/2: μ_k/ε ≈ 1.6, 3.8, 6.05, 8.4, …

## (C) Independent reproduction of λ₄–λ₆ with a different representation
Done for all eight (P03_DOSSIER.md Table 2).

**The global Newton solver** (`bq_global.py`) shares no numerical ingredient with the marching solver:
- unhatted Θ, Ω, X solved simultaneously;
- upwind-biased finite differences in s;
- the Biot–Savart law as a local elliptic equation, with no FFT;
- smoothness imposed through the local data, with λ as an eigenvalue;
- Newton with the exact sparse Jacobian;
- s_max extrapolation.

**Agreement with the marching solver.**
- λ₁–λ₅: 1e-7 to 6e-7.
- λ₀: 1.7e-6.
- λ₆: 4.6e-6.
- λ₇: 9.6e-6.

The λ₆ and λ₇ differences lie within the marching solver's own hs-uncertainty for these rungs.

## (D) Matched-asymptotic calculation with explicit error terms
See `programs/P03_boussinesq_ladder/ASYMPTOTICS.md`.

**Lemma (proved).** Write F = R(cos Θ + η) with R > 0, Θ ↑ ∞ and |η| < 1. Then there are infinitely many
resonances, and they satisfy Θ(λ_n) = (n + ½)π + o(1). If Θ = Cz + Θ₀ + o(1), then λ_n = 1 + C/(nπ + c + o(1)).

**Formal WKB derivation.**
- Eikonal: the local boundary-layer eigenproblem. Every dropped term is listed, with its order O(ε).
- Connection regions: the stagnation point (x ≲ ε) and the front.
- C = 2 Re ∫κ₀ ds > 0. Positivity is exact in closed form for Hou–Luo, where Ω > 0 holds throughout the layer; in
  2D it is checked pointwise (Re κ > 0 and Im κ < 0 on all profiles).

**The remainder, stated only as far as supported.**
- Hypothesis (U), that the stalled layer has a regular ε-expansion, is **tested**. It holds with O(ε) corrections
  for x ≲ 0.4 and fails in the dip region x ≈ 0.45–0.72.
- The dip scales as follows: depth ε^{0.18}, width ε^{0.69}, local root κ ~ D̂^{−3/2}. This gives a non-analytic
  phase correction ∝ z^{0.58}.
- So the supported form is Φ(λ) = C/(λ−1) + O((λ−1)^{−0.6}), not C/(λ−1) + Φ₀ + o(1). That is still sufficient for
  infinitely many resonances and for the limiting spacing π/C.
- Quantization Φ(λ_n) = nπ + δ + o(1): with the phase integrated through the front, Re Φ(λ_n) − nπ is constant to
  ±0.03 for n = 4–7.

**The spacing coefficient.** It is not derived analytically. Its extrapolation depends on the (unknown) correction
exponent, which seven rungs cannot fix:
- 1.478 if the corrections are analytic;
- 1.49–1.50 with the non-analytic correction implied by the dip scaling.

The honest statement is: π/C lies between 1.476 (the monotone lower bound) and ≈ 1.51, and 3/2 is a candidate.
The earlier "≈ 1.50" came from a single front cut-off (the cut-off study is in ASYMPTOTICS.md §5.2).

## What would settle the remaining open points
- **Asymptotic spacing.** Solve an inner problem for the dip/front region as ε → 0, or compute 3–4 more rungs at
  higher resolution (|m − 2| ≈ 1e-8 near n = 8).
- **The link between the rung number and the instability index.** It is observed, not derived. The ε-scaling of
  the lower eigenvalues suggests a WKB theory for the unstable spectrum too.
- **Literature re-check.** Re-check against Wang et al.'s tabulated values once arXiv is reachable.
