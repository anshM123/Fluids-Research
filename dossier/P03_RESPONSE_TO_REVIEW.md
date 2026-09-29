# P03 — point-by-point response to the internal review

Every item points to the file that holds the evidence. Numbers are from `programs/P03_boussinesq_ladder`.

## (A) Align the claims
| Request | Done |
|---|---|
| Do not claim proof | No document claims a proof. The one proved statement is an elementary lemma (ASYMPTOTICS.md §2). Every claim carries a tag: [N] numerical, [F] formal, [P] proved, [O] open. |
| Headline "Numerical evidence and asymptotic analysis reveal a quantized hierarchy …, with evidence for an infinite ladder accumulating at λ = 1" | Used verbatim in P03_DOSSIER.md, the top-level DOSSIER.md and the manuscript abstract. |
| Central concept: "discrete smooth blow-up profiles are regularity-selected points on a continuous singular branch" | P03_DOSSIER.md and manuscript §2. The manuscript title was changed to "Regularity-selected blow-up …". |
| The 8 core claims | P03_DOSSIER.md §1 lists them one by one, each with its status and evidence. |
| Call profiles 5–8 "next in continuation ordering" until their instability indices are verified | The indices are now verified. The march method gives 4, 5, 6, 7 for λ₄…λ₇. The independent method confirms every eigenvalue of λ₄–λ₆ and the six largest of λ₇. λ₇'s lowest (0.0699) lies in method 2's artifact band; it is confirmed by method 1 with two truncations. |
| 3/2 only as a candidate unless derived | 3/2 is a candidate, not derived. Our estimate is now more careful (see D). |

## (B) Stability spectra of at least λ₄, λ₅, λ₆
Done for all eight profiles. The counts come from method 1 (real scans plus argument-principle checks); the eigenvalue values come from both methods (P03_DOSSIER.md Table 3, `fig_bq_stability.png`).

**Methods.**
- Method 1 is the eigen-condition ν(μ) = 1 of the linearised march + Biot–Savart map T_μ: real-axis scans, Brent
  refinement, and argument-principle counts for complex μ.
- Method 2 is a shift-invert eigen-solve of the independent global discretisation.

**Instability indices.** 0, 1, 2, 3, 4, 5, 6, 7 for λ₀…λ₇: n unstable modes for the n-th profile, as Wang et al.
reported for n ≤ 3.

**Agreement between methods.** Method 2 reproduces all 21 unstable eigenvalues of λ₁–λ₆ and the six largest of λ₇
to ≤ 5e-5. Its own discretisation artifacts (a truncation-dependent dense family, inflow-boundary modes, a
pseudospectral cloud) are why it is not used for counting. Examples:
- λ₁: μ = 0.373789 (method 1) and 0.37379 (method 2);
- λ₂: 0.55418 / 0.22048 and 0.55419 / 0.22048.

**Checks.**
- The trivial time-translation mode μ = 1 is recovered by both methods.
- Truncation artifacts near μ ≈ 0.66/|s_start| move when the truncation is moved (0.052 → 0.033 → 0.022 for
  s_start = −12, −20, −30), and are discarded.
- Genuine eigenvalues do not move; for example λ₇'s lowest is 0.0699 for both s_start = −20 and −30.
- Complex modes: the argument-principle counts on [0.06, 1.5] × [−1.5, 1.5] give exactly real + trivial (λ₁: 2;
  λ₂: 3; λ₃: 4; λ₄: 5; λ₅: 5.97 ≈ 6, with one unresolved phase step at the left edge; λ₆, λ₇ not checked).

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
- λ₇: 9.6e-6 at the global solver's default origin truncation s_min = −12, and 2.1e-7 with s_min = −20.

The remaining λ₆ difference lies within the marching solver's own hs-uncertainty for that rung.

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

---------------------------------------------------------------------------------------------------------

# Second review ("push this to Nature"): point-by-point

**1. Tighten "exactly n unstable modes".**
- The text now says: "the n-th computed profile has n resolved real unstable modes (0 ≤ n ≤ 7)".
- Complete right-half-plane counts now exist for n = 5, 6, 7 (`stab_contour2.py`, k = 16, adaptive):
  5.952, 6.942, 7.921 zeros, i.e. n + 1 including the trivial mode. No oscillatory unstable spectrum was found.
- Domain variants (origin truncation −30, contour moved to x_lo = 0.075/0.07) give 5.954 (n = 5) and 6.945 (n = 6),
  against 5.952 and 6.942. The grid variants (hs = 0.0125) give 7.003 (n = 6) and 7.973 (n = 7).
- As you proposed, "exactly n" is therefore restored as a numerical statement: rung n has exactly n unstable
  eigenvalues, all real, for 0 ≤ n ≤ 7.
- A control with 11 rungs of the Hou–Luo model (`hl_stability.py`) gives exactly n unstable eigenvalues for
  n = 0…10, unchanged under a domain change. There, complex pairs appear from n = 8 while the total stays n. This
  is exactly why "real" and "exactly n" must be claimed separately.

**2. The word "infinite".** Every document now carries the hierarchy table:
| level | statement |
|---|---|
| Numerically established | 8 rungs |
| Strongly supported | the singular branch and regularity selection |
| Asymptotically predicted | unbounded phase and an infinite ladder |
| Not proved | infinitely many smooth profiles; C > 0 in 2D |

The abstract says explicitly that the infinite sequence is an asymptotic prediction. New Hou–Luo data down to
ε = 0.012 (λ − 1 ≈ 0.025) support the leading-order form: the phase coefficient converges, with a non-analytic correction
∝ (λ−1)^{2/3}, so Φ = C/(λ−1) + O((λ−1)^{−1/3}) there.

**3. Independent solver front and center.** It is now §2 of the manuscript, with both solvers side by side in
Table 1.

**4. Title, abstract, evidence chain.** Adopted: "Quantization and accumulation of unstable self-similar blow-up
profiles in the 2D Boussinesq equations". The abstract leads with the discovery, then the redundancy, then the
mechanism.

**5. The μ/ε limiting spectrum** (priority 2; INSTABILITY_LADDER.md §3–7, `fig_instability_ladder.png`,
`fig_spectral_flow.png`). This is now the main new result.
- **The ladder.** In both models the lower unstable eigenvalues form an arithmetic ladder
  μ_{n,k} ≈ (λ_n − 1)(k + c). Its spacing extrapolates to 0.996 times λ_n − 1 in both models. Its offset at
  the rungs is a model constant: 0.65–0.67 (Hou–Luo), 0.70–0.72 (2D).
- **Spectral flow along the continuous branch.** We followed the spectrum through nine Hou–Luo branch points between
  the rungs 7 and 10, and five 2D Boussinesq branch points between the rungs 3 and 5. Each state was re-converged at
  its exact λ.
  - The whole lower ladder moves up rigidly by one spacing per rung interval, linearly in z = 1/(λ−1).
  - Its offset equals the rung value plus the fraction of the interval, to within 0.011 (Hou–Luo) and 0.03 (2D).
  - One new member enters through μ = 0 per interval, at the predicted phase and position. That is why each rung has
    one more unstable mode than the last.
- **A formal eigen-condition from the same phase.** Two ingredients:
  - The growth rate enters the stalled-layer WKB root only as κ(μ) = κ(0) + iμ/D̂. This is exact for Hou–Luo and
    verified to 1e-10 for the 2D local eigenproblem. So the profile-ladder phase Re Φ is the same for every μ.
  - The stagnation point, where the transport is an Euler operator, contributes the Mellin/Γ-function phase
    π(μ/(λ−1) − 1).

  Together: Re Φ(λ) − πμ/(λ−1) + θ₀ = (k + ½)π. This predicts:
  - a spacing of exactly λ − 1;
  - the phase-locked flow;
  - index = number of half-turns of Φ, i.e. n.

  All three agree with the data. The highest eigenvalue lies 0.06–0.10 above the uniform ladder, which the
  leading-order condition does not describe.
- **Your conjecture.** "One singular limiting operator generates both the profile ladder and its instability ladder"
  thus becomes a concrete, tested statement: one WKB phase quantizes both.
- **Not derived.** The offsets θ₀ and c, and the 2D stagnation-point problem itself. The 2D spectral flow and spacing,
  however, are those the eigen-condition predicts.

**6. The front/dip inner problem** (priority 3).
- Measured in Hou–Luo to ε = 0.012: depth ∝ ε^{0.31}, width ∝ ε^{0.34}, and the dip-to-front distance ∝ ε^{1.23}.
- The 2D exponents measured at ε ≥ 0.044 are pre-asymptotic.
- An analytic inner solution is not yet available; the asymptotic 2D spacing stays a range, 1.476–1.51.

**7. Validated-numerics certificate** (priority 4). Not attempted. It is listed as the natural next project: the
global solver's exact sparse Jacobian and polynomial (quadratic) nonlinearity make a Newton–Kantorovich
certificate for the discretized problem straightforward. A certificate for the PDE needs rigorous
discretization-error bounds.

**Literature.** The reviewer's check through September 2026 has been added to AUDIT_I:
- the follow-up arXiv:2511.22819 (improved accuracy and an additional IPM unstable solution) has no Boussinesq
  rungs beyond λ₃ or λ₄;
- Chen–Huang–Li (2026) treats a different class of singular two-stage profiles;
- the stable profile has computer-assisted backing, while the unstable family is numerical only.
