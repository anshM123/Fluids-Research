# Validation status of the new profiles and their indices (task C, scoped)

This file states what is validated, and how, for the profiles and instability indices reported in the paper. It
also sets out what a computer-assisted proof (CAP) of one new rung would require. No CAP has been carried out.

## 1. What is validated numerically

**Existence and position of the smooth profiles (rungs).**
- **Solver A** (log-polar, marching Newton–Krylov).
  - Newton converges quadratically, to residual 1e-13.
  - The rung equation m(λ) = 2 is solved to |m − 2| < 1e-12.
- **Solver B** (global sparse Newton, `bq_global.py`) shares no numerical ingredient with solver A: grid, velocity
  evaluation, march and linear algebra are all different. It reproduces all eight 2D rungs to ≤ 1e-5 (≤ 1e-6 for
  five of them).
- **Grid convergence.** Refinement (h_s 0.025 → 0.0125; Nb 32 → 48, 64; s_start −20 → −12, −30; s_max 100 → 130)
  moves the rungs as follows:

  | model | rung | shift |
  |---|---|---|
  | 2D Boussinesq | n ≤ 4 | ≤ 1e-6 in λ |
  | 2D Boussinesq | n = 7 | 7e-4 in λ |
  | IPM | λ₁ | 9e-9 (h_s 0.025 → 0.0125) |

  Details are in README "Results" and `cross_l*` / `ipm_rung1_hs0125.out`.
- **Hou–Luo.** Two grid sizes (N = 65,536 and 131,072) and two domain sizes.

**The index.**
- **Real modes.** Parity changes of the number of real eigenvalues ν > 1 of T_μ, refined by bisection. They are
  reproduced by an independent linearization (method 2, global Jacobian) in 2D.
- **All modes (argument principle).** det(I − T_μ) is evaluated on [x_lo, 1.5] × [−1.5, 1.5], with adaptive
  refinement to |Δarg| < π/6.

  | model / rungs | count | stable under |
  |---|---|---|
  | 2D, n = 5–7 | exactly n + 1 zeros (index n plus the trivial μ = 1) | changes of origin truncation, contour and grid |
  | Hou–Luo, n ≤ 10 | index n | — |
  | IPM, n = 1–4 | 1.947, 3.082, 4.077, 4.902 → n + 1 (index n plus the trivial mode) | box height ±1.5 (n ≤ 2) and ±2.5 (n = 3, 4) |

  The deviation from the integer is ≤ 0.08 in all counts except IPM U₄ (0.098). There, two adjacent unresolved
  jumps of opposite sign (+1.017 and −1.031 rad at Im μ = 2.5, Re μ ≈ 0.161) cancel; min|1 − ν| = 0.7–0.9 at both, so
  they are an eigenvalue swap at the truncation, not a zero near the contour.
- **Large |Im μ|.** Along Re μ = x_lo the spectral radius of T_μ decays roughly like 2/|Im μ|:

  | case | ρ(T_μ) |
  |---|---|
  | IPM U₁ | 1.21, 0.84, 0.61, 0.48, 0.33, 0.25, 0.17, 0.12 at Im μ = 1.5, 2, 3, 4, 6, 8, 12, 24 |
  | 2D n = 7 | < 1 for Im μ ≥ 3 at every sampled point (max 0.84 near Im μ = 6), but 1.27 at Im μ = 2.5 |

  So eigenvalues with |Im μ| ≥ 3 are excluded wherever ρ < 1. The strip 1.5 ≤ |Im μ| ≤ 3.5 is closed by an
  argument-principle count on the rectangle [x_lo, 1.5] × [1.5, 3.5] (`contour_rect.py`). For 2D n = 7 the count
  is **0.0000**, i.e. no eigenvalues in the strip (`crect_bq7_strip.out`, 72 evaluations; one segment with an
  eigenvalue swap at the truncation, min|1 − ν| = 0.5 there). IPM U₃ and U₄ are counted directly on
  [0.04, 1.5] × [±2.5]: U₃ gives 4.077 → trivial + 3 and U₄ gives 4.902 → trivial + 4.
- **Small Re μ.** 0 < Re μ < x_lo contains only the dilation mode, which is split by the origin truncation to
  ±0.6/|s_start|. Its position moves as predicted when s_start changes (−20 → −30).

## 2. What is not validated
- **Rigour.** Nothing above is interval arithmetic. The statements are numerical, with convergence evidence.
- **The asymptotic statements.** The infinite hierarchy, the eigen-condition and the derived spacing are formal
  (matched asymptotics). Every hypothesis is checked numerically.

## 3. What a computer-assisted proof of one new rung would require
The natural target is the Hou–Luo rung n = 1 (λ₁ = 1.44767467): 1D, smooth, with one unstable direction.

**Framework.**
- A Cayley map x = tan(θ/2) makes the Hilbert transform diagonal in Fourier modes, and turns x d/dx into
  sin θ d/dθ. The self-similar equations Ω + VΩ' = Θ', VΘ' = (λ−1)Θ, U' = HΩ then become polynomial equations for
  Fourier coefficients.
- λ is an unknown, fixed by the smoothness (m = 2) condition. Equivalently, Θ is even and analytic at the origin.

**The obstacle.** The profile is not analytic at infinity (θ = π). Ω has the far-field exponents 1/(1+λ) and
2/(1+λ), and these depend on the unknown λ.
- The singular part must be split off, Ω = Σ_j a_j (1+x²)^{−p_j(λ)/2} + (smooth), with the p_j(λ) handled as
  functions of λ.
- Alternatively, use weighted Hölder/Sobolev estimates, as in the Chen–Hou proof for the stable profile.

**The certificate.** A radii-polynomial / Newton–Kantorovich bound in a weighted ℓ¹ space for the smooth part, with
explicit tail bounds, plus rigorous enclosures of the singular coefficients. python-flint (arb ball arithmetic) is
available in this environment.

**Effort.** Weeks to months of work. The approximate profile (accurate to 1e-10) and its linearization are
available here as starting data.

**The index.** A rigorous index needs validated enclosures of the unstable spectrum of the linearization. For
example: a rigorous argument-principle count for det(I − T_μ) on a contour, with interval evaluation of a
finite-rank approximation of T_μ, plus a bound on the remainder.

This is stated as future work in the paper (Discussion). The paper claims no proof.
