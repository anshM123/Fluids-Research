# Regularity-selected blow-up: a quantized hierarchy of self-similar singularities in the 2D Boussinesq equations

*Draft. Every number comes from programs/P03_boussinesq_ladder, and each value has a log file there. Claims are
tagged in `dossier/P03_DOSSIER.md`: numerical, formal-asymptotic, proved, or open.*

## Abstract
Unstable self-similar singularities of the two-dimensional Boussinesq equations with a boundary, the standard proxy
for axisymmetric three-dimensional Euler blow-up, were recently discovered with physics-informed neural networks.
Four smooth profiles were found, with blow-up rates λ_n following an empirical law.

We show that these profiles are not isolated objects. They are regularity-selected points on a single continuous
branch of least-singular self-similar solutions. Along this branch the local exponent m(λ) at the stagnation point
equals 2 exactly at a smooth profile.

**Numerics.**
- Continuation along the branch with a classical solver locates eight such regularity resonances, from
  λ₀ = 1.9205593 to λ₇ = 1.0883384. Four of them lie beyond those previously reported.
- An independent global Newton solver reproduces all eight resonances to ≤ 1e-5 (≤ 6e-7 for n = 1–5).
- Linear stability finds n unstable modes for the n-th profile for all eight, which extends the pattern reported for
  the first four.
- A second, independent method reproduces all 21 unstable eigenvalues of λ₁–λ₆ to about 1e-5.

**Mechanism.**
- As λ → 1⁺ the profiles develop a stalled boundary layer next to a front that tends to a square-root cusp.
- In this layer the smoothness defect is carried by exponentially small oscillatory solutions of a boundary-layer
  eigenproblem. The resonances therefore obey a phase-quantization condition Re Φ(λ_n) = nπ + δ + o(1), with
  Φ ∼ C/(λ−1) and C > 0.
- The computed phase gain per resonance approaches π (3.12 at n = 7). This is evidence for an infinite ladder
  accumulating at λ = 1, with 1/(λ_n − 1) ≈ (π/C) n. Here π/C lies between 1.476 and about 1.51, and 3/2 is a
  candidate.

**Contrast.** The Hou–Luo boundary model shows the same mechanism. In the Córdoba–Córdoba–Fontelos model, by
contrast, the cusp forms at a finite λ* ≈ 0.4536 and the ladder ends after three profiles.

## 1. Introduction (outline)
**Background.**
- Blow-up for 3D Euler with boundary: the Luo–Hou scenario; stable blow-up (Chen–Hou); the role of unstable
  singularities for the Navier–Stokes question (Wang et al. 2025).
- The empirical ladders of Wang et al.: are they finite or infinite, and what sets their law?

**Our viewpoint.** Smooth profiles are the regularity-selected members of a continuous family of singular
profiles. The search for smooth profiles, one at a time, becomes root finding along one branch; the branch is then
analysed asymptotically.

**Analogy.** Phase quantization of a WKB wave in a stalled layer plays the role of a Bohr–Sommerfeld condition
for blow-up rates. The analogy is used for intuition only.

## 2. The least-singular branch and its regularity resonances
- **Local analysis.** At the stagnation point Θ ≈ −|y₁|^m with m = (λ−1)/(1+λ−A). Smooth ⇔ m = 2 ⇔
  λ = −3 − 2∂₁U₁(0).
- **Solvers.**
  - (i) A marching solver: log-polar, exact local structure factored out, Newton–Krylov.
  - (ii) An independent global Newton solver: different unknowns, finite differences, and the elliptic
    Biot–Savart law; smoothness imposed and λ treated as an eigenvalue.
  - Methods and resolution tables are in the SI.
- **Fig. 1a.** m(λ) − 2 against z = 1/(λ−1): an oscillation about zero with geometrically decaying amplitude. The
  eight zeros are the smooth profiles (Table 1). The two solvers agree (Table 2).
- **Fig. 2.** Stability. Method (a) uses the eigen-condition ν(μ) = 1 for the linearised march–Biot–Savart map
  T_μ: real-axis crossings, plus an argument-principle count that includes complex μ. Method (b) is a shift-invert
  eigen-solve of the global linearisation. Both give n unstable modes for profile n (Table 3).

## 3. Anatomy of the λ → 1 limit: a stalled layer next to a limiting cusp
- **Fig. 1c.** Along the boundary, the radial self-similar speed is O(λ−1) up to a front x_c ≈ 0.72, then O(1).
- **Cusp.** Outside the front D ∝ (x − x_c)^{1/2}. This is the square-root cusp that terminates the CCF branch,
  reached here only as λ → 1. The dip in D closes only as ε → 0 (D̂_min ∝ ε^{0.18}).
- **Why m − 2 is exponentially small.** A non-smooth profile induces a wall velocity ∝ (m−2), which is
  incompatible with the O(ε) stalled flow at every order of the formal expansion.

## 4. Phase quantization (ASYMPTOTICS.md)
**Local problem.** Perturbations exp(i∫κ ds/ε) with vertical scale εx satisfy a transport–Biot–Savart
eigenproblem; all dropped terms are O(ε). The Hou–Luo analogue has the closed form iD̂κ² + κ − Ω = 0, and Re κ > 0
whenever Ω > 0.

**Representation.**
- m − 2 = |K| e^{−Im Φ}[cos(Re Φ + arg K) + O(ε)].
- Re Φ = C/(λ−1) + E, with C = 2 Re∫κ₀ ds > 0.
- E = o(1/(λ−1)) in general, and E = Φ₀ + o(1) under a stated uniformity hypothesis.

**Elementary lemma (proved).** With this representation:
- the resonances are infinitely many;
- they satisfy Re Φ(λ_n) = nπ + δ + o(1);
- they obey λ_n = 1 + C/(nπ + c + o(1)) whenever E = Φ₀ + o(1).

**Numerical verification.**
- Phase gain per resonance 2.98, 3.03, 3.07, 3.07, 3.115, 3.12 (2D) and 3.02–3.09 (Hou–Luo).
- The linear law Re Φ = Cz + Φ₀ holds to ±0.024 over n = 4–7 for every front cut-off.
- Integrated through the front, Re Φ(λ_n) − nπ is constant to ±0.03.
- The directly measured spacings increase monotonically, 1.4646 → 1.4758, and bound π/C from below.
- The extrapolated spacing depends on the correction exponent, which seven rungs cannot fix:
  - analytic 1/z phase corrections give π/C = 1.478;
  - the non-analytic correction suggested by the measured dip scaling (∝ z^{0.6}) gives 1.49–1.50;
  - the finite-z WKB slopes give 1.478–1.493.
- 3/2 (C = 2π/3) is a candidate consistent with the non-analytic extrapolation. It is not derived.

**Heuristic (not a derivation).** The n-th resonance carries n half-wavelengths of the layer oscillation, and its
instability index is n. No link between the two is derived here. The lower unstable eigenvalues scale with
ε = (λ_n−1)/2 and form an almost equally spaced sequence, μ_k/ε ≈ 1.6, 3.8, 6.05, 8.4, …

## 5. Finite versus infinite ladders
- **CCF (P02).** The sonic depth vanishes at λ* = 0.4535843. The branch ends in a square-root cusp with a
  log-periodic approach and carries three profiles. This is consistent with the unsuccessful search for a third
  profile in λ ∈ [0.455, 0.4713] (Wang, Léger, Lai, Buckmaster, arXiv:2511.22819).
- **Boussinesq and Hou–Luo.** The cusp is approached only as λ → 1, and the stalled layer in front of it supports
  the quantized oscillation. The evidence therefore points to an unbounded ladder.

## 6. Methods and verification (SI)
- **Solver details and pitfalls.** FFT/Gibbs contamination of Newton directions; roundoff-limited Jacobian-vector
  products.
- **Resolution tables.**
  - Marching solver: Nb 32/48/64; hs 0.025/0.0125; s_start −12/−20; s_max 100/130.
  - Global solver: hs, Nb, s_min, s_max with extrapolation.
  - Hou–Luo: N = 8192–131072.
- **Stability.**
  - Validation: |T₁v − v|/|v| = 3e-6 on the time-translation mode.
  - The truncation-induced modes near μ ≈ 0.66/|s_start| move with the truncation and are discarded.
  - The argument-principle count on [0.06, 1.5] × [−1.5, 1.5] returns n + 1 zeros (trivial mode included) ⟨for the
    profiles checked⟩.
- **Limitations.**
  - Numerical evidence plus formal asymptotics, not a proof.
  - The deepest resonances rely on |m − 2| ~ 1e-7.
  - The published λ₂, λ₃ of Wang et al. could not be accessed for a digit-by-digit comparison.

## Table 1 — resonances (marching solver)
| n | λ_n | 1/(λ_n − 1) |
|---|---|---|
| 0 | 1.9205593 | 1.08630 |
| 1 | 1.3990961 | 2.50566 |
| 2 | 1.2523487 | 3.96277 |
| 3 | 1.1842533 | 5.42733 |
| 4 | 1.1449857 | 6.89723 |
| 5 | 1.1194738 | 8.37004 |
| 6 | 1.1015817 | 9.84429 |
| 7 | 1.0883384 | 11.32010 |

## Table 2 — independent global solver (s_max-extrapolated; hs = 0.0125 for n ≥ 3)
| n | marching | global | difference |
|---|---|---|---|
| 0 | 1.9205593 | 1.9205610 | +1.7e-6 |
| 1 | 1.3990961 | 1.3990960 | −1.1e-7 |
| 2 | 1.2523487 | 1.2523481 | −5.6e-7 |
| 3 | 1.1842533 | 1.1842530 | −3.2e-7 |
| 4 | 1.1449857 | 1.1449853 | −4.2e-7 |
| 5 | 1.1194738 | 1.1194739 | +1.1e-7 |
| 6 | 1.1015817 | 1.1015771 | −4.6e-6 |
| 7 | 1.0883384 | 1.0883480 | +9.6e-6 |

## Table 3 — unstable eigenvalues μ_k (perturbations ∝ e^{μτ})
Method 2 (global eigen-solver) is shown. Method 1 (march-based) brackets each value within 0.05 and refines
λ₁, λ₂ to 1e-5.

| n | index | μ_k |
|---|---|---|
| 0 | 0 | — |
| 1 | 1 | 0.37379 |
| 2 | 2 | 0.55419, 0.22048 |
| 3 | 3 | 0.63437, 0.37541, 0.15492 |
| 4 | 4 | 0.68004, 0.45973, 0.28629, 0.11908 |
| 5 | 5 | 0.70985, 0.51143, 0.36869, 0.23121, 0.09654 |
| 6 | 6 | 0.73139, 0.54481, 0.42577, 0.30749, 0.19415, 0.08154 |
| 7 | 7 | 0.74632, 0.56676, 0.46713, 0.36192, 0.26371, 0.16678, 0.0699 (method 1) |
