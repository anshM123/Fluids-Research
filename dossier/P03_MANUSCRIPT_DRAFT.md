# An infinite ladder of unstable self-similar blow-ups in the 2D Boussinesq equations, and its WKB origin

*Draft. The numbers come from programs/P03_boussinesq_ladder, and every value has a log file there.*

## Abstract
Unstable self-similar singularities of the 2D Boussinesq equations with boundary are the proxy for 3D
axisymmetric Euler blow-up. They were recently discovered with high-precision physics-informed neural networks,
which found four smooth profiles whose blow-up rates λ_n follow an empirical law 1/(λ_n − 1) ≈ an + b.

Here we show with a classical, neural-network-free continuation method that:
- all these profiles lie on a single continuous branch of least-singular self-similar solutions, and are the
  zeros of one scalar smoothness function;
- the branch carries at least eight smooth profiles (λ₀ = 1.9205593 … λ₇ = 1.0883384), four more than
  previously known.

As λ → 1, the profiles develop a quasi-stagnant boundary layer, bounded by a front that tends to a
square-root cusp. Inside the layer, the smoothness condition becomes a WKB quantisation condition for
exponentially small, oscillating perturbations, with a complex local wavenumber given by a boundary-layer
eigenproblem. The computed WKB phase gains π between consecutive smooth profiles. The ladder is therefore
infinite and accumulates at λ = 1, with an asymptotic spacing of 1/(λ_n − 1) equal to π/(2 Re a) ≈ 1.50 rather
than the fitted 1.42.

The Hou–Luo boundary model has the same structure, which we follow through twelve profiles. By contrast, the
Córdoba–Córdoba–Fontelos ladder is finite, because there the cusp forms at finite λ. This gives a criterion for
when unstable-singularity ladders terminate.

## 1. Introduction (outline)
- Blow-up for 3D Euler with boundary: Luo–Hou scenario; Chen–Hou stable blow-up; the importance of unstable
  singularities for the Navier–Stokes question (Wang et al. 2025).
- The empirical ladders of Wang et al. (IPM, Boussinesq): is the ladder infinite? What sets its law?
- Our approach: turn the discrete search for smooth profiles into root finding along one branch, then analyse
  the branch asymptotically. The same approach gave the finite CCF ladder (companion result, P02).

## 2. The least-singular branch
- Self-similar equations; local analysis at the stagnation point gives Θ ≈ −|y₁|^m with m = (λ−1)/(1+λ−A).
  Smooth profiles ⇔ m = 2 ⇔ λ = −3 − 2∂₁U₁(0).
- Log-polar solver with the exact local structure factored out, Newton–Krylov, continuation in λ. Resolution
  studies are in Table S1.
- Fig. 1a: m(λ) − 2 versus z = 1/(λ−1). Oscillation about zero with geometrically decaying amplitude; the
  eight zeros are the smooth profiles (Table 1).

## 3. Anatomy of the λ → 1 limit
- Fig. 1c: the radial self-similar speed along the boundary is O(ε) up to a front x_c ≈ 0.72, then O(1).
- Outside the front D ∝ (x − x_c)^{1/2}. This is the square-root cusp that terminates the CCF branch, reached here
  only as λ → 1. The dip in D closes only as ε → 0 (D̂_min ∝ ε^{0.18}).
- Why m − 2 is beyond all orders in ε: non-smooth velocity components are incompatible with the O(ε) stagnant
  flow.

## 4. WKB quantisation
- Local problem: perturbations exp(i∫κ ds/ε) with vertical scale εx satisfy a transport–Biot–Savart
  eigenproblem. The Hou–Luo analogue has the closed form iD̂κ² + κ − Ω/D̂ = 0.
- Quantisation: Re Φ(λ_n) = Φ₀ + nπ with Φ = ε⁻¹∫κ ds. Verified: ΔRe Φ = 2.98, 3.03, 3.07, 3.07, 3.115, 3.12
  (2D) and 3.02–3.09 (Hou–Luo).
- Consequences:
  - 1/(λ_n − 1) ≈ z_* + n π/(2 Re a);
  - smoothness defect m − 2 ≈ C(λ−1)^{3/2} e^{−Im Φ} cos(Re Φ + φ₀), with the exponent from the Hou–Luo fit
    p = 1.58.
- Interpretation: the n-th smooth profile carries n half-wavelengths of a standing wave in the stagnant layer.
  This matches the observation that it has n unstable modes.

## 5. Finite versus infinite ladders
- CCF: the sonic depth vanishes at λ* = 0.4536; the branch ends in a cusp with log-periodic approach, so
  there are three profiles.
- Boussinesq and Hou–Luo: the cusp is approached only as λ → 1, and the quasi-stagnant region in front of it
  supports WKB standing waves, so the ladder is infinite.

## 6. Methods and verification (to SI)
- Solver details and two numerical pitfalls (FFT/Gibbs contamination of Newton directions; roundoff-level
  Jacobian-vector products).
- Resolution tables:
  - Nb 32/48/64; hs 0.025/0.0125; s_start −12/−20; s_max 100/130.
  - Hou–Luo N = 8192–131072.
- Independent residual check with finite differences in s.
- Limitations:
  - numerical evidence plus asymptotics, not a proof;
  - the deepest crossings rely on errors ≤ 1e-9 in m;
  - the published λ₂, λ₃ of Wang et al. could not be accessed for a digit-by-digit comparison.

## Table 1
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
