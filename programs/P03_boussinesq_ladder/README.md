# P03 — The ladder of self-similar blow-up profiles of the 2D Boussinesq equations with boundary

**Question.** DeepMind (Wang et al., arXiv:2509.14185) found with PINNs a stable profile and three unstable
self-similar blow-up profiles (plus a candidate fourth) for 2D Boussinesq with boundary, the standard proxy for
3D axisymmetric Euler with boundary (Hou–Luo scenario). Their blow-up rates follow the empirical law
λ_n ≈ 1 + 1/(1.4187 n + 1.0863). Is this ladder infinite (λ_n → 1), or finite?

**Approach.** A classical, non-neural solver computes the whole one-parameter family of *least-singular*
self-similar profiles. For every λ the family has Θ ≈ −|y₁|^m at the stagnation point, with
m(λ) = (λ−1)/(1+λ−A(λ)) and A = −∂₁U₁(0) the strain there. Smooth profiles are exactly the crossings m(λ) = 2,
which is the 2D analogue of the CCF branch analysis in P02.

## Formulation (see the docstring of `bq_logpolar.py`)
- The self-similar equations Ω + ((1+λ)y+U)·∇Ω = ∂₁Θ and (1−λ)Θ + ((1+λ)y+U)·∇Θ = 0, with −ΔΨ = Ω and
  U = ∇^⊥Ψ, are posed on the quarter plane (symmetry).
- They are solved in log-polar variables s = ln r, β ∈ [0, π/2] (β = 0 is the boundary).
- The exact stagnation-point structure is factored out: Θ = cos^m β Θ̂ and Ω = cos^{m−1} β Ω̂.
- All characteristics leave the origin, so (Θ̂, Ω̂) are marched in s from the exact local solution at s = −20.
- The unknown is X = Ψ/r² on s ≥ −20 (velocity-gradient scale). Newton–Krylov is applied to X − BS(march(X)).
- The Biot–Savart solve uses an FFT in s of e^{(2−a)s}Ω (decaying at both ends) and Chebyshev collocation in β.

## Solver versions
| file | content |
|---|---|
| `bq_logpolar.py` | original solver: FFT-based velocity, explicit RK4 march |
| `bq_solver.py` | **production solver**: 8th-order finite-difference velocity from X, and an implicit 2-stage Gauss–Legendre march for s < s_sw (RK4 beyond) |
| `bq_newton.py` | Newton–Krylov (GMRES). `fd='central'` uses scaled central-difference matvecs, which gives quadratic convergence (1e-7 → 2e-11 → 1.6e-13) |
| `bq_scan.py` / `bq_scan2.py` | λ-continuation (secant predictor, adaptive step) with the original / production solver |
| `bq_crossing.py` | secant iteration on m(λ) = 2 (a smooth profile), with Newton at each λ |
| `bq_regrid.py` | transfer of states between grids (Nb, hs) for resolution studies |

Two numerical pitfalls were found and removed:
1. **Gibbs contamination of the Jacobian.** The FFT derivative of Ψ̃ = e^{(2−a)s}X assumes periodicity in s.
   Newton directions that do not decay at large s grow like e^{17} at s = 100, and the resulting ringing
   corrupted the Jacobian in the uniform-strain direction (residual floor 3.6e-6).
2. **Roundoff-level Jacobian-vector products.** Forward differences with ε = 1e-7 on unit-2-norm Krylov vectors
   perturb each entry by only ~1e-10; the matvecs are then roundoff-limited (floor ~1e-7). Central differences
   with the largest perturbed entry set to 1e-6 restore quadratic convergence.

Stiffness: near the stagnation point the angular transport rate is w/(V_r/r) ≈ 2Aβ/ε with ε = 1+λ−A → 0 as
λ → 1 on smooth profiles. The explicit RK4 march (h = 0.025) is unstable for λ ≲ 1.3 at Nb = 32. The implicit
Gauss–Legendre march is A-stable.

## Results (Nb = 32, hs = 0.025 unless stated)
**The ladder.** Smooth profiles are the zeros of m(λ) − 2 along one branch (`bq_crossing.py`, `cross_l*_nb32.log`):

| n | λ_n | z_n = 1/(λ_n − 1) | notes |
|---|---|---|---|
| 0 | 1.9205593 | 1.08630 | stable (Chen–Hou; Wang et al. 1.9205) |
| 1 | 1.3990961 | 2.50566 | Nb 32/48/64: 1.39909599/1.39909609/1.39909609; hs 0.0125: 1.39909601 |
| 2 | 1.2523487 | 3.96278 | Nb 48 agrees to 9e-8 |
| 3 | 1.1842533 | 5.42733 | |
| 4 | 1.1449864 | 6.89718 | candidate 4th unstable profile in Wang et al. (their line gives 1.1479) |
| 5 | 1.1194818 | 8.37040 | new |
| 6 | 1.1015817 | 9.84428 | new; hs 0.0125 value (hs 0.025 gave 1.1015235); Nb 48 changes m by −4.6e-10; hs 0.00625 running |
| 7 | 1.0890078 | 11.2350 | new; amplitude of m − 2 here ≈ 1e-7, hs refinement running |

- Extrema of m − 2: +7.41e-2, −6.80e-3, +7.24e-4, −8.34e-5, +1.0e-5, …
- For λ > λ₀, m decreases monotonically (1.34 at λ = 4.1), so there are no further profiles.
- The branch continues at least to λ = 1.069 (scan D), with no sign of termination.

**Mechanism** (THEORY.md): a quasi-stagnant boundary region x < x_c ≈ 0.72, where D = V₁/x = O(ε), bounded by a front
that becomes a square-root cusp as λ → 1. The WKB phase Φ = ε⁻¹∫κ ds, from the local 2D eigenproblem
(`bq_local_eig.py`, `bq_wkb2d.py`), gains 2.98, 3.03, 3.07, 3.11 between consecutive profiles, tending to π.
The asymptotic spacing is π/(2 Re a) ≈ 1.50 with a ≈ 1.05 − 0.6i.

**Hou–Luo** (`hl_solver.py`, `hl_scan_*.log`): 12 crossings up to z = 15.
- λ_n^{HL} = 1.99871, 1.44767, 1.28676, 1.21092, 1.16668, 1.13772, 1.11731, 1.10215, 1.09046, …
- Phase gained per profile 3.02–3.09, tending to π.
- Amplitude |m − 2| e^{−Im Φ} ∝ z^{−1.58}.
- The dip closes only as ε → 0 (D̂_min = 0.285 at z = 30.8).

**Figure:** `fig_bq_ladder.png`.

**Independent residual check** (`check_residual.py`: equations evaluated with 6th-order finite differences in s):
- Near the origin (r < 0.14) the relative residual is ~1e-7.
- Across the front region at λ₁ it is 1.6e-4 (Nb 32) and 1.7e-6 (Nb 48).
- The outer region (r ≳ 3) has angular structure that Nb = 32–48 under-resolves (residual 1e-3–1e-2).
- The λ_n themselves are insensitive to this (Nb 32/48 agree to ≤ 1e-7 in λ₁, λ₂ and 5e-10 in m at λ₆). The strain
  at the origin sees only low angular modes, and the boundary transport is resolved to 1e-11.
