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

## Results so far (Nb = 32, hs = 0.025)
- λ₀ (stable): m = 2 at λ = 1.92056 (DeepMind: 1.9205).
- λ₁: m = 2 at λ ≈ 1.39903 (DeepMind fit: 1.3992).
- m(λ) has a maximum ≈ 2.0740 near λ ≈ 1.68 between λ₀ and λ₁. For λ > λ₀, m decreases monotonically
  (1.533 at λ = 3); a scan to λ = 8 is running.
- Below λ ≈ 1.9 a dip develops in V_r/r along the boundary at r ≈ 0.03–0.3 (min/ε = 0.83 at λ = 1.40). This is the
  2D analogue of the sonic-point mechanism that terminates the CCF branch (P02).

(Scans continuing; see `scan2_*.log`.)
