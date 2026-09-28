# DOSSIER P03 — The infinite ladder of self-similar blow-ups of the 2D Boussinesq equations

**Result.** Eight smooth self-similar blow-up profiles of 2D Boussinesq with boundary (the Hou–Luo scenario of
3D axisymmetric Euler) were computed with a classical solver. All of them lie on **one continuous branch** of
least-singular profiles; the smooth ones are the zeros of a single scalar function m(λ) − 2.
- λ₀ is the known stable profile; λ₁, λ₂, λ₃ (and presumably λ₄) correspond to the unstable profiles found with
  PINNs by Wang et al. (2025); λ₅ and λ₆ are new.
- As λ → 1 the profiles develop a quasi-stagnant boundary region. The front bounding it becomes a CCF-type
  square-root cusp only in the limit.
- In the stagnant region the smoothness condition becomes a **WKB quantisation**: the complex phase
  Φ = ε⁻¹∫κ ds, computed from a local 2D eigenproblem, gains π between consecutive smooth profiles.
- Hence the ladder is infinite and accumulates at λ = 1, with 1/(λ_n − 1) ≈ const + n·1.50.
- This explains the empirical law found with PINNs and corrects its slope (1.4187 → 1.50 asymptotically).
- The 1D Hou–Luo model has the same structure: 8+ profiles and the same WKB law.

Status keys: VERIFIED (converged in all discretisation parameters, ≥2 routes), DERIVED (asymptotic argument +
numerical confirmation), PRELIMINARY.

## 1. Question and gap
- Wang et al. (arXiv:2509.14185) found for 2D Boussinesq a stable profile (λ₀ ≈ 1.9205, also Chen–Hou), three
  unstable profiles and a candidate fourth. The n-th profile has n unstable modes.
- They report λ_n ≈ 1 + 1/(1.4187 n + 1.0863), which is the line through λ₀ and λ₁.
- Whether the ladder is infinite, and why such a law should hold, was open. The same holds for the 1D CCF model:
  P02 in this repository shows that ladder is finite.
- Nearest work:
  - Chen–Huang–Li (arXiv:2604.01868): singular (unbounded) self-similar profiles of Hou–Luo/Boussinesq.
  - Huang–Qin–Wang–Wei (CMP 2025): existence of the Hou–Luo stable profile.
  - Chen–Hou: stable Boussinesq/Euler blow-up.
- None of these computes the unstable ladder with a classical method or addresses its accumulation.

## 2. Method (programs/P03_boussinesq_ladder)
- **Least-singular family.** For every λ, look for the profile with Θ ≈ −|y₁|^m at the stagnation point, where
  m = (λ−1)/(1+λ−A) and A = −∂₁U₁(0). Smooth profiles ⇔ m = 2 ⇔ λ = −3 − 2∂₁U₁(0), the identity used by
  Wang et al.
- **Solver.**
  - Log-polar coordinates (s = ln r, β), with the exact local structure factored out: Θ = cos^m β Θ̂.
  - All characteristics leave the origin, so (Θ̂, Ω̂) are marched in s from the local solution. The march is
    implicit Gauss–Legendre near the origin (stiff as λ → 1) and RK4 beyond.
  - Spectral Biot–Savart solve.
  - Newton–Krylov on X = Ψ/r² with central-difference matvecs, converging quadratically to 1e-13.
  - Two numerical pitfalls were removed; see README.
- **Continuation and crossings.** Continuation in λ with a secant predictor gives m(λ). Crossings are located by
  safeguarded regula falsi with a full Newton solve at each λ.
- **Cost.** One solve takes 20–60 s on one core (Ns × Nb = 8800 × 33 unknowns). No GPU and no neural network.

## 3. Results
### 3.1 The ladder (λ₀–λ₃ VERIFIED in Nb/hs; λ₄–λ₇ at hs = 0.0125, hs = 0.00625 check of λ₆ running)
| n | λ_n (this work) | z_n = 1/(λ_n − 1) | z_{n}−z_{n−1} | two-point law of Wang et al. |
|---|---|---|---|---|
| 0 | 1.9205593 | 1.08630 | — | 1.9206 (stable; Chen–Hou, Wang et al. 1.9205) |
| 1 | 1.3990961 | 2.50566 | 1.41937 | 1.3992 |
| 2 | 1.2523487 | 3.96277 | 1.45711 | 1.2549 |
| 3 | 1.1842533 | 5.42733 | 1.46455 | 1.1872 |
| 4 | 1.1449857 | 6.89723 | 1.46991 | 1.1479 (candidate in Wang et al.) |
| 5 | 1.1194738 | 8.37004 | 1.47281 | 1.1223 |
| 6 | 1.1015817 | 9.84429 | 1.47426 | 1.1042 |
| 7 | 1.0883384 | 11.32010 | 1.47581 | 1.0908 |

- **Resolution** (m at fixed λ):
  - λ₁: Nb 32/48/64 give 1.39909599/1.39909609/1.39909609; hs 0.025 → 0.0125 moves it by 2e-8; s_start −20 → −12 by 2e-6 (an O(e^{s_start}) effect; −20 is used).
  - λ₂: Nb 32/48 agree to 9e-8.
  - λ = 1.15: hs change 7e-11, Nb change 1.3e-10.
  - Deep crossings are hs-sensitive because the front sharpens: 0.025 → 0.0125 moves λ₄, λ₅, λ₆, λ₇ by
    −7e-7, −8e-6, +6e-5, −7e-4. The table uses hs = 0.0125.
- **Branch shape.**
  - m(λ) − 2 alternates in sign between consecutive crossings. Extrema: +7.41e-2, −6.80e-3, +7.24e-4, −8.34e-5, …
  - For λ > λ₀, m decreases monotonically (m = 1.34 at λ = 4.1), so there are no further crossings.

### 3.2 Structure as λ → 1 (DERIVED + numerics)
- **Quasi-stagnant region.** The radial speed along the boundary is D = V₁/x ≈ ε D̂ (ε = 1+λ−A → 0) for
  x < x_c ≈ 0.7. The front x_c is fixed as λ → 1.
- **Square-root cusp.** Outside the front, D → k(x − x_c)^{1/2}, the same square-root cusp that terminates the CCF
  branch (P02). Here it is reached only in the limit λ → 1.
- **No sonic point.** The dip inside the front closes like D̂_min ∝ ε^{0.18}, so D_min → 0 only at λ = 1.

### 3.3 Mechanism: WKB quantisation (DERIVED; THEORY.md)
- **Why m → 2.** A non-smooth profile (m ≠ 2) carries a velocity component ∝ (m−2) x^{m−1} on the boundary. The
  stagnant region only tolerates O(ε) velocity deviations, so m − 2 is beyond all orders in ε. The smoothness
  defect is carried by exponentially small WKB waves.
- **Local problem.** Perturbations ∝ exp(i∫κ ds/ε) with a vertical structure on the scale εx obey a local
  eigenproblem (transport + Biot–Savart; `bq_local_eig.py`). Its complex root κ is continued from the Hou–Luo
  closed form κ = (−1 + √(1+4iΩ))/(2iD̂).
- **Quantisation.**
  - m(λ) − 2 ≈ Re[K e^{iΦ(λ)}], so Re Φ(λ_n) = Φ₀ + nπ.
  - Computed phase gains between consecutive smooth profiles: 2.98, 3.03, 3.07, 3.07, 3.115, 3.12 (→ π).
  - εΦ → a ≈ 1.05 − 0.6i, so the asymptotic spacing is π/(2 Re a) ≈ 1.50. The observed spacings increase
    monotonically, 1.457 → 1.4758.
- **Amplitude.** The oscillation amplitude decays like e^{−Im Φ} times an algebraic factor (≈ z^{−1.4}).

### 3.4 Cross-check: Hou–Luo boundary model (hl_solver.py; VERIFIED)
- **Same solver design, 1D.** 8th-order quadrature, Gauss–Legendre march and exact Mellin symbol. Resolution
  N = 8192 → 65536 changes m by 1.4e-7 → 2e-10 at z = 7.45.
- **Ladder.**
  - λ_n^{HL} = 1.99871, 1.44767, 1.28676, 1.21092, 1.16668, 1.13772, 1.11731, 1.10215, 1.09047, …
  - Spacing in z → 1.267, WKB prediction 1.27.
  - Phase gain per interval 3.02–3.09 (→ π).
- **Dip.** D̂_min ∝ ε^{0.35}, no termination down to λ = 1.04 (z = 24).

## 4. Limitations
- Numerical evidence plus asymptotics, not a proof. "Infinite" rests on:
  - 7 computed crossings and the continuation of the branch to λ ≈ 1.08 (2D) / 1.04 (HL);
  - the WKB mechanism, whose ingredients (stagnant region, front, local eigenproblem) are verified on the
    computed profiles.
- Deep crossings have tiny amplitudes (|m − 2| ~ 1e-6 near λ₆, ~1e-7 near λ₇). Their existence (the sign change) is
  robust, but their location is hs-sensitive (see 3.1). λ₇ is uncertain at about 5e-5 even at hs = 0.0125.
- **Exact values of Wang et al.** were not accessible (arXiv blocked in this environment). The comparison uses
  their stated λ₀ and their two-point law. Our λ₂ and λ₃ lie 2.6e-3 and 2.9e-3 below the law, which is consistent
  with the law being a line through λ₀ and λ₁.
- **Instability counts** (n unstable modes for the n-th profile) are taken from Wang et al. for n ≤ 3 and not
  recomputed here.
- **The amplitude prefactor** (z^{−p}) is fitted, not derived.

## 5. Why this matters
- It gives the first mechanism for the infinite ladder of unstable blow-ups in an incompressible fluid model
  (the Hou–Luo scenario for 3D Euler with boundary), and a first-principles version of the empirical PINN law.
- Together with P02 it gives a unified picture of when ladders are finite or infinite:
  - **CCF:** the sonic cusp forms at finite λ*, so the branch terminates.
  - **Boussinesq / Hou–Luo:** the cusp forms only as λ → 1, and the quasi-stagnant region in front of it supports
    WKB standing waves, so the ladder is infinite.
- The classical branch method finds every smooth profile on the branch, whereas a PINN search can only target them
  one at a time. It also supplies high-precision profiles (7 digits in λ in seconds to minutes on one core) for
  computer-assisted proofs.
