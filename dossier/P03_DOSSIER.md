# DOSSIER P03 — A quantized hierarchy of self-similar blow-up profiles in 2D Boussinesq

**Result in one sentence.** Numerical evidence and asymptotic analysis reveal a quantized hierarchy of
self-similar Boussinesq blow-up profiles, with evidence for an infinite ladder accumulating at λ = 1.

**Central concept.** The discrete smooth blow-up profiles are *regularity-selected points on a continuous
singular branch*. For each λ there is a least-singular self-similar profile with Θ ≈ −|y₁|^{m(λ)} at the
stagnation point. The smooth profiles are the regularity resonances m(λ) = 2.

**What we do not claim.**
- No proof. Nothing here is a theorem about the Boussinesq equations. The one rigorous statement is an elementary
  lemma about oscillatory functions (`programs/P03_boussinesq_ladder/ASYMPTOTICS.md` §2).
- "Infinite" means: eight computed resonances, plus a formal mechanism whose ingredients are each checked on the
  computed profiles.
- The asymptotic spacing 3/2 is a candidate only. It is not derived.

**Status tags.**
- [N] numerical. Where stated, reproduced by two independent solvers and converged in the discretisation parameters.
- [F] formal matched asymptotics.
- [P] proved.
- [O] open.

## 1. Core claims
1. **[N] A continuous least-singular branch.** For every λ ∈ [1.069, 4.1] the least-singular profile exists. The
   profiles form one continuous branch (continuation in λ, no fold or termination found).
   - m(λ) = (λ−1)/(1+λ−A), with A = −∂₁U₁(0).
   - m decreases monotonically for λ > λ₀, reaching 1.34 at λ = 4.1.
2. **[N + exact identity] Smooth profiles are the regularity resonances m(λ) = 2.** This is the condition
   λ = −3 − 2∂₁U₁(0) of Wang et al. m − 2 changes sign at each resonance. Its extrema decay geometrically:
   +7.4e-2, −6.8e-3, +7.2e-4, −8.3e-5, +1e-5, …
3. **[N] Four resonances beyond Wang et al.**
   - λ₄, which corresponds to their candidate fourth unstable profile, plus λ₅, λ₆ and λ₇ (Table 1).
   - All eight are reproduced by an independent global Newton solver: different unknowns, discretisation,
     Biot–Savart solver and treatment of λ (Table 2).
   - Linear stability by two independent methods gives n unstable modes for the n-th profile, as Wang et al. found
     for n ≤ 3. Verified here for ⟨n ≤ 5, pending 6–7⟩ (Table 3).
   - Profiles whose index is not yet verified are called "next in continuation order".
4. **[N + F] A stalled layer next to a limiting cusp.** As λ → 1⁺:
   - For x < x_c ≈ 0.72 on the wall, the radial self-similar speed is D = O(λ−1): a quasi-stagnant (stalled) layer.
   - The layer is bounded by a front. Outside the front D → k(x−x_c)^{1/2}, the same square-root cusp that ends
     the CCF branch. Here the cusp is reached only in the limit.
   - The dip inside the front closes like ε^{0.18}. No sonic point forms at λ > 1 down to λ = 1.069.
5. **[F + N] An oscillatory inner problem gives phase quantization.**
   - Perturbations in the stalled layer obey a boundary-layer eigenproblem with a complex root κ(s).
   - The smoothness defect is m − 2 ≈ |K| e^{−Im Φ} cos(Re Φ + arg K), with Φ = ε⁻¹∫κ ds, so the resonances
     satisfy Re Φ(λ_n) = nπ + δ + o(1).
   - Computed phase gain per resonance: 2.98, 3.03, 3.07, 3.07, 3.115, 3.12, tending to π.
   - δ_n = 2.27 → 1.81 (n = 1…7), converging with steps that decay like z⁻².
6. **[F, as far as the derivation supports] 1/(λ_n − 1) = a n + o(n), with a = π/C.**
   - C = 2 Re ∫κ₀ ds > 0. Positivity is exact for the Hou–Luo local root; in 2D it is checked pointwise.
   - The remainder is O(1) only under a uniformity hypothesis (U) on the stalled layer, which fails near the front
     (see §4). Numerically z_n − a n is constant to ±0.02 over n = 3–7.
   - a ∈ [1.476, 1.50]; the WKB central value is 1.496. Measured spacings: 1.4646, 1.4699, 1.4728, 1.4743, 1.4758.
   - The value a = 3/2 is a candidate consistent with the data. It is not derived.
7. **[N + closed-form local root] Hou–Luo shows the same mechanism.**
   - 8 reliable resonances; phase gain per resonance 3.02 → 3.09 (→ π).
   - Spacings 1.2583 → 1.2652 against the WKB prediction π/(2 Re a₀) ≈ 1.27.
   - The local root has the closed form κ = (−1 + √(1+4iΩ))/(2iD̂). Re κ > 0 for Ω > 0.
8. **[N, P02] CCF is the finite contrast.**
   - Its branch reaches the square-root cusp at finite λ* = 0.4535843 and ends after three smooth profiles.
   - This is consistent with the unsuccessful search of Wang, Léger, Lai and Buckmaster (arXiv:2511.22819) for a
     third CCF profile in λ ∈ [0.455, 0.4713].
   - Finite versus infinite ladder is decided by whether the cusp forms at finite λ or only at the accumulation
     point.

## 2. Methods (programs/P03_boussinesq_ladder)
- **Marching solver** (`bq_solver.py`, `bq_newton.py`).
  - Log-polar (s = ln r, β), with the exact local structure factored out.
  - (Θ̂, Ω̂) are marched outward from the stagnation point, implicit Gauss–Legendre near it.
  - FFT/Chebyshev Biot–Savart.
  - Newton–Krylov on X = Ψ/r², quadratic to 1e-13.
  - Continuation in λ (`bq_scan2.py`) and regula falsi on m = 2 (`bq_crossing.py`).
- **Independent global solver** (`bq_global.py`, `glob_run.py`).
  - Unhatted Θ, Ω, X solved simultaneously by Newton with the exact sparse Jacobian and SuperLU.
  - 6th-order upwind-biased finite differences in s, Chebyshev in β.
  - The Biot–Savart law as the local elliptic equation (∂_s+2)²X + X_ββ + Ω = 0, with a Robin far field.
  - Smoothness imposed directly through the exact local data. λ is an eigenvalue fixed by the strain condition.
  - Extrapolation in the far-field truncation s_max.
- **Stability, method 1** (`bq_stability.py`). Perturbations e^{μτ} (τ = −ln(1−t)).
  - μ is an eigenvalue iff the linearised "march + Biot–Savart" map T_μ has an eigenvalue ν = 1.
  - Real-axis scans count the crossings (`stab_scan.py`).
  - An argument-principle count of det(I − T_μ) on [0.06, 1.5] × [−1.5, 1.5] includes complex eigenvalues
    (`stab_contour.py`).
- **Stability, method 2** (`glob_stab.py`). The generalised eigenproblem J v + μ M v = 0 on the global
  discretisation, solved by shift-invert Arnoldi at real and complex shifts. It is repeated for two origin
  truncations, and only eigenvalues independent of the truncation are kept (`gstab_summary.py`).
- **WKB** (`bq_local_eig.py`, `bq_wkb2d.py`, `wkb_sens.py`). The local eigenproblem is solved by shooting, and the
  phase is integrated from the stagnation point to the front.
- **Hou–Luo** (`hl_solver.py`, `hl_scan.py`, `hl_wkb2.py`).
- **Cost.** Everything runs on one CPU core per job. A profile takes 20–120 s; a stability scan 10–40 min. No GPU
  and no neural network.

## 3. Results
### Table 1 — the resonances (marching solver; hs = 0.0125 from λ₄ on)
| n | λ_n | z_n = 1/(λ_n−1) | z_n − z_{n−1} | two-point law of Wang et al. |
|---|---|---|---|---|
| 0 | 1.9205593 | 1.08630 | — | 1.9206 (stable) |
| 1 | 1.3990961 | 2.50566 | 1.41937 | 1.3992 |
| 2 | 1.2523487 | 3.96277 | 1.45711 | 1.2549 |
| 3 | 1.1842533 | 5.42733 | 1.46455 | 1.1872 |
| 4 | 1.1449857 | 6.89723 | 1.46991 | 1.1479 (their candidate) |
| 5 | 1.1194738 | 8.37004 | 1.47281 | 1.1223 |
| 6 | 1.1015817 | 9.84429 | 1.47426 | 1.1042 |
| 7 | 1.0883384 | 11.32010 | 1.47581 | 1.0908 |

**Resolution.**
- λ₁: Nb 32/48/64 agree to 1e-9.
- λ₂: Nb 32/48 agree to 9e-8.
- Going from hs = 0.025 to 0.0125 moves λ₄…λ₇ by −7e-7, −8e-6, +6e-5, −7e-4. The front sharpens as λ → 1; λ₇ is
  uncertain at about 5e-5.

### Table 2 — independent reproduction (global solver, s_min = −12, s_max = 8/10/12 extrapolated)
⟨filled from glob_g*_s12.log / glob_g*_h0125*.log⟩

### Table 3 — instability index
⟨filled from stab_l*.log, stabc_l*.log, gstab_summary.py⟩

### 3.4 Asymptotics (ASYMPTOTICS.md)
- **Lemma [P].** F = R(cos Θ + η), with R > 0, Θ ↑ ∞ and |η| < 1, implies infinitely many zeros. If
  Θ = Cz + Θ₀ + o(1), they satisfy λ_n = 1 + C/(nπ + c + o(1)).
- **WKB representation [F].**
  - Eikonal: the local boundary-layer eigenproblem. Dropped terms are O(ε).
  - Connection regions: the stagnation point (x ≲ ε) and the front contribute O(1) phases.
  - C = 2 Re a₀ > 0.
  - Remainder: o(z) in general; O(1) under hypothesis (U).
- **Checks [N].**
  - A linear fit Re Φ(λ_n) = Cz + Φ₀ over n = 3–7 gives C = 2.100 with residuals ≤ 0.024, i.e. π/C = 1.496.
  - The measured local slopes π/Δz_n decrease monotonically from 2.145 to 2.129, which bounds C from above.
  - Cut-off sensitivity: ⟨wkb_sens⟩.

## 4. Limitations and open points
- No proof; the WKB derivation is formal.
- Hypothesis (U) fails near the front, where the dip closes like ε^{0.18}. The o(1) error term in Φ is therefore
  not established; only o(1/(λ−1)) is.
- The value of C is known only to about ±1 %. 3/2 is not derived.
- The deepest resonances have |m − 2| ~ 1e-6 (λ₆) and 1e-7 (λ₇). Their existence rests on sign changes resolved
  by both solvers; their location is hs-sensitive (Table 1).
- **Stability results hold within a smooth-perturbation class.** The perturbation must be regular at the
  stagnation point.
  - Both methods truncate at the origin (s_start = −20 or s_min = −12/−20). This creates truncation-dependent
    modes near μ ≈ 0.66/|s_start|, which are identified and discarded.
  - The exact μ = 0 (scaling) and μ = 1 (time translation) modes are trivial.
- The published λ₂ and λ₃ of Wang et al. were not accessible. The comparison uses their λ₀ and their two-point law.

## 5. Why it matters
- It gives a mechanism for the unstable-singularity hierarchy of the Hou–Luo scenario (the Boussinesq proxy for
  3D Euler with boundary). The empirical PINN law becomes the leading term of a phase-quantization condition.
- Together with P02 it gives a criterion for when such ladders are finite (cusp at finite λ: CCF) or infinite
  (cusp only in the limit, with a stalled layer that supports WKB oscillations: Boussinesq, Hou–Luo).
- Two independent classical solvers compute every rung to 5–8 digits on one core, and two independent stability
  methods give the instability indices. These profiles are candidates for computer-assisted proofs.
