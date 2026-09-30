# A common quantization mechanism organizes self-similar fluid singularities and their stability

*Main-text draft (Article format, four figures). Every number is produced by a script in
`programs/P03_boussinesq_ladder` (P02 for CCF) and logged there. The claim levels in Box 1 bind the text. The long
technical draft `P03_MANUSCRIPT_DRAFT.md`, together with ASYMPTOTICS.md, INSTABILITY_LADDER.md and
PREDICTIONS_IPM.md, is the Supplementary Information.*

**Status of this draft:** the IPM paragraphs are filled in as the blind test (PREDICTIONS_IPM.md) proceeds.
Placeholders are marked ⟨IPM⟩.

## Abstract
Self-similar blow-up in incompressible fluids comes in hierarchies. The stable profile of the Hou–Luo scenario is
accompanied by unstable profiles, the n-th with n unstable directions, and these have recently been computed in
several models. It has been unclear why such profiles exist, whether their number is finite, and why the n-th has
exactly n instabilities.

We show that a single phase organizes both the profiles and their stability.
- **Profiles.** The smooth profiles are resonances of one continuous branch of singular self-similar solutions.
  Along the branch, a wave trapped in a stalled boundary layer accumulates a phase C/ε, where ε is the exponent
  of the transported scalar, and every half-turn of this phase produces a smooth profile.
- **Stability.** The same phase, shifted by πμ/ε at the stagnation point, quantizes the growth rates μ of the
  instabilities. They form a lattice of step ε, and one new unstable mode enters per half-turn, so the n-th profile
  has exactly n.
- **Spacing constant.** For the Hou–Luo model we derive it from a limit problem: π/C = 1.279, against 1.276–1.284
  extrapolated from eleven computed profiles.
- **New profiles.** Four new 2D Boussinesq profiles (n = 4–7), eleven Hou–Luo profiles and ⟨IPM⟩ porous-medium
  profiles, each with exactly n unstable modes.
- **Blind test.** In the porous-medium equation the predictions were written down before any computation.
  ⟨IPM outcome⟩.
- **Finite versus infinite.** When the branch instead ends at a sonic cusp at finite parameter, as in the
  Córdoba–Córdoba–Fontelos model, the hierarchy is finite.

## Box 1 — what is proved, derived, measured
| level | statement |
|---|---|
| proved (elementary) | If the smoothness defect has the form R[cos Θ + η] with Θ → ∞ and \|η\| < 1, there are infinitely many smooth profiles, one per half-turn of Θ (ASYMPTOTICS §2). |
| derived (formal matched asymptotics) | The WKB form of the defect; the eigen-condition Re Φ − πμ/ε + θ₀ = (k+½)π; unit spacing, phase-locked spectral flow and index = n; the Hou–Luo spacing π/C_∞ from the λ → 1 limit problem. |
| measured | 2D Boussinesq rungs n = 0–7 (two independent solvers); Hou–Luo n = 0–10; exact index n (argument-principle counts, grid and domain variants); spectral flow at 9 (HL) and 5 (2D) branch points; ⟨IPM⟩. |
| not proved | Infinitely many profiles for any of these PDEs; existence of any unstable profile in the sense of a computer-assisted proof; the 2D spacing constant; the ladder offsets. |

## 1. Introduction
[~300 words.]
- Blow-up for 3D Euler with boundary: the Luo–Hou scenario, and the Chen–Hou computer-assisted proof of stable,
  nearly self-similar blow-up for 2D Boussinesq and axisymmetric Euler with boundary.
- Unstable self-similar solutions:
  - they are relevant where the stable one is inaccessible;
  - they are thresholds between blow-up scenarios;
  - they are candidates for viscous or boundary-free settings.
- Wang et al. (2025) computed the first few by physics-informed neural networks, in Boussinesq, IPM and CCF. They
  found empirical rules: 1/(λ_n − 1) and 1/λ_n linear in n, and n unstable modes for the n-th profile.
- Open questions: are the hierarchies infinite, where do the linear laws come from, and why is the index n?
- This paper answers all three with one mechanism, and tests it blindly in a third model.

## 2. Smooth profiles are phase resonances of one singular branch (Fig. 1, Fig. 2)
- **The branch.** For every λ there is a least-singular self-similar profile, with local behaviour
  Θ ≈ −|y₁|^m at the stagnation point. m(λ) = ε_λ/(1+λ−A), where ε_λ = λ−1 (Boussinesq, HL) or λ (IPM), and A is
  the strain. The profile is smooth iff m = 2.
- **Computation.** We continue this branch with a marching Newton–Krylov solver and find its intersections with
  m = 2.
  - 2D Boussinesq: eight, λ₀…λ₇. n ≤ 3 reproduces Wang et al.; n = 4–7 are new.
  - Hou–Luo: eleven.
  - IPM: ⟨IPM⟩.
- **Independent check.** A second solver, a global sparse Newton method sharing no numerical ingredient,
  reproduces every Boussinesq rung to 1e-5.
- **What the branch looks like.** Between rungs, m − 2 oscillates with an amplitude that decays by a factor
  ≈ 10–20 per half period (Fig. 2a).

## 3. The phase: a wave in a stalled boundary layer (Fig. 2)
- **The stalled layer.** As the scalar exponent ε → 0, the self-similar flow along the boundary stalls between
  the stagnation point and a front at x_c: the radial speed is D = O(ε).
- **The wave.** Linear waves in the stalled layer behave as exp(i∫κ ds/ε). κ is the root of a local eigenproblem,
  in closed form for Hou–Luo: κ = (−1 + √(1+4iΩ))/(2iD̂).
- **The smoothness defect** is carried by this wave, F = m − 2 ∝ Re[K e^{iΦ}] with Φ = ε⁻¹∫κ ds ≈ C/ε.
  Smooth profiles occur once per half-turn of Re Φ:

      1/ε_n = (π/C) n + b + o(1).

- **Why the oscillation is centred at the smooth value.** m ≠ 2 would induce a boundary velocity ∝ (m−2)x^{m−1},
  which the stalled layer forbids at every order in ε. The defect is therefore beyond all orders and purely
  oscillatory.

**The spacing constant is computable.**
- At ε = 0 the Hou–Luo profile solves a limit problem:
  - an airfoil equation H[Ω] = −2 on the stalled layer, with the vorticity slaved to the buoyancy gradient;
  - transport of the shed vorticity beyond the front.
- The finite-ε profiles shed a front vorticity that grows without bound. The corresponding limit gives
  C_∞ = 2.4565, so π/C_∞ = 1.279 ± 0.001, with no fit to the rung data.
- The measured spacings, 1.2325 … 1.2699 for n = 1–10, increase monotonically and extrapolate to 1.276–1.284
  (Fig. 2d).
- The 2D analogue of the limit problem remains open. The 2D spacing is measured as 1.483–1.495.

## 4. The same phase quantizes the instabilities (Fig. 3)
Two facts connect the spectrum to the profile phase.
1. **In the stalled layer the growth rate does not change the phase.** κ(μ) = κ(0) + iμ/D̂.
   - Exact in Hou–Luo.
   - To 1e-10 in the 2D local eigenproblem.
   - It is the WKB form of the identity L_μ(Θφ) = Θ L_{μ+ε}φ.
2. **At the stagnation point a growth rate μ changes the local exponent** from m to m(1 − μ/ε). The Mellin phase
   of the wave therefore shifts by πμ/ε.

Matching gives one eigen-condition:

    Re Φ(λ) − πμ/ε + θ₀ = (k + ½)π.

Its consequences, all confirmed:
- (i) The unstable growth rates form a lattice of step ε. The lower spacing extrapolates to 0.996ε in both
  Boussinesq and Hou–Luo.
- (ii) Along the branch the lattice moves rigidly with the phase, one step per rung interval. Offsets agree with
  prediction to 0.011 (HL, 9 branch points) and 0.03 (2D, 5).
- (iii) One eigenvalue crosses μ = 0 per half-turn, so index(U_n) = index(U₀) + n = n.
  - This is an oscillation theorem for blow-up profiles: the n-th profile has n unstable directions, as the n-th
    Sturm–Liouville eigenfunction has n nodes.
  - Argument-principle counts confirm that the index is exactly n, with no complex modes: n ≤ 7 in 2D and n ≤ 10
    in HL.

## 5. A blind test: incompressible porous media (Fig. 1, Fig. 3)
⟨IPM⟩
- **The predictions** (PREDICTIONS_IPM.md, committed before any IPM continuation):
  - P1: infinite ladder;
  - P2: rung law with monotone spacings;
  - P3: index n;
  - P4: spectral lattice of step λ_n;
  - P5: rigid spectral flow;
  - P6: front-localized modes.
- **Stage 2:** numerical predictions for U₃ and U₄ from the lower rungs.
- Outcome: ⟨IPM⟩.

## 6. Why some hierarchies are finite (Fig. 4)
- **The lemma's criterion.** Infinitely many zeros require the defect to oscillate about zero: η → 0 in
  F = R[cos Θ + η].
- **A stalled-layer ending (ε → 0) guarantees this,** because the formal expansion is smooth at every order.
- **The CCF branch ends differently.** It ends at a finite λ* = 0.4536, at a profile with an interior sonic
  square-root cusp. It approaches the cusp log-periodically, but about p* = 2.00577 ≠ 2. Its phase is unbounded,
  yet the hierarchy stops after three profiles.
- **So an unbounded phase is not enough.** What decides finiteness is the nature of the ending:
  - a stalled layer: infinite (Boussinesq, Hou–Luo, ⟨IPM⟩);
  - a sonic cusp: finite (CCF).

## 7. Discussion
- **Euler.** 2D Boussinesq with boundary is equivalent, away from the axis, to axisymmetric Euler with swirl near
  the wall. The hierarchy found here is therefore a hierarchy of candidate Euler singularities of increasing
  codimension, all sharing the Hou–Luo geometry.
  - We make no claim about Euler without boundary, or about Navier–Stokes.
  - Nor do we claim existence in the rigorous sense for n ≥ 1.
- **Outlook.** A computer-assisted proof for one unstable rung is the natural next step. Our profiles and spectra
  (released with the paper) are accurate enough to serve as its approximate solutions.
- **Generality.** The mechanism needs only three ingredients:
  - a transported scalar whose exponent can vanish;
  - a stagnation point;
  - a boundary on which the flow can stall.

## Figures
1. **The hierarchies.**
   - 1/ε_n against n for 2D Boussinesq, Hou–Luo, IPM and CCF. Literature profiles are open symbols, new ones
     filled; CCF terminates at its cusp.
   - Inset: instability index = n.
2. **One branch, one phase.**
   - (a) m(λ) − 2 along the branch, with the rungs as zeros.
   - (b) the stalled layer and the front.
   - (c) the phase Re Φ against 1/ε.
   - (d) the spacing constant from the limit problem against the measured spacings.
3. **The same phase quantizes the instabilities.**
   - Spectral flow along the branch (HL, 2D, IPM), with the eigen-condition lines.
   - Index staircase.
   - Spacing → ε.
4. **Infinite versus finite.**
   - Stalled-layer endings (defect oscillating about 0) against the CCF sonic cusp (oscillating about p* − 2 > 0).
   - Schematics of the two end structures.

## Methods (draft)
**Self-similar equations.** Perturbations of all three models are written in self-similar variables
y = x/(1−t)^{1+λ} and τ = −ln(1−t).
- **2D Boussinesq.** Ω + V·∇Ω = ∂₁Θ, V·∇Θ = (λ−1)Θ, with V = (1+λ)y + U, U = ∇^⊥Ψ and −ΔΨ = Ω, on the half-plane with
  no penetration.
- **Hou–Luo.** The boundary restriction of Boussinesq, with U' = HΩ.
- **IPM.** V·∇R = λR with Ω = −∂₁R.

In each case the least-singular profile at given λ is Θ (or R) ≈ c|y₁|^m at the stagnation point, with
m = ε/(1+λ−A), where A is the strain and ε the exponent of the transported scalar.

**Solver A (2D, IPM).**
- Log-polar coordinates (s = ln r, β) with the exact stagnation-point structure factored out: Θ = cos^m β Θ̂.
- Θ̂ (and Ω̂) are marched outward in s from the exact local solution at s = −20. The march is implicit
  Gauss–Legendre near the origin and RK4 beyond.
- Biot–Savart: 8th-order finite differences in s and Chebyshev collocation in β (Nb = 32, h_s = 0.025).
- Newton–Krylov on X − BS(march(X)) with X = Ψ/r², using central-difference matvecs (quadratic convergence).
- Continuation in λ, and secant iteration on m(λ) = 2 for the smooth profiles.

**Solver B (2D check).** A global sparse Newton method on a different grid and representation. It shares no
numerical ingredient with solver A and reproduces every Boussinesq rung to ≤ 1e-5.

**Hou–Luo.** A log-grid march with the Mellin multiplier of Ω ↦ U/ξ, N = 65,536–131,072.

**Linear stability.**
- The eigen-condition is that the linearised "march + Biot–Savart" map T_μ has eigenvalue ν = 1. It is validated
  on the exact time-translation mode μ = 1 (|T₁x − x|/|x| ≤ 5e-8).
- Real unstable eigenvalues: parity changes of the number of real ν > 1 on a grid in μ/ε, refined by bisection.
- Complete right-half-plane counts: the argument principle for det(I − T_μ) on [x_lo, 1.5] × [−1.5, 1.5], with
  adaptive refinement, domain and grid variants, and a deeper origin truncation (s_start = −30).
- Exclusion of eigenvalues with large |Im μ|: the spectral radius of T_μ along Re μ = x_lo for |Im μ| up to ⟨Y⟩
  (Supplementary).

**The λ → 1 limit problem (Hou–Luo).**
- Stalled layer: odd Chebyshev-weighted expansion Ω = Σb_jT_j/√(1−ξ²), 160–240 modes.
- Outer region: graded grid y = 1 + σ^p plus a logarithmic grid to 10⁸, with principal-value quadrature with local
  subtraction.
- Solution: a fixed-point iteration, then Newton–Krylov continuation in Ω_f.
- Phase integral: evaluated with θ = u² quadrature at the edge.

**IPM blind protocol.**
- PREDICTIONS_IPM.md (P1–P6) was committed (git 4dc82af) before any IPM continuation.
- The stage-2 numerical predictions were committed before the rungs they concern were computed, each from the
  lower rungs by a fixed procedure (`ipm_stage2_predict.py`).

**Data and code.** All scripts, logs, profiles and spectra are in the repository. `README.md` in
programs/P03_boussinesq_ladder lists the command, and the approximate time, for every number in the paper.
