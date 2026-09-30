# A common quantization mechanism organizes self-similar fluid singularities and their stability

*Main-text draft (Article format, four figures). Every number is produced by a script in
`programs/P03_boussinesq_ladder` (P02 for CCF) and logged there. The claim levels in Box 1 bind the text. The long
technical draft `P03_MANUSCRIPT_DRAFT.md`, together with ASYMPTOTICS.md, INSTABILITY_LADDER.md and
PREDICTIONS_IPM.md, is the Supplementary Information.*

**Status of this draft:**
- All numbers are final for this version.
- Before submission, compare our IPM λ₂–λ₄ with the published neural-network values (Wang et al. 2025;
  Wang–Léger–Lai–Buckmaster 2025). They were not accessible from this environment.

## Abstract
Self-similar blow-up in incompressible fluids comes in hierarchies. The stable profile of the Hou–Luo scenario is
accompanied by unstable profiles, the n-th with n unstable directions, now computed in several models. Why such
profiles exist, whether their number is finite, and why the n-th has exactly n instabilities has been unclear.

Here we show that a single phase organizes both the profiles and their stability.
- **Profiles.** The smooth profiles are resonances of one continuous branch of singular self-similar solutions.
  Where the boundary flow stalls, a trapped wave accumulates a phase ∝ 1/ε, and every half-turn of it produces a
  profile.
- **Stability.** The same phase, shifted by πμ/ε at the stagnation point, quantizes the growth rates μ. They form a
  lattice of step ε, and one new unstable mode enters per half-turn, so the n-th profile has exactly n.
- **The spacing constant.** For the Hou–Luo model we derive it from a limit problem, a Kutta-closed airfoil
  equation: π/C = 1.279, against 1.276–1.284 from eleven computed profiles.
- **New profiles.** We report four new 2D Boussinesq profiles (index 4–7), with exact instability counts.
- **Blind test.** In the incompressible-porous-media equation we committed predictions before computing. The
  spectral predictions held, with growth rates predicted to 0.2 % and the spectral flow to 1e-3. The predicted
  accumulation of profiles as λ → 0 failed: the IPM hierarchy follows the same linear law about a finite
  accumulation point, where its boundary flow turns sonic.
- **Finite versus infinite.** How the branch ends — a stalled layer, or a sonic cusp as in the
  Córdoba–Córdoba–Fontelos model — decides whether the hierarchy is infinite.

## Box 1 — what is proved, derived, measured
| level | statement |
|---|---|
| proved (elementary) | If the smoothness defect has the form R[cos Θ + η] with Θ → ∞ and \|η\| < 1, there are infinitely many smooth profiles, one per half-turn of Θ (ASYMPTOTICS §2). |
| derived (formal matched asymptotics) | The WKB form of the defect; the eigen-condition Re Φ − πμ/ε + θ₀ = (k+½)π; unit spacing, phase-locked spectral flow and index = n; the Hou–Luo spacing π/C_∞ from the λ → 1 limit problem. |
| measured | 2D Boussinesq rungs n = 0–7 (two independent solvers); Hou–Luo n = 0–10; IPM n = 0–6; exact index n (argument-principle counts, grid and domain variants; IPM n ≤ 4); spectral flow at 9 (HL), 5 (2D) and 4 (IPM, blind) branch points. |
| not proved | Infinitely many profiles for any of these PDEs; existence of any unstable profile in the sense of a computer-assisted proof; the 2D spacing constant; the ladder offsets; the IPM endpoint (extrapolated, not reached). |

## 1. Introduction
Whether the incompressible Euler equations can form a singularity in finite time from smooth initial data is a
central open problem of mathematical fluid dynamics. Near a solid boundary it has been settled. Luo and Hou
identified a scenario of explosive growth at the wall of a rotating cylinder, and Chen and Hou proved, with
computer assistance, that it gives stable, nearly self-similar blow-up of the axisymmetric Euler equations and of
their two-dimensional proxy, the Boussinesq equations with boundary.

Stable blow-up is not the whole story. Unstable self-similar solutions are thresholds between qualitatively
different evolutions. They are the natural candidates wherever the stable scenario is unavailable, and they
organize the phase space near blow-up. Using physics-informed neural networks, Wang et al. recently found the first
few unstable profiles in three models:
- the Boussinesq equations with boundary;
- incompressible porous media (IPM);
- the one-dimensional Córdoba–Córdoba–Fontelos (CCF) equation.

In each, the profiles were ordered by their number of unstable directions, and the blow-up rates followed
empirical straight lines. This
leaves three questions. Why should such hierarchies exist? Are they finite or infinite? And why is the instability
index of the n-th profile exactly n?

Here we answer all three with a single mechanism, and test it in a model on which it was not built. The profiles
are resonances of one continuous branch of singular solutions. The resonance condition is set by the phase of a
wave trapped in a stalled boundary layer. The same phase, shifted at the stagnation point, quantizes the unstable
spectrum.

## 2. Smooth profiles are resonances of one singular branch (Figs 1, 2a)
**The branch.** Write a self-similar solution as a profile in y = x/(1−t)^{1+λ}. For each λ there is a
*least-singular* profile. Near the stagnation point on the boundary its transported scalar behaves as c|y₁|^m, with
the exponent fixed by the local strain A:
- m = (λ−1)/(1+λ−A) for the Boussinesq temperature;
- m = λ/(1+λ−A) for the IPM density.

Generically m is not an integer and the profile is singular at the stagnation point. It is smooth exactly when
m(λ) = 2.

**Following it.** We follow this branch continuously in λ with a marching Newton–Krylov method and locate its
intersections with m = 2 (Methods). In 2D Boussinesq we find eight: n ≤ 3 reproduce Wang et al., and n = 4–7 are
new. A second solver, a global sparse Newton method sharing no numerical component with the first, reproduces all
eight to 1e-5 or better. In the one-dimensional Hou–Luo model, the boundary restriction of Boussinesq, we find
eleven.

**Its shape.** Along the branch, m − 2 oscillates with an amplitude that decays by a factor of 10–20 per half period
(Fig. 2a). The smooth profiles are the zeros of this oscillation. With ε the exponent of the transported scalar
(ε = λ − 1 for Boussinesq and Hou–Luo, ε = λ for IPM), they are nearly equally spaced in 1/ε (Fig. 1). The
empirical linear laws are the leading behaviour of this resonance sequence.

## 3. A stalled layer carries the phase (Fig. 2)
**The stalled layer.** As ε decreases, the self-similar flow along the boundary stalls. Between the stagnation
point and a front at x_c, the radial speed D = V₁/x is O(ε), while it is O(1) beyond (Fig. 2b).
- In this layer the vorticity is slaved to the gradient of the transported scalar.
- Linear waves behave as exp(i∫κ ds/ε), with κ the root of a local eigenproblem. In Hou–Luo it is explicit:
  κ = (−1 + √(1 + 4iΩ))/(2iD̂), with D̂ = D/ε.

**The defect is a wave.** The smoothness defect m − 2 is carried by this wave. It has the form
R cos(Re Φ(λ) − δ), with the phase Φ = ε⁻¹∫κ ds ≈ C/ε. The smooth profiles therefore occur once per half-turn of
Re Φ:

    1/ε_n = (π/C) n + b + o(1).                                                             (1)

**Why it is centred on zero.** A non-smooth corner (m ≠ 2) induces a boundary velocity ∝ (m − 2)x^{m−1}, which is
incompatible with a stalled layer at every order in ε. So m − 2 is beyond all orders in ε, and it oscillates about
zero: every half-turn produces a new profile. An elementary lemma (Box 1) turns this form into infinitely many
profiles accumulating at ε = 0.

**The constant C is computable.** At ε = 0 the Hou–Luo profile obeys a limit problem.
- On the stalled layer, U = −2ξ becomes an airfoil equation, H[Ω] = −2.
- Beyond the front, the shed vorticity is transported by the flow it induces.
- The problem has a one-parameter family of solutions, labelled by the vorticity shed at the front. The
  finite-ε profiles shed an ever larger vorticity as ε → 0, because the dip of the stalled layer just behind
  the front deepens (Supplementary). The corresponding limit, in which the edge singularity of the layer vorticity disappears
  (a Kutta condition), gives C = 2.4565.

Hence π/C = 1.279 ± 0.001. The measured phase coefficients converge towards this value (Fig. 2c), and the eleven
Hou–Luo rungs, whose spacings increase monotonically to 1.270, extrapolate to 1.276–1.284 (Fig. 2d). The spacing
of the hierarchy is thus derived, not fitted. The 2D analogue of the limit problem, a stalled wall layer carrying
a vortex sheet, has not been solved; there the spacing is measured (1.48–1.50).

## 4. The same phase quantizes the instabilities (Fig. 3)
Perturbations of a profile grow as e^{μτ}, with τ = −ln(1−t). Two facts tie the growth rate μ to the profile
phase.
1. **In the stalled layer μ does not change the phase of the waves.**
   - The root shifts by exactly iμ/D̂: κ(μ) = κ(0) + iμ/D̂.
   - This is exact in Hou–Luo and holds to 1e-10 for the 2D local eigenproblem.
   - It is the WKB form of the identity L_μ(Θφ) = Θ L_{μ+ε}φ: growth only rescales the amplitude.
2. **At the stagnation point μ changes the local exponent,** from m to m(1 − μ/ε). The Mellin phase of the wave
   emerging there therefore shifts by πμ/ε.

Matching the two gives one eigen-condition:

    Re Φ(λ) − πμ/ε + θ₀ = (k + ½)π,   k ∈ ℤ.                                               (2)

Equations (1) and (2) share the phase Φ, and three consequences follow.

**(i) A lattice of step ε.** The unstable growth rates of every profile form a lattice of step ε, with finite-ε
corrections. The gap between the two lowest members is 1.10ε (2D, n = 7) and 1.07ε (Hou–Luo, n = 8), and it
extrapolates to 1.003ε and 1.002ε (Fig. 3a, d).

**(ii) The lattice moves with the phase.** Along the continuous branch between two profiles, (2) moves the whole
lattice up by exactly one step. A new member crosses μ = 0 once per interval, at a position fixed by the lattice
offset. We computed the spectra at 9 (Hou–Luo) and 5 (2D) branch points between profiles. The offsets follow the
phase to within 0.011 and 0.03, and each new member enters at its predicted position (Fig. 3b).

**(iii) The index is n.**
- The rungs (1) and the zero crossings of (2) are two lattices in Re Φ with the same step π, so they interlace.
- Exactly one unstable mode is added between consecutive profiles, and the n-th profile has index(U₀) + n = n.
- This is an oscillation theorem for blow-up profiles, analogous to the n nodes of the n-th Sturm–Liouville
  eigenfunction.
- Argument-principle counts of det(I − T_μ) confirm that the index is exactly n, with all modes real. They cover
  n ≤ 7 in 2D (stable under changes of grid, domain and origin truncation) and n ≤ 10 in Hou–Luo.
- Eigenvalues with large |Im μ| are excluded because the spectral radius of T_μ falls below one there.

## 5. A blind test: incompressible porous media (Figs 1, 3c)
The theory was built on Boussinesq and Hou–Luo. To test it, we applied it to a third model, IPM. IPM differs from
Boussinesq in having no vorticity transport: the vorticity is slaved to the density gradient, Ω = −∂₁R.

**Protocol.** Before computing anything beyond a single solver test, we wrote down six predictions and committed
them to a public repository (PREDICTIONS_IPM.md; git 4dc82af):
- P1: an infinite ladder with no fold or sonic point;
- P2: a linear rung law, with spacings increasing monotonically;
- P3: index n;
- P4: an unstable lattice of step λ_n;
- P5: rigid spectral flow between rungs;
- P6: modes localized at the front.

Numerical predictions for each next profile were then committed before that profile was computed, each generated
from the lower profiles by a fixed rule.

**Outcome (full scorecard in PREDICTIONS_IPM.md).**
- **Profiles.** We find seven IPM profiles: λ = 1.0285723, 0.4721297, 0.3149622, 0.2415661, 0.1987300, 0.1706181
  and 0.15092. The last two are computed on the finer grid, and the others are accurate to ≲1 % of a spacing.
  - The first two reproduce the published values. n ≤ 4 have been reported before, by neural-network methods.
  - Each rung from λ₃ on was predicted from the lower ones before it was computed. The errors were 6 %, 2.5 %, 2 % and
    0.1 % of a rung spacing, for λ₃ through λ₆.
  - The straight line through the first two profiles, the empirical law of the literature, is off by 28 % of a
    spacing already at λ₃, and more further on.
- **Spectrum: every parameter-free prediction held.**
  - The index is exactly n for n = 0–4 (argument-principle counts 1.95, 3.08, 4.08, 4.90 for n = 1–4, i.e. the
    trivial mode plus n).
  - All modes are real.
  - The growth rates lie on a lattice of step λ_n: U₃ has μ/λ₃ = 0.805, 1.879, 2.937, and U₄ has 0.796, 1.846,
    2.879, 3.868.
  - The lowest member of U₃ was predicted to within 0.2 %, and all four members of U₄ to within 0.2–1.2 %.
- **Spectral flow.** Between rungs every member moves linearly with the phase. The member that continues from the
  top of one rung to the top of the next was predicted to within 0.001–0.0014 at four branch points; the others
  to within 0.05. The unstable modes peak at the dip of the boundary flow, as predicted.
- **What failed: where the ladder accumulates.**
  - We predicted accumulation as λ → 0, with spacings in 1/λ tending to a constant. Instead the spacings contract
    by a nearly constant factor, 0.92 per rung (1.146, 1.057, 0.965, 0.892, 0.829, 0.765).
  - The rungs obey the one-phase linear law (1) in the shifted variable 1/(λ − λ_c), with λ_c = 0.037, the spacing
    (1.296) constant to 0.3 % over six intervals. The same fit returns the known accumulation point λ_c = 1.000 for Boussinesq and
    Hou–Luo.
  - The dip of the boundary flow deepens faster than λ and extrapolates to closure near λ ≈ 0.09.
  - Both point to an IPM branch that ends at *finite* λ. There the dip becomes sonic, as in the CCF reduction of
    IPM, rather than stalling as λ → 0.
- **What the failure teaches.** It sharpens rather than breaks the mechanism. The quantization and the spectral
  structure are common to the three models. The *location* of the accumulation point is not: it is wherever the
  boundary flow stalls. Our prediction took that to be the zero of the scalar exponent; in IPM it comes earlier.
- **Also failed.** Quantities we had calibrated on the other two models — the finite-λ correction to the lattice
  step and the displacement of its top member — were too large for IPM, which lies closer to the asymptotic lattice
  than either.

## 6. How a branch ends decides the hierarchy (Fig. 4)
An unbounded phase is not enough for an infinite hierarchy. The lemma behind (1) needs the smoothness defect to
oscillate about zero, F = R(cos Θ + η) with |η| < 1. Whether it does depends on how the branch ends.

**Stalled-layer ending.** In Boussinesq and Hou–Luo the boundary flow stalls as ε → 0. Every order of the formal
expansion is then smooth, so η → 0, and the profiles continue indefinitely, accumulating at ε = 0 (Fig. 4a).

**Sonic-cusp ending.** The CCF branch ends at finite λ* = 0.45358.
- There the characteristic speed vanishes at an interior point and the profile develops a square-root cusp
  (Fig. 4d).
- Near the cusp the defect oscillates log-periodically in the sonic depth δ, with frequency 2τ = 1.2988 from
  τ tanh(πτ/2) = ½.
- But it oscillates about p* − 2 = 0.0058 ≠ 0, with an amplitude ∝ δ. So η → ∞ and the profiles stop after
  three (Fig. 4b).

**IPM lies between.**
- Its boundary dip deepens faster than λ and appears to close at finite λ, as in its reduction CCF.
- But the cusp is reached much later than in CCF, after at least seven profiles.
- The ending type seems to be inherited from the one-dimensional reduction:
  - Hou–Luo (a stalled layer) for Boussinesq;
  - CCF (a sonic point) for IPM.
- Whether the IPM hierarchy terminates, or accumulates at its endpoint, is decided by the value of the defect
  there. That requires following the branch into the closing dip, as done for CCF.

## 7. Discussion
**What is established.** The profiles and instabilities of self-similar blow-up in these models are not a
catalogue. They are the two faces of one quantization. The same phase sets:
- where smooth profiles exist;
- how fast their perturbations grow;
- how many directions are unstable.

The mechanism uses only three ingredients: a transported scalar whose exponent can vanish, a stagnation point, and
a boundary along which the self-similar flow can stall. It should therefore apply beyond the three models
considered here, for example to other Euler-type reductions with boundary.

**Euler.** The Boussinesq system with boundary is equivalent, away from the axis, to axisymmetric Euler with swirl
near the wall. The hierarchy found here is therefore a hierarchy of candidate Euler singularities of increasing
codimension, all sharing the Luo–Hou geometry. We make no claim about Euler without boundary, or about
Navier–Stokes.

**Limitations.**
- The infinite hierarchies and the eigen-condition are derived by formal matched asymptotics. Every hypothesis is
  checked numerically, but none of this is a theorem.
- For n ≥ 1 no profile has a computer-assisted existence proof.
- The spacing constant is derived only for Hou–Luo.
- The ladder offsets are measured, not derived.

**Next step.** A rigorous proof for one unstable profile is natural. The profiles and spectra released with this
paper, which have Newton residuals of 1e-13 and are reproduced by independent solvers, can serve as its
approximate solutions.

## Figure legends
**Figure 1 | Hierarchies of self-similar blow-up profiles.**
- The axes are 1/ε_n against n, where ε is the exponent of the transported scalar: λ−1 for the temperature of 2D
  Boussinesq and Hou–Luo, λ for the density of IPM and for CCF. The n-th profile has exactly n unstable modes in
  every case computed.
- Open symbols: profiles reported before (Wang et al.; Chen–Hou for Hou–Luo n = 0). Filled: this work.
- Dashed lines continue the Boussinesq and Hou–Luo ladders with their asymptotic spacing, derived for Hou–Luo
  (π/C = 1.279).
- The CCF ladder terminates at a sonic cusp (star).

**Figure 2 | One branch, one phase.**
- **a.** Smoothness defect |m − 2| along the continuous branch of least-singular profiles. Smooth profiles (ticks)
  are its zeros, one per half-turn of the phase.
- **b.** The stalled layer: the radial self-similar speed along the boundary, normalized by its stagnation-point
  value, for deep Hou–Luo profiles. It is O(λ−1) up to the front at x_c ≈ 0.61, where it dips and then rises to
  O(1).
- **c.** The phase coefficient measured on finite-ε profiles converges to the value a₀ = 1.2283 obtained from the
  λ → 1 limit problem.
- **d.** The spacing of consecutive Hou–Luo profiles against 1/z (z = 1/(λ−1)), extrapolating to the derived
  constant π/(2a₀) = 1.279.

**Figure 3 | The same phase quantizes the instabilities.**
- **a.** Unstable growth rates μ/ε of the n-th profile for 2D Boussinesq, Hou–Luo and IPM: n real modes on a
  lattice of step ≈ ε.
- **b.** Spectral flow along the continuous Hou–Luo branch: at branch points between profiles the lattice has moved
  up with the phase, and a new member has entered through μ = 0.
- **c.** IPM, blind test. Spectra predicted before computation (open symbols) and computed (filled), at rungs and at
  branch points between them.
- **d.** The step of the lattice (gap between the two lowest members) in units of ε. It tends to 1, the value
  predicted by the eigen-condition (2).

**Figure 4 | Two ways a branch can end.**
- **a.** Stalled-layer endings (ε → 0). The defect oscillates about zero with exponentially small amplitude, so
  smooth profiles never stop.
- **b.** The CCF branch ends at a sonic cusp. As the sonic depth δ → 0 the defect oscillates log-periodically about
  p* − 2 = 0.0058 > 0, so no smooth profile exists beyond the third.
- **c.** The stalled layer (Hou–Luo).
- **d.** The sonic cusp (CCF): the characteristic speed vanishes at an interior point.

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
- Complete right-half-plane counts: the argument principle for det(I − T_μ) on [x_lo, 1.5] × [−1.5, 1.5] (IPM U₃,
  U₄: × [−2.5, 2.5]), with adaptive refinement, domain and grid variants, and a deeper origin truncation
  (s_start = −30).
- Exclusion of eigenvalues with large |Im μ|:
  - argument-principle counts on the strips 1.5 ≤ |Im μ| ≤ 3.5 (zero found);
  - a spectral radius of T_μ below one for larger |Im μ|, decaying roughly as 2/|Im μ| (Supplementary).

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
