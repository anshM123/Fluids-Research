# A single phase quantizes self-similar singularities and their instabilities in incompressible flow

*Article draft for Nature Physics (primary) or Nature Computational Science (see the two cover letters in
`P03_COVER_LETTERS.md`). Main text ≈ 3,600 words; abstract 150 words; five figures and one table; Methods.
Every number is produced by a script in `programs/P03_boussinesq_ladder` (P02 for CCF) and logged there.
Items marked **[PENDING]** depend on computations that are running (λ₇ holdout, IPM endpoint); they will be filled
in, or the claim removed, before submission.*

---

## Abstract
Incompressible flows can form self-similar singularities in hierarchies of unstable profiles, recently found by
physics-informed machine learning. Why such hierarchies exist, whether they end, and why the n-th profile has
exactly n unstable directions have remained open. Here we show that one phase, accumulated by a wave trapped in a
stalled boundary layer, answers all three questions. Smooth profiles are resonances of a single branch of singular
solutions, one per half-turn of the phase. The same phase, shifted at the stagnation point, quantizes the growth
rates and adds one unstable direction per profile. Computed from the profiles alone, the phase locates them to
0.2% of a spacing in a model it was not built on, although the defect it predicts shrinks tenfold per profile. A
pre-registered test in porous-media flow confirmed the spectral predictions and refuted the predicted endpoint,
showing that how the stalled layer ends decides where hierarchies accumulate or stop.

---

## Main text

Whether the incompressible Euler equations can develop a singularity in finite time from smooth data is among the
central open questions of mathematical physics. Near a solid boundary it has been answered. Luo and Hou identified
explosive growth at the wall of a rotating cylinder¹, and Chen and Hou proved, with computer assistance, that this
scenario gives stable, nearly self-similar blow-up of both the axisymmetric Euler equations and their
two-dimensional proxy, the Boussinesq system with a boundary². Stable blow-up, however, is only the bottom of a
hierarchy. Using physics-informed neural networks, Wang et al.³ found the first few *unstable* self-similar
profiles of the Boussinesq equations, of incompressible porous-media flow (IPM) and of the one-dimensional
Córdoba–Córdoba–Fontelos (CCF) model. These profiles are thresholds between qualitatively different evolutions and
the natural candidates wherever the stable scenario is unavailable. In every model they came ordered by their
number of unstable directions, with blow-up rates on empirical straight lines. Three questions follow. Why do such
hierarchies exist? Are they finite? Why does the n-th profile have exactly n unstable directions?

We answer all three with one mechanism, test it on a model it was not built on, and turn it into a computational
observable that remains accurate where the profiles themselves become exponentially hard to locate.

### Hierarchies are resonances of one singular branch
We write a self-similar solution as a profile in y = x/(1 − t)^{1+λ}, with blow-up rate λ. For each λ there is a
*least-singular* profile. Near the stagnation point on the wall its transported scalar (temperature Θ in
Boussinesq, density R in IPM) behaves as c|y₁|^m. The exponent is set by the local strain A:

    m(λ) = ε / (1 + λ − A),       ε = λ − 1 (Boussinesq, Hou–Luo),   ε = λ (IPM),

where ε is the exponent of the transported scalar. Generically m is not an integer and the profile is singular at
the stagnation point; it is smooth exactly when m = 2. We follow this branch continuously in λ by Newton–Krylov
continuation and locate its intersections with m = 2 (Methods).

The smooth profiles found by neural networks are such intersections, and there are more of them (Fig. 1):
- **2D Boussinesq: eight**, four of them new (n = 4–7). A second solver that shares no numerical component with the
  first reproduces all eight to ≤ 1e-5.
- **Hou–Luo** (the boundary restriction of Boussinesq): eleven.
- **IPM: seven**, two of them new (n = 5, 6). **[PENDING: λ₇.]**

Along the branch, m − 2 oscillates with an amplitude that falls by a factor of about ten per half-period, and the
smooth profiles are its zeros (Fig. 2a). A defect of this form, R(λ)[cos Θ(λ) + η] with Θ → ∞ and |η| < 1, has
infinitely many zeros (Methods). The questions are what Θ is, and what decides η.

### A trapped wave carries the phase
As ε decreases, the self-similar flow along the wall stalls. Between the stagnation point and a front, the wall
speed D = V₁/x is O(ε), compared with O(1) beyond (Fig. 2b). In this layer, perturbations travel as short waves,
exp(i∫κ ds/D₀), where s = ln x and D₀ = ε/m is the stagnation-point speed. The local wavenumber κ is the root of an
eigenproblem across the layer. It is closed-form for Hou–Luo, and an ordinary-differential eigenproblem for 2D
Boussinesq and IPM (Methods).

The smoothness defect is carried by this wave, and Θ is its phase,

    Φ₀(λ) = (1/D₀) ∫ κ(s; λ) ds          (stagnation point → front),

so the smooth profiles satisfy

    Re Φ₀(λ_n) = nπ + δ.                                                          (1)

Two properties make Φ₀ more than a fitting device.

**It is computable from the profile alone.** κ needs only O(1)-accurate wall data of the profile, never the
exponentially small defect.
- In IPM, a model on which the theory was not built, Re Φ₀(λ_n) − nπ = 2.033, 2.023, 2.032, 2.037, 2.028 for
  n = 1–5. That is constant to ±0.007, i.e. 0.2 % of a spacing, with no adjustable parameter beyond the front
  cut-off.
- At n = 6 the offset is 1.984 (Fig. 2c).
- In Hou–Luo the phase gains 3.02–3.09 per profile, tending to π.

**In IPM it takes an exact, transparent form.** Because the vorticity is slaved to the density gradient
(Ω = −∂₁R), rescaling the wall-normal coordinate by D̂ = D/D₀ removes the wall speed from the local problem. So
κ = K(G, μ, R_y)/D̂ exactly, with G and R_y the wall gradients of R and μ the wall-normal strain, and

    Φ₀ = ∫ K ds / D.                                                              (2)

The phase is the integral of a local wavenumber over the unnormalized wall speed: the slower the wall flow, the
more phase a wave accumulates.

**For Hou–Luo the spacing constant follows from a limit problem.**
- At ε = 0 the stalled layer obeys an airfoil equation, H[Ω] = −2.
- The vorticity shed at the front is carried by the flow it induces.
- The edge singularity of the layer vorticity must vanish, a Kutta condition. This fixes Φ₀ ≈ C/ε with C = 2.4565.
- Hence the asymptotic spacing π/C = 1.279 ± 0.001, against 1.276–1.284 extrapolated from eleven computed profiles
  (Fig. 2d).

The spacing of the hierarchy is thus derived, not fitted.

### The same phase quantizes the instabilities
Perturbations of a profile grow as e^{μτ}, with τ = −ln(1 − t). The growth rate enters the phase at only one place:
- **In the layer it does not change the phase.** The wave root shifts by exactly iμ/D̂, κ(μ) = κ(0) + iμ/D̂. This
  only rescales the amplitude, as the identity L_μ(Θφ) = Θ L_{μ+ε}φ of the transport operator requires.
- **At the stagnation point it does.** It changes the local exponent from m to m(1 − μ/ε), which shifts the Mellin
  phase of the emerging wave by πμ/ε.

Matching gives one eigen-condition with the same phase:

    Re Φ₀(λ) − πμ/ε + θ₀ = (k + ½)π,     k ∈ ℤ.                                    (3)

Equations (1) and (3) share Φ₀, and three consequences follow (Fig. 3).

**A lattice.** The unstable growth rates of each profile form a lattice of step ε. The gap between the two lowest
members is 1.10ε (2D, n = 7) and 1.07ε (Hou–Luo, n = 8), extrapolating to 1.003ε and 1.002ε.

**Rigid flow.** Between consecutive profiles, (3) moves the whole lattice up by one step, and one new member enters
through μ = 0. At 9 (Hou–Luo) and 5 (2D) branch points between profiles, the lattice offsets follow the phase to
within 0.011 and 0.03.

**An oscillation theorem.** The profiles (1) and the zero crossings of (3) are two lattices in Re Φ₀ with the same
step π, so they interlace. Exactly one unstable direction is added per profile, and the n-th profile has index n,
as the n-th Sturm–Liouville eigenfunction has n nodes. Argument-principle counts of det(I − T_μ) confirm index n
with all modes real:
- n ≤ 7 in 2D, stable under changes of grid, domain and origin truncation;
- n ≤ 10 in Hou–Luo;
- n ≤ 4 in IPM.

### Computing beyond the exponential wall
The defect that defines a smooth profile is exponentially small. In IPM its extrema fall from 8 × 10⁻² to
5 × 10⁻⁸ over seven half-periods (Fig. 4a). Locating a deep profile therefore requires solving for m to a precision
that grows tenfold per profile.

We audited every source of error at the deepest computed profiles (Fig. 4b, Methods):

| source | effect on m |
|---|---|
| Newton floor | about 5 × 10⁻¹⁰ |
| origin truncation | an offset that changes sign: 2 × 10⁻⁷ at s_start = −14, 5 × 10⁻⁹ at −16, ≤ 3 × 10⁻¹⁰ at −20 |
| angular resolution | ≤ 1 × 10⁻⁹ at s_start = −20 (but 2 × 10⁻⁸ when combined with a shallow truncation) |
| radial resolution | ≤ 1 × 10⁻⁹ |
| far field | 2 × 10⁻¹⁰ |

The converged production settings place the sixth IPM profile at λ₆ = 0.15092 ± 0.00005 (0.25 % of a spacing).
Two profiles deeper, the defect itself reaches this noise floor.

The phase does not suffer from this wall. Equation (2) needs the profile only to O(1) accuracy, so Φ₀ can be
evaluated on profiles that are far easier to compute than the defect they encode. That makes it a predictor for
profiles that direct solvers cannot locate. Spurious resonances at under-resolved depth are exactly what this
removes.
- **[PENDING]** Holdout test at the seventh IPM profile (below).

This is the computational counterpart of the physical mechanism: a resonance observable that turns an
exponentially ill-conditioned search into a well-conditioned quadrature.

### A blind test in porous-media flow
The theory was built on Boussinesq and Hou–Luo. To test it, we applied it to IPM, where the vorticity has no
transport equation of its own. Before computing anything beyond a single solver check, we committed six
predictions to a public repository: an infinite ladder accumulating as λ → 0, a linear law with increasing spacings,
index n, a lattice of step λ_n, rigid spectral flow, and localization at the front. Each further profile and
spectrum was then predicted from the lower ones by a fixed rule and time-stamped before it was computed (Table 1).

**The spectral predictions held.**
- The index is n for n = 0–4.
- The growth rates of the fifth profile were predicted to 0.2–1.2 %.
- The member continuing between profiles was predicted to 10⁻³ at four branch points.

**The predicted accumulation failed.** The spacings of 1/λ_n contract by a nearly constant factor, 0.92 per profile,
instead of tending to a constant.

**Why it failed.** The anatomy of the phase explains the failure (Fig. 5b):
- The contribution of the inner layer is independent of λ.
- The dip ahead of the front collapses as it deepens: its width scales like the cube of its depth, so its
  contribution stays bounded.
- All the growth comes from the front receding outwards, from x = 0.68 to x = 0.99 as 1/λ goes from 2 to 7.
- The stalled layer therefore lengthens as λ decreases, and Φ₀ grows faster than 1/λ.

**What it fixes.** A shifted law, 1/(λ_n − λ_c) linear in n with λ_c = 0.037 ± 0.003, fits the seven profiles, but
it is a description, not a mechanism. The receding front fixes the mechanism.
- **[PENDING]** Pre-registered holdout predictions for λ₇ (rung fits z₇ = 7.334–7.342; phase rule z₇ = 7.352) and
  its matched-resolution computation.
- **[PENDING]** The theory-based λ₈ prediction from the derived endpoint class.

### How a hierarchy ends
A divergent phase is necessary for an infinite hierarchy but not sufficient. The defect must also oscillate about
zero (η → 0). Together the two conditions classify how hierarchies end (Fig. 5).

**Stalled layer of fixed extent (Boussinesq, Hou–Luo).**
- The front converges to x_c ≈ 0.61–0.72 and Φ₀ ≈ C/ε.
- Every order of the formal ε-expansion is smooth, so the defect is beyond all orders and centred.
- The hierarchy is infinite, accumulates at ε = 0, and its spacings increase towards π/C.

**Sonic cusp (CCF).** The branch ends at finite λ* = 0.45358, where the characteristic speed vanishes inside the
domain.
- The phase still diverges, logarithmically in the sonic depth, with frequency 2τ = 1.2988 from τ tanh(πτ/2) = ½.
- But the defect oscillates about p* − 2 = 0.0058 ≠ 0, so η → ∞ and the hierarchy stops after three profiles.

**Receding stalled layer (IPM).**
- The defect is centred, as in Boussinesq, but the layer lengthens, so the spacings contract.
- **[PENDING]** Whether the front recedes indefinitely, so that profiles accumulate at λ = 0 with contracting
  spacings, or the dip closes into a sonic point at finite λ, so that the hierarchy terminates with a finite phase.
  The resolved deep branch decides this; the scaling (2) gives the divergence class in each case.

### Discussion
The profiles and instabilities of self-similar blow-up in these models are not a catalogue. They are two faces of
one quantization, and the same phase sets three things:
- where smooth profiles exist;
- how fast their perturbations grow;
- how many directions are unstable.

The mechanism needs only three ingredients: a transported scalar whose exponent can be small, a stagnation point,
and a wall along which the self-similar flow stalls. What is model-specific is how the stalled layer ends:
- a fixed front gives an infinite ladder with a derived spacing;
- a sonic cusp gives a finite one;
- a receding front gives contracting spacings.

**Euler.** Near the wall and away from the axis, axisymmetric Euler with swirl reduces at leading order to
Boussinesq with boundary. Chen and Hou's proof treats Euler as a perturbation of Boussinesq in this way². Each
profile found here is therefore a candidate Euler singularity at a wall, of codimension n. Carrying an unstable
profile over to Euler requires controlling that perturbation along its n unstable directions, which we have not
done.

**Computation.** The resonance phase is an observable that a computer can evaluate where the object it predicts is
exponentially small. It should be useful wherever unstable self-similar solutions are sought by neural or
classical solvers, including Navier–Stokes and other boundary-driven singularities. Our profiles and spectra,
released with this paper, have Newton residuals of 10⁻¹³ and serve as approximate solutions for computer-assisted
proofs.

**Limitations.**
- The infinite hierarchies and the eigen-condition are derived by formal matched asymptotics. Every hypothesis is
  checked numerically, but none is a theorem.
- No profile with n ≥ 1 has a computer-assisted proof.
- The spacing constant is derived only for Hou–Luo.

---

## Table 1 | The porous-media blind test (all predictions time-stamped in the repository before computation)
| prediction (commit) | outcome |
|---|---|
| P1 accumulation as λ → 0; no fold or sonic point (4dc82af) | **failed** in its accumulation clause: spacings contract by 0.92 per profile |
| P2 linear law, increasing spacings (4dc82af) | **failed**: spacings decrease |
| P3 index n, all modes real (4dc82af) | **held** for n = 0–4 (contour counts 0.99, 1.95, 3.08, 4.08, 4.90) |
| P4 lattice of step λ_n (4dc82af) | **held**: gaps 1.12 (U₂); 1.07, 1.06 (U₃); 1.05, 1.03, 0.99 (U₄) |
| P5 rigid spectral flow, entry at μ = 0 (4dc82af) | **held**: continuing member within 10⁻³ at four branch points |
| P6 modes localized at the front (4dc82af) | **held**: peak at the dip, x ≈ 0.61 |
| λ₃ … λ₆ from the lower profiles (stage 2) | errors of 6 %, 2.5 %, 2 % and 0.1 % of a spacing |
| U₄ growth rates (stage 2b) | 0.2–1.2 % |
| λ₇ (stage 3: H1–H3 c922751, H4 273ae28) | **[PENDING]** |

## Figure legends
**Figure 1 | Hierarchies of self-similar profiles.**
- Axes: 1/ε_n against n for 2D Boussinesq, Hou–Luo, IPM and CCF (ε = λ − 1 or λ).
- The n-th profile has n unstable modes in every case computed.
- Open symbols: reported before. Filled: this work.
- Dashed lines: Hou–Luo and 2D continued with the asymptotic spacing (derived for Hou–Luo).
- Inset: IPM in the shifted variable. The CCF ladder ends at a sonic cusp (star).

**Figure 2 | One branch, one phase.**
- **a**, Smoothness defect |m − 2| along the continuous branch; smooth profiles are its zeros.
- **b**, The stalled layer: wall speed D/D₀ for deep profiles.
- **c**, The phase computed from the profiles alone: Re Φ₀(λ_n) − nπ for IPM (constant to ±0.007 for n = 1–5) and
  the increments for Hou–Luo.
- **d**, Hou–Luo spacings against 1/z, extrapolating to π/C = 1.279 from the λ → 1 limit problem.

**Figure 3 | The same phase quantizes the instabilities.**
- **a**, Unstable growth rates μ/ε of the n-th profile in three models: n real modes on a lattice of step ≈ ε.
- **b**, Rigid spectral flow along the Hou–Luo branch.
- **c**, IPM, blind: predicted (open) and computed (filled) spectra.
- **d**, The lattice step tends to ε.

**Figure 4 | Computing beyond the exponential wall.**
- **a**, The IPM defect |m − 2| along the branch (log scale) against the noise floors of the coarse and production
  solvers. The amplitude falls tenfold per half-period.
- **b**, The audit at the deepest profiles: the effect on m of the origin truncation, angular and radial resolution,
  far field and Newton floor.
- **c**, The phase observable at the same depth, insensitive to all of these.

**Figure 5 | How the stalled layer ends decides the hierarchy.**
- **a**, Fixed layer (Boussinesq, Hou–Luo): the front converges and the spacings increase to π/C.
- **b**, Receding layer (IPM): anatomy of the phase. The inner layer contributes a constant, the dip a bounded part,
  and the receding front all the growth.
- **c**, Sonic cusp (CCF): the defect oscillates log-periodically about p* − 2 > 0, so the hierarchy is finite.
- **d**, The three wall-speed geometries.

---

## Methods
**Self-similar equations.** In self-similar variables y = x/(1 − t)^{1+λ}, τ = −ln(1 − t):
- **2D Boussinesq:** Ω + V·∇Ω = ∂₁Θ and V·∇Θ = (λ − 1)Θ, with V = (1 + λ)y + U, U = ∇^⊥Ψ, −ΔΨ = Ω, on the half-plane
  with no penetration.
- **IPM:** V·∇R = λR with Ω = −∂₁R.
- **Hou–Luo:** the boundary restriction of Boussinesq, U′ = HΩ.

In each case the least-singular profile at given λ behaves as c|y₁|^m at the stagnation point, with
m = ε/(1 + λ − A).

**Solver A (2D, IPM).**
- Log-polar coordinates (s, β), with the exact stagnation-point structure factored out (Θ = cos^mβ Θ̂).
- March outward in s from the exact local solution at s_start = −20 (implicit Gauss–Legendre near the origin, RK4
  beyond).
- Biot–Savart by eighth-order differences in s and Chebyshev collocation in β (Nb = 32; h_s = 0.025 and 0.0125).
- Newton–Krylov on X − BS(march(X)) with central-difference Jacobian products (quadratic convergence to 10⁻¹³).
- Continuation in λ, and secant iteration on m(λ) = 2.

**Solver B** (2D check) is a global sparse Newton method on a different grid and representation.

**Hou–Luo:** a log-grid march with the Mellin multiplier of Ω ↦ U/ξ (N = 65,536–131,072).

**Error audit.**
- Origin truncation s_start ∈ {−14, …, −30}.
- Nb ∈ {32, 48, 64}.
- h_s ∈ {0.025, 0.0125, 0.00625}.
- s_max ∈ {100, 130}.
- Newton tolerance down to the round-off floor.
- Scans on fixed grids of z = 1/λ at matched settings locate deep profiles below the scatter of single solves.
- All runs are logged (`PREDICTIONS_IPM.md`, stage 3).

**Linear stability.**
- The eigen-condition is that the linearized march + Biot–Savart map T_μ has eigenvalue 1. It is validated on the
  exact time-translation mode μ = 1 (relative residual ≤ 5 × 10⁻⁸).
- Real modes: parity of the number of real eigenvalues ν > 1 of T_μ.
- Complete counts: the argument principle for det(I − T_μ) on rectangles in the right half-plane, with adaptive
  refinement and a deeper origin truncation.
- Large |Im μ| is excluded by a spectral radius of T_μ below one, decaying as about 2/|Im μ|.

**Local eigenproblem and phase.** Near a wall point in the stalled layer we write the perturbation as
r̃(Y)e^{iϕ}, with Y = y/(εx) and ϕ = ε⁻¹∫κ ds, and keep leading order in ε.
- **IPM:**

      [iκ(D̂ + ĉY) + μY∂_Y] r̃ = −G∂_Yψ̃ + iκR_yψ̃,    −∂_Y²ψ̃ + κ²ψ̃ = −iκr̃,

  with ψ̃(0) = 0 and no growing mode, and ĉ = G.
- **Boussinesq:** the same, plus a vorticity transport equation.
- **Root:** found by shooting, and tracked along the wall from the Hou–Luo closed form
  κ_HL = (−1 + √(1 + 4iĉ))/(2iD̂).
- **Phase:** Φ₀ = D₀⁻¹∫κ ds up to the front cut-off D/D₀ = 2, interpolated between grid points.
- **Exact scaling (IPM):** Y → Y/D̂, κ → κD̂, ψ̃ → ψ̃/D̂ leaves the problem invariant, so κ = K(G, μ, R_y)/D̂.
  Numerically D̂κ is constant to four digits for D̂ from 0.6 to 0.012.

**Hou–Luo limit problem.**
- Stalled layer: odd Chebyshev-weighted expansion Ω = Σ b_j T_j/√(1 − ξ²) (160–240 modes).
- Outer region: graded and logarithmic grids to 10⁸, with principal-value quadrature with local subtraction.
- Solution: Newton–Krylov continuation in the shed vorticity Ω_f.
- Phase: C = 2 Re∫κ₀ dξ/ξ.

**Blind protocol.**
- `PREDICTIONS_IPM.md` (P1–P6) was committed (git 4dc82af) before any IPM continuation.
- Each numerical prediction was committed before the computation it concerns, by a fixed rule (stage 2:
  `ipm_stage2_predict.py`; stage 3: `ipm_stage3_predict.py`, H4 from the phase).
- Failed predictions are kept.

**Data availability.** All profiles (Newton residual ≤ 10⁻¹³), spectra, branch scans and logs are in the
repository. **[A DOI-stamped archive (Zenodo) will be created on submission.]**

**Code availability.** All solvers and analysis scripts are in the repository, with the command and run time for
every number (`programs/P03_boussinesq_ladder/README.md`). **[Zenodo DOI and, for Nature Computational Science,
a Code Ocean capsule on submission.]**

## References (to be completed)
1. Luo, G. & Hou, T. Y. Potentially singular solutions of the 3D axisymmetric Euler equations. PNAS 111, 12968 (2014).
2. Chen, J. & Hou, T. Y. Stable nearly self-similar blowup of the 2D Boussinesq and 3D Euler equations with smooth
   data I: analysis. arXiv:2210.07191; II: rigorous numerics (2023).
3. Wang, Y. et al. Discovery of unstable singularities. arXiv:2509.14185 (2025).
4. Wang, Y., Léger, T., Lai, C.-Y. & Buckmaster, T. Resolving sharp gradients of unstable singularities to machine
   precision via neural networks. arXiv:2511.22819 (2025).
5. Córdoba, A., Córdoba, D. & Fontelos, M. A. Formation of singularities for a transport equation with nonlocal
   velocity. Ann. Math. 162, 1377 (2005).
6. Eggers, J. & Fontelos, M. A. Singularities: Formation, Structure, and Propagation (CUP, 2015).
