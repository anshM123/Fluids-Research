# Cover letters (drafts; not sent — the authors decide where and when to submit)

*Before sending: complete the [PENDING] items in the manuscript, compare our IPM λ₂–λ₄ with
Wang–Léger–Lai–Buckmaster (arXiv:2511.22819; not reachable from our environment), add authors, affiliations and
suggested referees, and create the Zenodo/Code Ocean archives.*

---

## A. Nature Physics

Dear Editor,

We submit our Article "A single phase quantizes self-similar singularities and their instabilities in
incompressible flow" for consideration in Nature Physics.

**The problem.**
- Unstable self-similar singularities of incompressible flow were recently discovered by physics-informed machine
  learning in three models: the Boussinesq equations with boundary (the standard proxy for axisymmetric Euler at a
  wall), porous-media flow and the Córdoba–Córdoba–Fontelos model.
- In each they form a hierarchy: the n-th profile has n unstable directions, and the blow-up rates lie on
  empirical lines.
- Why such hierarchies exist, whether they end, and why the instability count is n were open.

**The result.** One physical mechanism answers all three.
- A wave trapped in a stalled boundary layer accumulates a phase.
- Smooth profiles are resonances of a single branch of singular solutions, one per half-turn of that phase.
- The same phase, shifted at the stagnation point, quantizes the instability growth rates into a lattice that
  moves rigidly between profiles. This makes the n-th profile exactly n-fold unstable: an oscillation theorem for
  blow-up profiles, analogous to Sturm–Liouville nodes.

**Why the evidence should persuade.**
- **The phase is a computable physical observable.** Evaluated on the profiles alone, it locates them to 0.2 % of a
  spacing in a model on which the theory was not built.
- **The spacing constant is derived.** For the Hou–Luo model it follows from an explicit limit problem (an airfoil
  equation with a Kutta condition): π/C = 1.279, against 1.276–1.284 extrapolated from eleven profiles.
- **The test was pre-registered.** We time-stamped predictions in a public repository before computing a third
  model.
  - The spectral predictions held: index n, and growth rates to 0.2–1.2 %.
  - Our predicted accumulation point failed.
  - Explaining the failure revealed the general principle: how the stalled layer ends decides the fate of the
    hierarchy. A fixed front gives an infinite ladder with a derived spacing; a sonic cusp gives a finite one; a
    receding front gives contracting spacings.

**Why Nature Physics.** The work replaces an empirical catalogue of singularities, obtained by machine learning,
with a universal mechanism that predicts them, and it classifies how such hierarchies end. It connects WKB
quantization, oscillation theorems and singularity formation, and we expect it to interest readers in fluid
dynamics, nonlinear waves, mathematical physics and the physics of finite-time singularities.

**What we do not claim.**
- The hierarchies and the eigen-condition are derived by formal matched asymptotics with numerical verification of
  every hypothesis. They are not theorems.
- We make no claim about Euler without boundary, or about Navier–Stokes.

All code, profiles, spectra and the time-stamped prediction record are openly available.

Sincerely,
⟨corresponding author, on behalf of all authors⟩

---

## B. Nature Computational Science

Dear Editor,

We submit our Article "A single phase quantizes self-similar singularities and their instabilities in
incompressible flow" for consideration in Nature Computational Science.

**The computational problem.** Unstable self-similar singularities are among the hardest objects to compute in
fluid dynamics. Physics-informed neural networks recently found the first few in three models, but each further
profile is harder: the quantity that defines a smooth profile is exponentially small, and it shrinks tenfold per
profile. In our solvers it reaches the round-off-limited noise floor within two profiles of the deepest ones known.
The same difficulty limits neural approaches; the most recent work documents residual-hidden profile errors at the
fourth unstable IPM solution.

**Our contribution is a computational framework that sidesteps this wall.**
1. **Branch continuation with resonance detection.** Instead of training a solver for each profile, we continue a
   single branch of singular solutions and find the smooth profiles as its resonances. This found new profiles in
   two models, and a second, independent solver reproduces the first to 10⁻⁵.
2. **A resonance phase observable.** It is computed from O(1)-accurate profiles through a local eigenproblem in the
   boundary layer, and it predicts where the exponentially small defect vanishes. In a model it was not built on, it
   locates the profiles to 0.2 % of a spacing. In porous-media flow it takes an exact form, Φ = ∫K ds/D, through a
   scaling symmetry of the local problem.
3. **Certified instability counts.** Argument-principle counts of a Fredholm-type determinant, with exclusion of
   large growth rates, give the exact index of every profile.
4. **A pre-registered computational experiment.** All predictions for a new model were time-stamped in a public
   repository before computation, and the failed ones are kept.
5. **An error audit for exponentially sensitive quantities.** Origin truncation, angular and radial resolution, far
   field and solver floor are quantified separately at the deepest profiles.

**The scientific payoff.**
- The framework explains why these hierarchies exist and why the n-th profile is n-fold unstable.
- It derives the spacing constant in one model.
- It shows that how a stalled boundary layer ends decides whether a hierarchy is infinite.

**Why Nature Computational Science.** The work shows a computational strategy — continuation plus a resonance
observable plus pre-registration — that turns an exponentially ill-conditioned search into a well-conditioned
one. Its ingredients transfer to other searches for unstable self-similar solutions, including those now done with
neural networks.

**Reproducibility.** Every number in the paper has a logged command and run time. Profiles, spectra and branch scans
are released, and a Code Ocean capsule and a Zenodo archive will accompany the submission.

Sincerely,
⟨corresponding author, on behalf of all authors⟩
