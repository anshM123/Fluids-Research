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
- In porous-media flow the endpoint of the hierarchy is left open, and the seventh-profile holdout was not decisive
  by the direct method. Both are documented in a companion computational paper submitted to Nature Computational
  Science.

All code, profiles, spectra and the time-stamped prediction record are openly available.

Sincerely,
⟨corresponding author, on behalf of all authors⟩

---

## B. Nature Computational Science

Dear Editor,

We submit our Article "A resonance phase locates unstable self-similar singularities beyond an exponential
precision wall" for consideration in Nature Computational Science. A companion Article on the underlying physics
(phase quantization, spectral flow and the instability index) is being submitted to Nature Physics. The two papers
share no claims that depend on each other's review.

**The computational problem.**
- Unstable self-similar blow-up profiles of incompressible flow are now sought with physics-informed neural networks
  and classical solvers, and they come in hierarchies.
- We show that following such a hierarchy meets an exponential precision wall. Each profile is a zero of a
  smoothness defect whose amplitude falls about tenfold per profile.
- A solver's residual does not reveal this. Our profiles have residuals of 10⁻¹³, while their defects are decided at
  10⁻⁹.

**What the paper shows.**
1. **The wall and its anatomy.** We give a complete error audit of a converged solver at the deepest profiles of
   incompressible porous-media flow. It finds:
   - a front–grid locking error that we measure and cancel by grid shifting (∝ h_s^7.1);
   - an angular-resolution error that stops converging;
   - a round-off cost of deeper origin truncation.
2. **A pre-registered negative result.** Five predictions for the seventh profile and a decision rule were
   time-stamped in a public repository before the computation. The defect could not decide between them
   (z₇ = 7.35 ± 0.08 against a decisiveness threshold of 0.003), and we report it as not decisive, as the rule
   requires.
3. **An observable that bypasses the wall.** The resonance phase is a quadrature over O(1)-accurate profile data.
   - Under the same numerical changes it moves 10 to more than 1,000 times less than the defect.
   - It locates the seventh profile eleven times more precisely.
   - Its imaginary part predicts the decay rate of the defect, so the precision a direct method needs can be
     computed in advance.
4. **A validation lesson.** A failure of the phase's root tracker passed every residual check, and its failure
   counter was logged but not gated. It was caught by a geometric invariant and repaired with an exact scaling of
   the local problem. We document the failure in full.
5. **An honest boundary.** The deep branch develops a steepening front that a uniform grid cannot follow.
   - We leave the endpoint of the hierarchy open.
   - We register three deeper profiles, with bands and the precision a decisive test needs, as a benchmark for
     neural and classical solvers.

**Why Nature Computational Science.** The paper turns a recurring difficulty, exponentially small quantities that
decide qualitative structure, into an explicit computational lesson with a constructive remedy. Its practices are
cheap and transfer to other searches for unstable self-similar solutions:
- pre-registration with fixed decision rules;
- grid-shift measurement of periodic discretization errors;
- invariant-based tracking and geometric cross-checks.

**Reproducibility.**
- Every number has a logged command.
- The runs are restartable, and the full pre-registration record, including its corrections, is versioned.
- Profiles, spectra and scans will be released with a Zenodo archive and a Code Ocean capsule.

Sincerely,
⟨corresponding author, on behalf of all authors⟩
