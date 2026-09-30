# Presubmission enquiry — Nature (draft; not sent)

*For the authors to edit and send. Nothing has been submitted anywhere.*

**Before sending, check the following against the literature.**
- The published IPM blow-up rates λ₂–λ₄ (Wang et al. 2025; Wang–Léger–Lai–Buckmaster 2025), against ours:
  0.3149622, 0.2415661, 0.1987300.
- Whether any of the explanations claimed here has appeared since.

Neither literature source was reachable from this environment.

---

**Title.** A common quantization mechanism organizes self-similar fluid singularities and their stability

**Authors.** ⟨names, affiliations⟩

Dear Editor,

We would like to ask whether Nature would consider an Article on why unstable self-similar blow-up profiles in
incompressible fluids come in hierarchies, and why the n-th profile has exactly n unstable directions.

**Background.**
- Whether the 3D Euler equations can develop a finite-time singularity from smooth data is one of the central open
  problems of mathematical physics.
- Near a solid boundary it has been settled by a computer-assisted proof of stable, self-similar blow-up
  (Chen & Hou), following the Luo–Hou scenario.
- Last year, physics-informed neural networks uncovered the first few *unstable* self-similar profiles in the 2D
  Boussinesq equations (the standard proxy for axisymmetric Euler with boundary), in the incompressible-porous-media
  (IPM) equation and in the Córdoba–Córdoba–Fontelos (CCF) model (Wang et al. 2025). These solutions are thought to
  be the relevant candidates wherever the stable scenario is unavailable, for example with viscosity or without a
  boundary.
- Their blow-up rates followed an empirical straight line and their instability counts rose by one per profile.
  Nobody knew why, or whether the sequence continues.

**What we found.** One phase organizes both the profiles and their instabilities.
1. **Profiles.** The smooth profiles are resonances of a single continuous branch of generically singular
   self-similar solutions. As the exponent ε of the transported scalar decreases, the flow along the boundary stalls
   between the stagnation point and a front. A wave trapped in this stalled layer accumulates a phase C/ε, and every
   half-turn produces a smooth profile. This explains the linear laws, and predicts an infinite hierarchy
   accumulating at ε = 0.
2. **Instabilities.** The same phase, shifted by πμ/ε at the stagnation point, quantizes the unstable growth rates
   μ. They form a lattice of step ε that moves rigidly with the phase along the branch, admitting one new unstable
   mode per half-turn. The index of the n-th profile is therefore n, a Sturm-type oscillation theorem for blow-up
   profiles.
3. **Spacing constant.** For the Hou–Luo model the constant follows from an explicit λ → 1 limit problem, an airfoil
   equation on the stalled layer plus transport of the shed vorticity. It gives π/C = 1.279 against 1.276–1.284
   extrapolated from eleven computed profiles.
4. **Blind test.** For IPM we wrote down and time-stamped the predictions before computing anything: infinite
   ladder, index n, spectral lattice of step λ_n, and rigid spectral flow. Each further profile was predicted from
   the lower ones before it was computed.
   - **What held.** The spectral predictions held. The index is n; the four unstable growth rates of the fifth
     profile were predicted to 0.2–1.2 %; the spectral flow between profiles was predicted to 1e-3; and the
     profile positions to 0.1–6 % of a spacing.
   - **What failed.** The predicted accumulation at λ → 0 failed. The IPM profiles obey the same linear law, but
     about a finite accumulation point λ_c = 0.037 (constant spacing to 0.3 % over six intervals), where the
     boundary flow turns sonic.
   - **What it teaches.** The quantization and the spectral structure are common to the models. The location of
     the accumulation point is set by where the boundary flow stalls.
5. **Finite versus infinite.** When the branch ends instead at a sonic cusp at finite parameter, as in CCF, the
   phase still diverges but the defect oscillates about a non-zero value. The hierarchy then stops, after three
   profiles.

**Evidence.**
- Eight 2D Boussinesq profiles, four of them new (n = 4–7), eleven Hou–Luo profiles and seven IPM profiles, of
  which n = 5, 6 are new.
- Every Boussinesq profile is reproduced by an independent solver.
- Exact instability counts (argument principle, with grid and domain variants, plus exclusion of large |Im μ|):
  n for n ≤ 7 in 2D, n ≤ 10 in Hou–Luo and n ≤ 4 in IPM, all real.
- Spectral flow at 9 (Hou–Luo), 5 (2D) and 4 (IPM, predicted blind) branch points between profiles.

**Scope.**
- The infinite hierarchy and the eigen-condition are derived by formal matched asymptotics, and every hypothesis is
  checked numerically. They are not theorems.
- We make no claim about Euler without boundary or about Navier–Stokes.
- A computer-assisted existence proof of an unstable profile is the natural next step. The profiles and spectra
  we release are accurate enough to serve as its approximate solutions.

**Why Nature.** The work turns an empirical catalogue of singularities, obtained by machine learning, into a
predictive physical mechanism. The mechanism is shared by three different fluid models, it was tested blindly on
one of them, and it explains which singularity hierarchies are infinite and which terminate. We expect it to
interest researchers in fluid dynamics, mathematical analysis, nonlinear waves and scientific machine learning.

**Suggested reviewers.** ⟨to be chosen by the authors; avoid direct collaborators⟩

Sincerely,
⟨corresponding author⟩

---

## Abstract (≤ 200 words)
Unstable self-similar singularities have recently been discovered in several incompressible-fluid models, in
hierarchies whose n-th member has n unstable directions. Here we show that a single phase organizes both the
profiles and their stability. The smooth profiles are resonances of one continuous branch of singular self-similar
solutions: a wave trapped in a stalled boundary layer accumulates a phase C/ε as the exponent ε of the transported
scalar decreases, and each half-turn produces a new profile. The same phase, shifted by πμ/ε at the stagnation
point, quantizes the instability growth rates μ into a lattice of step ε that admits one new unstable mode per
half-turn, so the n-th profile has exactly n. For the Hou–Luo model we derive the asymptotic spacing from a limit
problem (π/C = 1.279), in agreement with eleven computed profiles. We find four new 2D Boussinesq profiles. In a
blind test on the porous-media equation, the predicted spectra held to about 1 %, while the predicted
accumulation point did not: the hierarchy accumulates where the boundary flow stalls. Where the branch instead ends
at a sonic cusp, as in the Córdoba–Córdoba–Fontelos model, the hierarchy is finite. The results turn an empirical catalogue of fluid
singularities into a predictive mechanism.
