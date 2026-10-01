# A resonance phase locates unstable self-similar singularities beyond an exponential precision wall

*Article draft for Nature Computational Science. It is the computational companion of the Nature Physics draft
(`P03_MANUSCRIPT_NPHYS.md`), which carries the physics (phase quantization, spectral flow, index n).*
- Length: main text ≈ 4,300 words; abstract 150 words; six figures and two tables; Methods.
- Every number comes from a script in `programs/P03_boussinesq_ladder` and its logged output (file names in square
  brackets).
- Figures: `fig_ncs.py` → `fig_ncs1.png` … `fig_ncs6.png`.
- **[FILL]** marks values still to be inserted from runs that are finishing.

---

## Abstract
Hierarchies of unstable self-similar blow-up profiles in incompressible flow are now computed by neural and
classical solvers, but how deep they can be followed is unknown.
- **The wall.** Direct computation meets an exponential precision wall: each profile is a zero of a smoothness
  defect whose amplitude falls about tenfold per profile.
- **The audit.** In incompressible porous-media flow we audit every source of error. An erratic floor of 10⁻⁹
  leaves the seventh profile unidentifiable by its defect, as a pre-registered holdout test confirmed.
- **The phase.** A resonance phase is a quadrature over O(1)-accurate profile data, and its imaginary part predicts
  how fast the defect falls. It moves 10 to more than 1,000 times less under the same numerical changes, and it locates the
  seventh profile ten times more precisely.
- **Validation.** One tracking failure passed every residual check; it was caught by a geometric invariant and
  repaired with an exact scaling. This shows how such observables must be validated.

We register three deeper profiles, with the precision a test requires, as a benchmark.

---

## Main text

Whether the incompressible Euler equations can form a singularity from smooth data is among the central open
problems of mathematical physics. Near a boundary it has been settled with computer assistance.
- Hou and Luo found explosive growth at the wall of a rotating cylinder¹.
- Chen and Hou proved stable, nearly self-similar blow-up for axisymmetric Euler and its two-dimensional proxy, the
  Boussinesq system with a boundary².

The stable profile is only the bottom of a hierarchy.
- Physics-informed neural networks found the first few *unstable* self-similar profiles of the Boussinesq
  equations, of incompressible porous-media flow (IPM) and of the Córdoba–Córdoba–Fontelos model³.
- These have since been refined towards machine precision⁴.
- Unstable profiles are the thresholds between qualitatively different evolutions. They are also the candidates
  for blow-up wherever the stable scenario is unavailable. Computer-assisted proofs require them to high accuracy².

How far down such a hierarchy can be computed, and with what certainty a computed profile is the profile it is
claimed to be, are therefore computational questions with mathematical consequences.

We show that these hierarchies carry an intrinsic obstruction, which we call an exponential precision wall, and we
give an observable that removes it.
- **The obstruction.** The profiles are isolated members of a continuous branch of singular self-similar
  solutions. Each sits where a smoothness defect vanishes, and that defect is exponentially small in the profile
  number. In IPM its amplitude falls about tenfold per profile, so each new profile needs ten times more precision.
- **The audit.** We audit the full error budget of a converged solver. Deep in the hierarchy it is dominated by
  three things that are neither small nor smooth:
  - a periodic locking of a receding front to the grid;
  - an erratic, non-convergent angular-resolution error;
  - round-off from deep origin truncation.
- **The holdout.** A holdout test, pre-registered with a fixed decision rule, shows the consequence: the seventh
  IPM profile cannot be identified by its defect.
- **The observable.** The resonance phase follows from the physical mechanism that creates the hierarchy (the
  companion paper). It is a quadrature over O(1)-accurate data, its imaginary part predicts the wall, and it moves
  10 to more than 1,000 times less than the defect under every numerical change we made.
- **Validating it.** A failure of the phase tracker passed every residual check. It was caught by a geometric
  invariant and repaired with an exact scaling of the local problem.
- **Where the method stops.** Deeper in the branch the profiles develop a steepening front that a uniform grid
  cannot follow, so we leave the endpoint of the IPM hierarchy open. We register three deeper profiles, with the
  precision a test would need, as a benchmark.

### Smooth profiles are zeros of an exponentially small defect
We seek self-similar solutions in y = x/(1 − t)^{1+λ}, with blow-up rate parameter λ (Methods).
- For every λ there is a least-singular profile. At the stagnation point on the wall its transported scalar (the
  density R in IPM) behaves as c|y₁|^m, with
  - m(λ) = ε/(1 + λ − A);
  - A the local strain;
  - ε = λ for IPM (ε = λ − 1 for Boussinesq).
- The profile is smooth exactly when m = 2.
- We follow this branch continuously in z = 1/λ, by Newton–Krylov continuation of a log-polar marching solver
  (Methods), and locate the zeros of the defect m − 2.

In IPM the defect oscillates along the branch with an amplitude that falls geometrically (Fig. 2a).
- Its extrema are 8.05 × 10⁻², 7.70 × 10⁻³, 7.71 × 10⁻⁴, 7.55 × 10⁻⁵, 7.04 × 10⁻⁶, 5.45 × 10⁻⁷ and
  5.07 × 10⁻⁸. Successive ratios are 10.0–12.9 per half-period.
- At its zeros, the smooth profiles λ₁ … λ₇, the defect crosses with slopes |dm/dz| = 5.4 × 10⁻², 5.9 × 10⁻³,
  6.4 × 10⁻⁴, 6.6 × 10⁻⁵, 6.6 × 10⁻⁶, 4.6 × 10⁻⁷ and 2.5 × 10⁻⁸ [`crossing_slopes()` in `fig_ncs.py`;
  `ipm_lambda7_final.out`].

A root located from a defect known to ±σ_m is uncertain by σ_z = σ_m/|dm/dz|.
- So any fixed error floor is overtaken: σ_z grows by a factor 10–20 per profile (Fig. 2c).
- The precision needed to locate profile n is not a matter of convergence order. It is set by the profile number.

The same mechanism is generic. A smoothness defect that is beyond all orders in a small parameter (here ε) is the
rule for hierarchies generated by a slow layer, as in the 2D Boussinesq and Hou–Luo hierarchies of the companion
paper.

The pipeline itself was validated on 2D Boussinesq (Methods).
- Two nonlinear solvers that share no numerical component agree on all eight profiles to 10⁻⁵.
- Two independent linearizations give the instability index n for n ≤ 7, and agree on 27 of 28 unstable eigenvalues
  to 5 × 10⁻⁵.

### Anatomy of the error floor
We audited every numerical parameter at the deepest profiles (Fig. 3; Table 1). Production settings are:
- angular collocation with Nb = 32;
- radial spacing h_s = 0.0125 in s = ln r;
- origin truncation s_start = −20;
- far field at s_max = 100;
- Newton tolerance near its floor.

Three findings matter beyond this problem.

**Front–grid locking.** As λ decreases, a front on the wall recedes across the radial grid.
- m is modulated with the period of one grid cell in the front's position.
- We measure the modulation directly by solving the same λ on grids shifted by δ = ¼, ½ and ¾ of a cell.
- Its half-range is 1.5 × 10⁻⁶, 1.2 × 10⁻⁸ and 8 × 10⁻¹¹ at h_s = 0.025, 0.0125 and 0.00625, i.e. it falls as
  h_s^7.1 (Fig. 3a).
- Averaging over the four shifts cancels the first three harmonics.
- On the coarsest grid it produced a spurious zero of m − 2 at z = 6.905 and displaced the sixth profile by 0.037
  in z, 5 % of a spacing (stage 2d of the record). A plausible-looking but spurious profile is the characteristic
  failure of direct location.

**Angular resolution does not converge at the floor.** At the seventh profile (z = 7.346, h_s = 0.00625), m − 2
moves:
- by −1.1 × 10⁻⁹ from Nb = 32 to 48;
- by a further −2.3 × 10⁻⁹ from 48 to 64.

At the sixth profile, Nb 48 → 64 changed m by only 9 × 10⁻¹¹. The floor is therefore erratic, not a convergent
tail that extrapolation could remove.

**Deeper truncation is not better.** Starting the march from the exact local solution at s_start leaves an offset
that decays roughly as e^{s_start} (2.1 × 10⁻⁸ at −15 to 2.2 × 10⁻¹⁰ at −19; Fig. 3c). Beyond about −20, two
things grow instead:
- the converged Newton residual: 1–2 × 10⁻¹³ at −16 to −17, 5 × 10⁻¹³ at −20, 8 × 10⁻¹³ at −24 **[FILL: −26, −30]**
  (Fig. 3d);
- a λ-independent round-off offset (−1.5 × 10⁻⁹ at −24, identical at the fifth and seventh profiles).

Together these give an erratic floor of 1–3 × 10⁻⁹ in m at the seventh profile. The Newton floor (≲ 10⁻¹²) and the
far field (s_max 100 → 160: 6 × 10⁻¹² after a path-dependent 6 × 10⁻¹⁰ at 130) are smaller.

### A pre-registered holdout that the defect cannot decide
The deepest profile the error budget permits was used as a holdout.

**Registration.** Before any fine-grid computation below the sixth profile, we committed predictions for z₇ = 1/λ₇
to a public repository (Table 2; Fig. 5):
- three rung-fit rules on the lower profiles: shifted law, geometric contraction, three-point law (git c922751);
- one rule from the resonance phase on a coarse branch (H4, 273ae28).

The decision rule, fixed before the crossing was known (f0b3f06), has three parts:
- a hypothesis is disfavoured if it misses by more than 2σ;
- the comparison counts as decisive only if σ_z ≤ 0.003, small enough to separate the two families of rung-fit
  predictions (0.007 apart) at 2σ;
- the measured value comes from matched scans at two resolutions, extrapolated if they disagree.

**Execution.**
- The committed procedure could not run as written: its h_s = 0.0125 scan had been paused to test the locking
  error. This deviation is recorded [`ipm_lambda7_analysis_committed_rule.out`].
- On the converged grid, h_s = 0.00625 with four shifts per point, the defect gives z₇ = 7.3451 ± 0.0016
  (statistical) [`ipm_lambda7_final.out`].

**Systematics.** The systematic variants at the crossing move z₇ by:
- −0.045 (Nb 48) and −0.137 (Nb 64);
- −0.021 and −0.060 (s_start −22, −24);
- +0.024 (s_max 130);
- −0.083 (all three together).

So z₇ = 7.35 ± 0.08 from the defect.

**Verdict.** Every registered prediction lies within 0.2σ, and σ_z exceeds the decisiveness threshold about
25-fold. As
the rule requires, the holdout is recorded as **not decisive**, and no hypothesis is preferred.

The negative result is the point. Five predictions spanning 0.018 in z cannot be told apart, even though each solve
reaches a Newton residual below 10⁻¹². The defect has run out of significant digits one profile after the sixth,
where σ_z was 0.002.

### A phase observable bypasses the wall
The hierarchy has a physical origin (companion paper).
- As ε decreases, the self-similar flow along the wall stalls in a layer between the stagnation point and a front
  (Fig. 1a).
- Perturbations travel through the layer as short waves, exp(i∫κ ds/D₀), with s = ln x and D₀ = ε/m.
- The local wavenumber κ(s) is the root of an ordinary-differential eigenproblem across the layer at each wall
  point (Methods).
- Smooth profiles are resonances of this wave (Fig. 1c):

      Re Φ₀(λ_n) = nπ + δ,      Φ₀ = (1/D₀) ∫ κ ds   (stagnation point → front cut-off D/D₀ = 2).        (1)

Three properties make Φ₀ a computational observable.

**It is a quadrature over O(1) data.**
- κ depends only on the local wall data of the profile: the wall speed D, its strain μ, and the density gradients
  G = ∂ₓR and R_y.
- None of these is exponentially small. Φ₀ therefore inherits the relative accuracy of the profile, not that of the
  defect.
- Fig. 1b shows the integrand at the seventh profile: the inner layer, x < 0.5, contributes 0.886 of
  I = D₀ Re Φ₀ = 1.631, and the receding outer structure contributes the rest.

**It has an exact scaling.**
- In IPM the vorticity is slaved to the density gradient, Ω = −∂ₓR. Rescaling the wall-normal coordinate by
  D̂ = D/D₀ removes the wall speed from the local problem, so κ = K(G, μ, R_y)/D̂ exactly.
- Numerically, K is invariant to 10⁻⁸ when D̂ is varied from 0.012 to 2 at fixed (G, μ, R_y) (Fig. 4d).
- K ≈ G for small G. On the approach and the dip K levels off at 0.8–1.05 for G ≳ 3, and it is lower on the front
  side (Fig. 4c). Φ₀ is therefore close to a weighted travel time ∫ds/D, a geometric property of the stalled layer.
- Along the branch, I = 0.76 T + const for the shallow profiles and 0.83 T + const for the deep ones, where
  T = ∫ds/D̂ [`ipm_endpoint_geometry.out`].

**It is stable where the defect is not.** At the seventh profile we recomputed Φ₀ on every variant state of the
audit (Fig. 4b) [`ipm_phase_variants.out`].
- The changes that move z₇ by up to 0.137 through the defect move it by at most 0.0012 through the phase.
- Nb 32 → 48 shifts Re Φ₀ once by 0.005, after which Nb 48 and 64 agree to 10⁻⁴. Meanwhile the defect keeps moving.
- s_start and s_max change Re Φ₀ by less than 5 × 10⁻⁵ (Δz < 10⁻⁵), against 0.02–0.06 through the defect.
- Grid shifts change it by ≤ 0.003 (Δz ≤ 0.0006), against ≤ 0.006 through the defect.
- Overall the phase moves 10 to more than 1,000 times less.

**Accuracy against the ladder.** δ_n = Re Φ₀(λ_n) − nπ is 2.0331, 2.0227, 2.0320, 2.0360 and 2.0301 for n = 1–5,
i.e. constant to ± 0.005 (0.2 % of a spacing), with no adjustable parameter beyond the cut-off (Fig. 4a).
- At n = 6 it is 1.974 ± 0.009, a drift of −0.057 that we carry as an uncertainty.
- At the seventh profile the phase reads z₇ = 7.354 ± 0.007, the range spanned by δ = 1.974 and 2.031. That is
  eleven times tighter than the defect and consistent with it.
- The grid dependence of Re Φ₀ at that depth is 0.010 (Δz = 0.002).

This reading is not one of the registered predictions. The registered phase rule H4 used a coarser branch and an
earlier tracker, gave 7.352, and lies inside the range.

**The imaginary part of the same quadrature predicts the wall** (Fig. 2b).
- If the defect behaves as |A| e^{−|Im Φ₀|} cos(Re Φ₀ + arg A), its envelope decays by e^{−|ΔIm Φ₀|} per
  half-period.
- From the computed phase, −π dImΦ₀/dReΦ₀ = 1.73, 1.85, 1.95, 2.05 and 2.15 e-folds per half-period at z = 3.0–6.5.
  The extrema decay by 2.30, 2.32, 2.37, 2.56 and 2.38.
- The deficit is (0.7–1.4) λ times the prediction, the size of the O(λ) corrections that the leading-order phase
  omits.
- The precision a direct method needs at profile n is therefore itself computable, from the same O(1) quadrature
  that locates the profile.

### A failure that passed every residual check
The phase requires tracking one root κ(s) of the local eigenproblem along the wall, and it failed in a way that
residuals could not reveal.
- **Where it failed.** On the deep states (z > 7.5), D̂ rises from about 0.3 at the dip to 2 at the cut-off within
  a few tracking steps, and κ changes by the same factor.
- **What the original tracker did.** It predicted each root from the previous κ, missed it, and on failure held κ
  fixed. That inflates K = κD̂ in proportion to D̂ (Fig. 6a): at z = 7.946, K reached 5.4 − 5.6i at the cut-off,
  where the continued root is 0.33 − 0.10i.
- **The size of the error.** I was 14 % high, 3.7 in Re Φ₀, more than one profile spacing.
- **Why residuals missed it.** Every number returned was a genuine root of the local problem or a held value, two
  roots lie close together on the front side, and every profile had a Newton residual below 10⁻¹².
- **Why the counter missed it.** The tracker counted its failed points: 0 on all shallow states, 2–4 on every
  deeper one. But it carried on, because on shallow states such points had been harmless.

The error was exposed by geometry.
- The old I(z) was non-monotone (1.659, 1.723, 1.746, 1.921, 1.747, 1.856, …), while the travel time T of the
  same profiles grows smoothly (Fig. 6b).
- The repair uses the exact scaling. Beyond the dip the tracker continues K rather than κ, with predictor
  K_prev/D̂ and a quarter of the step. It accepts a root only if K is continuous to 25 %, holds K (not κ) at a failed
  point, and excludes any state with a failed point.
- It reproduces the shallow values (Re Φ₀ = 11.4568 vs 11.4566 at n = 3) and removes the jump. Of 17 deep states
  one (z = 9.386, 11 failed points) is excluded.

Three lessons generalize:
- residuals certify roots, not branches;
- exact invariances provide both predictors and validators;
- failure counters must gate results, not only log them.

The failure and its repair are part of the pre-registration record.

### The deep branch: a second wall, and an open endpoint
Following the branch to z = 9.87 (λ = 0.101) on h_s = 0.0125 shows why the endpoint of the IPM hierarchy is
beyond the present method (Fig. 6c,d) [`ipm_deep_e5.out`, `ipm_travel_time.out`, `ipm_endpoint_geometry.out`].
- **The geometry.** The dip in the wall speed deepens ever more slowly: D̂_min falls from 0.36 to 0.20 between
  z = 7.47 and 9.51, at 0.107 per unit z early and 0.045 over the last unit. The dip and the front recede, x_dip
  from 0.88 to 1.00.
- **The front steepens.** max ∂ₓR rises from 7.8 at z = 7.47 to 21 at z = 9.87.
- **Resolution.** The violation of the exact identity Ω_b = −∂ₓR, our resolution indicator, grows from 3 % to
  7–12 % at h_s = 0.0125.
- **Fine-grid check.** Halving h_s at z = 8.066 and 8.786 [`ipm_deepres_z*.out`] gives:
  - D̂_min, the dip and front positions and T unchanged to 0.3 %;
  - the indicator down from 5.1 % and 9.1 % to 1.4 % and 2.5 %;
  - max ∂ₓR up by 4–6 %;
  - Re Φ₀ up by 0.037 and 0.017 (Δz ≤ 0.007);
  - m − 2 down from 2.4 × 10⁻⁸ and 3.3 × 10⁻⁷ to 1.0 × 10⁻⁹ and 2.0 × 10⁻⁹. Between these two depths, about one
    profile spacing apart, the error of the defect at fixed h_s grows fourteenfold while the signal falls about
    fifteenfold.

  So the geometry is resolved, and the slowing of the dip is real. The steepening front is what demands resolution,
  and it does so faster than the uniform grid supplies it.

**What the phase shows.** Along the branch the phase coefficient I grows linearly, with slope 0.098 per unit z and
a weak positive curvature (+0.005 ± 0.002) over 7.47 ≤ z ≤ 9.27.
- A saturating I, the fixed-layer class of Boussinesq and Hou–Luo, is excluded.
- Three continuations remain consistent with z ≤ 9.3:
  - I keeps growing linearly, so profiles accumulate at λ = 0 with λ_n ∝ n^{−1/2};
  - I diverges at a finite z_c, so profiles accumulate at λ_c = 1/z_c > 0;
  - the branch terminates at a singular front with a finite phase.
- The data do not decide between them:
  - pole and logarithmic fits to the full range prefer z_c ≈ 10–13, but with χ² ≈ 35 for 20 degrees of freedom;
  - the deep range alone allows any z_c from 9.5 to infinity;
  - linear extrapolations of D̂_min and 1/max ∂ₓR vanish near z = 14 and 12, but both decelerate.
- We therefore record the IPM endpoint class as **unresolved**. Deciding it needs a mesh that follows the front.

**Benchmark.** Within the computed range the phase still predicts profiles.
- We registered λ₈, λ₉ and λ₁₀ from (1) before any test (git 251bf1e; Table 2; `ipm_stage4_predict.py`).
- Each comes with a band that covers the possible drift of δ, the grid correction and the smoothing error.
- With each we give the precision a decisive test (σ_z ≤ 0.01) needs, from the decay rate predicted by Im Φ₀:
  - z₈ = 7.999 [7.970, 8.017] needs σ_m ≲ 2 × 10⁻¹¹;
  - z₉ = 8.628 [8.590, 8.645] needs σ_m ≲ 1 × 10⁻¹²;
  - z₁₀ = 9.201 [9.158, 9.214] needs σ_m ≲ 6 × 10⁻¹⁴, below what a double-precision solver of this kind delivers
    (our Newton residual floor is 5 × 10⁻¹³).
- The band at n = 8 contains the rung-fit predictions of stage 3 and excludes the earlier phase rule H4, which
  was made with the faulty tracker.
- The m − 2 values we had computed in this range lie at or below the floor. They are disclosed in the record and
  carry no information.

A test of λ₁₀ therefore needs extended precision as well as a front-adapted mesh.

### Discussion
**The wall is generic.** Hierarchies of unstable self-similar solutions generated by a slow layer have smoothness
defects beyond all orders in the layer parameter. Any method that locates them as zeros of the defect pays a
precision cost that grows exponentially with depth. This holds for neural solvers as much as classical ones. The
residual a solver reports does not measure this cost: our profiles have residuals of 10⁻¹³ while their defects are
decided at 10⁻⁹.

**The remedy is to change the observable, not the solver.**
- An asymptotic invariant of the mechanism, here a WKB phase, converts an exponentially ill-conditioned root search
  into a well-conditioned quadrature.
- Its imaginary part prices the wall in advance.
- Its exact symmetries supply the predictors and validators that make it robust.

The same construction applies wherever a hierarchy arises from a stalled or slow region, for example in boundary-
driven blow-up and in self-similar solutions of dispersive and kinetic equations.

**Validation practice.** We used three practices that are cheap and that we recommend.
- **Pre-register with fixed decision rules.** It turned an inconclusive computation into a reportable negative
  result instead of a post-hoc fit.
- **Measure periodic discretization errors by grid shifting.** This exposes locking that convergence studies at
  fixed grid offsets cannot see.
- **Cross-check derived observables against an independent geometric proxy.** This caught the failure that
  residuals missed.

**Limitations.**
- The phase is derived by formal matched asymptotics, not proved.
- δ drifts by 0.057 at the sixth profile, which caps the phase's accuracy at about 1 % of a spacing deep in the
  ladder.
- The IPM computations use one solver family with extensive internal variants. In 2D Boussinesq two independent
  solvers, which share no numerical component, agree to 10⁻⁵ on all eight profiles (Methods).
- The IPM endpoint is open.

---

## Methods

**Equations.** In self-similar variables y = x/(1 − t)^{1+λ}, τ = −ln(1 − t), IPM reads V·∇R = λR with
V = (1 + λ)y + U, U = ∇^⊥Ψ, −ΔΨ = Ω = −∂₁R, on the half-plane with no penetration.
- The least-singular profile is normalized at the stagnation point.
- m = λ/(1 + λ − A), with A = −∂₁U₁(0).
- The 2D Boussinesq and Hou–Luo equations are treated in the companion paper.

**Solver A (IPM and 2D Boussinesq).**
- **Coordinates.** Log-polar (s, β), s = ln r ∈ [−120, 100], with the stagnation-point structure factored out
  (R = cos^m β R̂).
- **March.** Outward in s from the exact local solution at s_start = −20: implicit Gauss–Legendre near the origin,
  RK4 beyond s_sw = 12.
- **Biot–Savart.** Eighth-order differences in s and Chebyshev collocation in β (Nb = 32). The unknown is the wall
  data X, and the fixed point X = BS(march(X)) is solved by Newton–Krylov with central-difference Jacobian products.
- **Continuation.** In z = 1/λ, with secant predictors.
- **Scans.** Matched scans on fixed z grids locate profiles. The production residual floor is 4–8 × 10⁻¹³
  (tolerance 7 × 10⁻¹³ at h_s = 0.0125 and 1.5 × 10⁻¹² at 0.00625).
- **Cost.** One solve takes 3–8 min (h_s = 0.0125) or 11–23 min (0.00625) on one core. The deep IPM campaign
  (λ₆ to λ₁₀) used about 13 core-hours of logged solver time on a four-core machine, plus about one hour of phase
  and geometry evaluations.

**Solver B (2D Boussinesq check)** is a global sparse Newton method on a different grid, with different unknowns,
Biot–Savart evaluation and treatment of λ. It reproduces all eight 2D profiles: λ₀ … λ₇ to between 1.1 × 10⁻⁷ and
9.6 × 10⁻⁶, and to ≤ 4.6 × 10⁻⁶ with matched truncation.

**Spectral solvers.** Instability indices are counted by two independent linearizations.
- **Method 1.** The linearized march and Biot–Savart map T_μ, with:
  - real-mode parity;
  - argument-principle counts of det(I − T_μ) on right-half-plane rectangles, validated on the exact
    time-translation mode μ = 1 (relative residual ≤ 5 × 10⁻⁸);
  - large |Im μ| excluded by a spectral-radius bound.
- **Method 2.** A shift-invert eigen-solve of solver B's discretization.
- **Results.** They give index n for n ≤ 7 in 2D, with counts stable under domain and grid variants, and agree on
  27 of 28 unstable eigenvalues to 5 × 10⁻⁵. In IPM, method 1 gives index n for n = 0–4.

**Error audit and grid shifting.** At z = 7.346 the production state (h_s = 0.00625) was re-solved with:
- Nb = 48 and 64;
- s_start = −22 and −24;
- s_max = 130 and 160;
- the combination of Nb 48, s_start −22 and s_max 130;
- the radial grid shifted by δh_s (s_min → −120 + δh_s, δ = ¼, ½, ¾; `ipm_shift.py`).

The same shifts were applied at z = 7.256, 7.316 and 7.376, and at h_s = 0.025 and 0.0125 (`rich7_*.out`, run by
`rich7_run*.sh`). z₇ is the root of a weighted cubic through the shift means (`ipm_lambda7_final.py`).

**Local eigenproblem and phase.** At each wall point we solve, with ψ̃(0) = 0 and no growing mode:

    [iκ(D̂ + ĉY) + μY∂_Y] r̃ = −G ∂_Y ψ̃ + iκ R_y ψ̃,      −∂_Y² ψ̃ + κ² ψ̃ = −iκ r̃,

where ĉ = G (IPM) and Y is the scaled wall-normal coordinate.
- **Root.** Found by shooting (DOP853) with a regular start at Y = 10⁻⁸ and a secant iteration on the growing-mode
  coefficient (`ipm_local_eig.py`).
- **Integration.** κ is tracked on a grid ds = 0.025 from s = −10, with the Hou–Luo closed form as a fallback
  guess. Φ₀ = D₀⁻¹∫κ ds runs to the first point beyond the dip where D̂ = 2, interpolated between grid points.
- **Repaired tracker** (`wkb_phase3`). Beyond the dip it uses step ds/4, predictor K_prev/D̂, acceptance on
  25 % continuity of K, and a held K at failed points. States with any failed point are excluded.
- **Cost.** About 30–60 s per state.

**Geometry diagnostics.**
- Travel time T = ∫_{−10}^{cut} ds/D̂, split into inner layer, approach, dip core (half depth) and front
  (`ipm_travel_time.py`, `ipm_phase_anatomy3.py`).
- Resolution indicator max|Ω_b + ∂ₓR| / max ∂ₓR over the dip-to-front window (`ipm_deep.py`,
  `ipm_deepres_check.py`).
- Class fits (`ipm_endpoint_geometry.py`).

**Pre-registration.** Every prediction was committed to the public repository before the computation it concerns,
and nothing is edited afterwards. Outcomes, deviations and corrections are appended (`PREDICTIONS_IPM.md`).
- **Stage 1** (4dc82af): six structural predictions for IPM, made before any IPM continuation.
- **Stage 2:** λ₃–λ₆ and the U₄ spectrum by fixed extrapolation rules.
- **Stage 3:** λ₇, λ₈ (c922751, 273ae28; decision rule f0b3f06).
- **Stage 4:** λ₈–λ₁₀ from the repaired phase (251bf1e), with a decision rule for future tests.

**Reproducibility.**
- **Restartable runs.** Every continuation and scan appends its state and log line after each step. A restart
  reloads the rows and the last two states for the secant predictor (`ipm_scan.py`, `ipm_deep.py`).
- **Queues.** Fine-grid queues run under `xargs -P`, with dependency waits and skip-if-done (`rich7_run*.sh`).
- **Interruptions.** Container restarts during the campaign were absorbed without loss of computed states.
- **Record.** Each figure is regenerated from logged outputs by `fig_ncs.py`.

**Data availability.** All profiles, spectra, scans and logs are archived **[Zenodo DOI on submission]**. The
pre-registration record is in the repository (`PREDICTIONS_IPM.md`).

**Code availability.** All solvers and analysis scripts **[repository URL; Zenodo DOI and Code Ocean capsule on
submission]**.

---

## Table 1 | Error audit at the seventh IPM profile (z = 7.346, production h_s = 0.00625, Nb 32, s_start −20, s_max 100; m − 2 = +3.0 × 10⁻¹¹)
| change | Δ(m − 2) | Δz₇ via defect | Δ Re Φ₀ | Δz₇ via phase |
|---|---|---|---|---|
| grid shift ¼ / ½ / ¾ | −9.2e-11 / +2.5e-11 / −1.5e-10 | −0.004 / +0.001 / −0.006 | −0.0014 / −0.0028 / −0.0023 | ≤ 0.0006 |
| Nb 48 | −1.13e-9 | −0.045 | +0.0054 | +0.0012 |
| Nb 64 | −3.42e-9 | −0.137 | +0.0054 | +0.0012 |
| s_start −22 | −5.1e-10 | −0.021 | < 5e-5 | < 1e-5 |
| s_start −24 | −1.51e-9 | −0.060 | < 5e-5 | < 1e-5 |
| s_max 130 | +5.9e-10 | +0.024 | < 5e-5 | < 1e-5 |
| s_max 160 | +6e-12 | +0.0002 | < 5e-5 | < 1e-5 |
| Nb 48 + s_start −22 + s_max 130 | −2.06e-9 | −0.083 | +0.0054 | +0.0012 |
| h_s 0.0125 → 0.00625 (z = 7.316) | — | — | +0.0105 | −0.002 |

## Table 2 | Pre-registered predictions and outcomes (IPM; git commits in the public repository)
| stage (commit) | prediction | outcome |
|---|---|---|
| 1 (4dc82af) | accumulation as λ → 0, linear law, index n, lattice of step λ_n, rigid spectral flow, front localization | accumulation and linear law **failed** (spacings contract, 0.92 per profile); index n (n = 0–4), lattice, rigid flow, localization **held** |
| 2 | λ₃ … λ₆ from lower profiles; U₄ growth rates | errors 6 %, 2.5 %, 2 %, 0.1 % of a spacing; 0.2–1.2 % |
| 3 (c922751, 273ae28; rule f0b3f06) | z₇: 7.3415 (H1), 7.3404 (H1′), 7.3343 (H2), 7.3342 (H3), 7.352 (H4) | defect: 7.35 ± 0.08 — **not decisive**; phase reading (repaired, not registered) 7.354 ± 0.007 |
| 3 | z₈: 8.0087 (H1), 8.0063 (H1′), 7.9899 (H2), 7.9917 (H3), 7.91–7.92 (H4, low confidence) | untested (defect floor) |
| 4 (251bf1e) | z₈ = 7.999 [7.970, 8.017], z₉ = 8.628 [8.590, 8.645], z₁₀ = 9.201 [9.158, 9.214]; σ_m ≲ 2e-11, 1e-12, 6e-14 for σ_z = 0.01 | open benchmark |

## Figure legends
**Figure 1 | The object and the observable** (IPM).
- **a**, Wall speed D̂ = D/D₀ along the wall for profiles of increasing depth: a stalled layer (D̂ ≈ 1), a dip,
  and a front that recedes as z = 1/λ grows. Dotted: the front cut-off of the phase.
- **b**, The phase integrand Re κ per unit s at the seventh profile, on a logarithmic x axis so that area is
  contribution, with the partition of I = D₀ Re Φ₀.
- **c**, Re Φ₀ along the branch (repaired tracker) against the levels nπ + δ (δ = 2.031). Filled: profiles located
  by the defect (n = 1–6). Open: n = 7 with its defect uncertainty.

**Figure 2 | The exponential wall.**
- **a**, |m − 2| along the branch (coarse branch and production scans) and its extrema, which fall tenfold per
  half-period, against the production error floor at λ₇.
- **b**, Decay per half-period predicted by the imaginary part of the phase, −π dImΦ₀/dReΦ₀, against the decay of
  consecutive extrema. The gap is O(λ).
- **c**, Uncertainty of z_n from the defect at a fixed floor σ_m = 10⁻⁹ (band 0.5–3.4 × 10⁻⁹; measured crossing
  slopes for n ≤ 7, projected beyond) and from the phase (offset scatter and grid error over dRe Φ₀/dz).

**Figure 3 | Anatomy of the error floor.**
- **a**, Half-range of m over four grid shifts against h_s at z = 7.26–7.38: front–grid locking ∝ h_s^7.1.
- **b**, Changes of m at the seventh profile under each audited parameter, against the changes that would move z₇
  by 0.003 and 0.3.
- **c**, Offset of m against the origin truncation s_start, ∝ e^{s_start} down to −19; the λ-independent offset at
  −24 is round-off.
- **d**, Converged Newton residual against s_start.

**Figure 4 | The phase is stable where the defect is not.**
- **a**, Offsets δ_n of the resonance condition (1): 2.031 ± 0.005 for n = 1–5, and 1.974 ± 0.009 at n = 6.
- **b**, For every audited change at the seventh profile: the shift of z₇ implied by the phase against that implied
  by the defect. Points below 2 × 10⁻⁵ are drawn at that line.
- **c**, Re K = Re κD̂ against the wall gradient G along two profiles: K ≈ G for small G, saturating near 1.
- **d**, The exact scaling: K recomputed with D̂ replaced by 0.012–2 at fixed (G, μ, R_y), at three wall points.

**Figure 5 | The seventh profile: a pre-registered holdout.** Registered predictions (commits), the measurement from
the defect (statistical and systematic bars), the decisiveness threshold of the registered rule, and the
(unregistered) phase reading.

**Figure 6 | A failure that residuals missed, and the deep branch.**
- **a**, K along the front side of the deep state z = 7.946: the original tracker holds κ, so K grows with D̂;
  the repaired tracker continues K.
- **b**, I = D₀ Re Φ₀ from both trackers against the geometric travel time (0.80 T + const). Numbers: failed points
  the original tracker counted but did not stop on.
- **c**, Steepest density gradient on the front side along the deep branch.
- **d**, Violation of the exact identity Ω_b = −∂ₓR (resolution indicator) at h_s = 0.0125, with the h_s = 0.00625
  checks.

## References (to be completed)
1. Luo, G. & Hou, T. Y. Potentially singular solutions of the 3D axisymmetric Euler equations. PNAS 111, 12968
   (2014).
2. Chen, J. & Hou, T. Y. Stable nearly self-similar blowup of the 2D Boussinesq and 3D Euler equations with smooth
   data I: analysis; II: rigorous numerics. arXiv:2210.07191, arXiv:2305.05660 (2022–2023).
3. Wang, Y. et al. Discovery of unstable singularities. arXiv:2509.14185 (2025).
4. Wang, Y., Léger, T., Lai, C.-Y. & Buckmaster, T. Resolving sharp gradients of unstable singularities to machine
   precision via neural networks. arXiv:2511.22819 (2025).
5. Córdoba, A., Córdoba, D. & Fontelos, M. A. Formation of singularities for a transport equation with nonlocal
   velocity. Ann. Math. 162, 1377 (2005).
6. Eggers, J. & Fontelos, M. A. Singularities: Formation, Structure, and Propagation (CUP, 2015).
7. Boyd, J. P. The devil's invention: asymptotic, superasymptotic and hyperasymptotic series. Acta Appl. Math. 56,
   1 (1999).
8. Nosek, B. A., Ebersole, C. R., DeHaven, A. C. & Mellor, D. T. The preregistration revolution. PNAS 115, 2600
   (2018).
