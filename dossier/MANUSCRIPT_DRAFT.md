# The unstable-singularity ladder of a nonlocal transport model is finite and ends in a log-periodic cusp

*Draft manuscript (Nature Computational Science Article format). Numbers marked † are final values from
`programs/P02_singularity_ladders/` logs; see dossier section K for the validation table.*

## Abstract
Unstable self-similar singularities are blow-up solutions reachable only from infinitely fine-tuned initial
data. They were recently computed with physics-informed neural networks for several models of fluid blow-up,
and empirical laws suggested infinite "ladders" of increasingly unstable profiles. Here we show that for the
Córdoba–Córdoba–Fontelos (CCF) model of nonlocal transport the ladder is finite.

An adaptive spectral continuation method in logarithmic variables reproduces the known profiles to twelve
digits in seconds on a single processor. It reveals that the stable profile and both known unstable profiles
are the three crossings of a single connected branch of self-similar solutions through the smoothness
condition p = 2 on the local exponent at the origin; the instability order equals the crossing order.

Beyond the second unstable profile, the branch develops a self-induced sonic layer and terminates at a profile
with an interior square-root cusp. We derive that the linearisation about this cusp reduces to a universal
Hilbert-transform operator, whose imaginary exponents ±iτ, τ tanh(πτ/2) = 1/2, predict a log-periodic oscillation
of the branch; this is confirmed numerically to 1%. Because the oscillation centre has p* = 2.0058 > 2, no third
unstable profile exists on the branch.

## Main text

### Introduction
Whether smooth solutions of the incompressible Euler equations can blow up in finite time is a central open
problem of mathematical fluid dynamics. Much of the recent progress has come from self-similar blow-up
profiles found numerically and then confirmed by computer-assisted proofs [Chen–Hou; Wang et al.].

Stable profiles attract an open set of initial data. Unstable profiles, which require n fine-tuned conditions,
are believed to govern blow-up in more regular or more constrained settings (for example with dissipation or
boundaries). They are also the natural candidates for the Navier–Stokes problem. A 2025 study
[Wang et al., arXiv:2509.14185] used physics-informed neural networks (PINNs) with Gauss–Newton refinement to
discover families of unstable profiles in three models:
- the incompressible porous media (IPM) equation;
- the 2D Boussinesq equations with boundary;
- the 1D Córdoba–Córdoba–Fontelos (CCF) equation θ_t + (Hθ)θ_x = 0, H the Hilbert transform.

For IPM and Boussinesq the blow-up rates λ_n of the n-th unstable profile followed simple empirical laws in n,
suggesting infinite ladders. For CCF the authors found a stable profile and two unstable ones (λ₁ ≈ 0.6057,
first found by Eggers and Fontelos, and λ₂ ≈ 0.4703), and considered an asymptotic law premature. A follow-up
[arXiv:2511.22819] resolved λ₂ = 0.47132422 to machine precision with gradient-normalised networks. It
searched λ ∈ [0.455, 0.4713] for a third profile and found only solutions with "a non-smooth signal at the
origin".

Here we answer the question for CCF with classical computational mathematics. We show that the three known
profiles are the only smooth profiles on a single connected branch of self-similar solutions, identify how and
why the branch ends, and derive a universal law for its termination. The answer, "finite", is a positive
structural result: it predicts where searches will fail and why.

### Results

**One branch, three smooth profiles.** Write θ = (T−t)^λ Θ(ξ), ξ = x/(T−t)^{1+λ}. The profile equation
−λΘ + [(1+λ)ξ + HΘ]Θ' = 0 is a first-order ODE with a nonlocal coefficient:
d ln Θ/d ln ξ = λ/den(ξ), where den = 1 + λ + HΘ/ξ is the "sonic factor".

For every λ in a range there is a solution with Θ ~ ξ^p at the origin, where p = λ/(1 + λ + h₁) and
h₁ = −(2/π)∫Θ/ξ². Such a solution is analytic only if p = 2 (p ∈ 2ℕ in general; on our branch 1 < p < 2.2).
Smooth profiles are therefore level crossings p(λ) = 2 of a continuous family.

We compute this family in logarithmic variables, in which the Hilbert transform is an exact Fourier (Mellin)
multiplier, using an integral formulation and adaptive grids (Methods). The result is a single connected
branch (Fig. 1a):
- from λ → ∞, where p → 1⁺, p rises through 2 at λ₀ = 1.180777662899;
- it peaks at 2.197, returns through 2 at λ₁ = 0.6057337012 and dips to 1.926;
- it crosses again at λ₂ = 0.471324227767, reaches 2.0134 and then oscillates into a terminal point.

These are exactly the three known profiles. Our λ₂ agrees with the machine-precision PINN value to all
published digits. Linear stability analysis (three independent methods, Methods) gives 0, 1 and 2 unstable
modes (μ = 0.366; μ = 0.7229, 0.3431), so the instability order equals the crossing order along the branch.

**A self-induced sonic layer and a terminal cusp.** Beyond λ₂ the minimum δ of the sonic factor decreases
towards zero and the profile develops an internal layer. Following the branch into this regime defeats
standard continuation. The null vector of the Jacobian becomes a translation of the thin layer, to which
neither λ nor δ is sensitive at first order.

We remove this degeneracy by using the scaling symmetry of the problem, which acts as a translation in
ln ξ, to pin the layer to the centre of a geometrically graded grid. δ then becomes a regular parameter. The
branch can be followed to δ ≈ 1e−5, where the layer is 1e−9 wide in ln ξ, with Newton converging
quadratically at every step.

Near the end of the branch the sonic factor does not have a parabolic minimum. It develops a square-root
cusp, den ≈ δ + k|ln ξ − ln ξ_s|^{1/2}, and the layer width scales as w = 9.1 δ² (Fig. 1d). This is
explained by a self-consistent local structure. A square-root singularity Θ ≈ Θ_s + B sgn(ξ−ξ_s)|ξ−ξ_s|^{1/2}
is mapped by H onto B|ξ−ξ_s|^{1/2}, a square-root zero of den. Feeding that zero back through the profile ODE
returns the same singularity if and only if B² = 2λΘ_s.

The measured cusp coefficient and jump of ln Θ agree with this relation to 1%. The branch therefore ends at a
continuous profile with an interior square-root singularity at (λ*, p*) = (0.45358, 2.00577)†.

**A universal log-periodic approach.** Linearising the profile equation about the cusp, Θ = Θ_cusp(1+u), all
profile-dependent constants cancel and u obeys the universal nonlocal equation u' = −H[u]/(2|x|). Power
solutions u = |x|^σ and sgn(x)|x|^σ give two exponent families:
- even: σ = ½ tan(πσ/2), with roots 0 (scaling) and ±½ (translation);
- odd: σ = −½ cot(πσ/2), whose smallest roots are purely imaginary, σ = ±iτ with τ tanh(πτ/2) = ½, τ = 0.64942.

Matching to the inner layer, which sits at scale w ∝ δ², excites the marginal pair with amplitude O(δ) and
phase τ ln w. Hence
  p(δ) − p* = δ [B + C cos(2τ ln δ) + D sin(2τ ln δ)] + o(δ),
and the same holds for λ. The branch oscillates into the cusp in the (λ, p) plane, with successive extrema a
factor e^{π/2τ} = 11.23 apart in δ (Fig. 1b,c).

The computed branch shows exactly this. Extrema of p are 2.0134, then 2.00501 at δ ≈ 1.3e−3, then 2.00584 at
δ ≈ 1.5e−4. Fits with a free frequency return 1.28–1.29, against the predicted 2τ = 1.2988, with no adjustable
parameter in the prediction.

**Consequences.**
1. *No third unstable profile on the branch.* Beyond λ₂, p stays in (2.0050, 2.0134) and converges to p* > 2.
   The oscillation amplitude decays like δ, so p = 2 is never reached again. The branch is closed at both ends,
   so it carries exactly three smooth profiles. The "non-smooth signal at the origin" found by
   arXiv:2511.22819 in [0.455, 0.4713] is precisely this arc, where p = 2.005–2.013.
2. *A capped ladder of dissipative thresholds.* For CCF with dissipation (−Δ)^{α/2}, self-similar blow-up
   with the n-th profile is compatible only with α < 1/(1+λ_n). The thresholds on this branch are 0.4586,
   0.6228 and 0.6797, and the sequence stops.
3. *Why the problem is hard for neural solvers.* Near the terminal cusp the relevant structures have widths
   1e−6–1e−9 in ln ξ, five or more orders below what was resolved by PINNs. Classical adaptive discretisations
   handle them at negligible cost.

### Discussion
Two mechanisms generate ladders of self-similar singularities:
- *Vanishing order* (CLM, generalised CLM, Burgers). The smoothness condition is met by a local exponent that
  grows monotonically along the family (q = 1/c_l for CLM), so the ladder can be infinite.
- *Self-induced sonic layers* (CCF). The local exponent oscillates along the branch and the branch ends at a
  singular profile, so the ladder is finite. Its length is set by where the oscillation centre p* lies relative
  to the smoothness values.

The terminal law we derived is local and universal. It applies to any transport equation with Hilbert-
transform velocity whose self-similar branch reaches a sonic point, and fractional generalisations follow
from the Mellin symbol of the velocity operator. Our fractional-CCF solver, which interpolates between CCF and
the explicit infinite Burgers ladder, gives a route to mapping how finite ladders unfold into infinite ones as
the nonlocality is varied.

Methodologically, the combination of an exact Mellin multiplier, an integral formulation, geometric grading
and symmetry-based layer pinning computes unstable profiles to 12 digits in seconds on a single CPU core. It
follows singular branches through structures five orders thinner than neural solvers resolved. We see this
as complementary to machine-learning discovery: neural methods excel at finding candidate solutions in higher
dimensions, while classical continuation can then establish the global solution structure, including negative
answers.

Limitations:
- Finiteness is established for the connected branch containing all known profiles; disconnected branches
  cannot be excluded numerically.
- The asymptotic theory is formal. A computer-assisted proof of the cusp termination and of p* > 2 is a
  natural next step.

### Methods (summary; full details in SI)
- *Log-variable formulation.* η = ln ξ, Θ = e^{cη}Ψ, β < c < 1, where β = λ/(1+λ) is the far-field
  exponent. HΘ/ξ = e^{(c−1)η}(K_c ∗ Ψ) with multiplier −i tanh(π(k − ic)/2).
- *Integral form.* φ = ln Ψ satisfies φ(η) − φ(0) = ∫_0^η (λ/den − c). This avoids the spurious constraint
  that square spectral discretisations of the index-1 operator introduce. 8th-order cumulative quadrature.
- *Normalisation.* The scaling symmetry is fixed either by the far-field amplitude, or by pinning
  den_{c+1} = den_{c−1} at the grid centre.
- *Grids.*
  - Uniform FFT grids (N ≤ 2¹⁸) for layer-free profiles.
  - Geometric map dη/ds = cosh u/(cosh u + k) with centre spacing down to 1e−11, and coordinates accumulated
    by Gauss–Legendre integration to keep kernel arguments exact to 1e−15.
  - Hilbert convolution by the Sidi–Israeli alternating-point rule; operator error ≤ 1.4e−14 on all grids.
- *Solvers.* Newton–Krylov or dense LU; bordered systems for the smoothness constraint p = 2, pseudo-
  arclength and sonic depth. Re-gridding after every step, to 24–48 points per layer width.
- *Stability.* Perturbations in self-similar time: μδ = λδ − den δ_η − (Hδ/ξ)Θ_η. Three methods:
  - the eigen-condition 1 ∈ spec T_μ of the smooth-solution operator;
  - linearised evolution;
  - the dense method-of-lines spectrum, with finite-domain loop modes removed.
- *Error control.*
  - λ₂ stable to 12 digits across grid spacing, map parameters, domain lengths (L₁ ∈ {30, 45},
    L₂ ∈ {120, 170}) and weight c ∈ [0.62, 0.78].
  - The branch agrees to 1e−10 between resolutions.
  - Re-gridding changes λ by ≤ 1e−8 (Table S2).

### Figure captions
- **Fig. 1 | The connected branch of CCF self-similar profiles and its terminal cusp.**
  - (a) Local exponent p at the origin versus λ along the whole branch. Blue dots mark the smooth profiles
    (p = 2): the stable λ₀ and the unstable λ₁, λ₂. The star marks the terminal cusp.
  - (b) Beyond λ₂ the branch oscillates into the cusp without reaching p = 2 again (inset: zoom).
  - (c) p versus sonic depth δ: computed points and the asymptotic law with the predicted frequency 2τ,
    τ tanh(πτ/2) = ½.
  - (d) Square-root cusp of the sonic factor, den/δ − 1 ∝ X^{1/2} with X = |η − η_s|/w, on both flanks.
- **Fig. 2 | Profiles.** Θ(ξ) normalised by the far-field amplitude, local slope d ln Θ/d ln ξ = λ/den, and
  sonic factor, for λ₀, λ₁, λ₂ and a near-terminal profile.
- **Fig. 3 (planned) | Instability spectra of λ₀, λ₁, λ₂** from three methods.

### Key references
- Córdoba, Córdoba, Fontelos, Ann. Math. 162 (2005)
- Eggers, Fontelos, Nonlinearity 33 (2020)
- Wang et al., arXiv:2509.14185 (2025)
- Wang, Léger, Lai, Buckmaster, arXiv:2511.22819 (2025)
- Chen, Hou, Huang (Hou–Luo, gCLM; 2021–2025)
- Elgindi, Ghoul, Masmoudi, Anal. PDE 14 (2021)
- Huang, Tong, Wang, arXiv:2603.25104 (2026)
- Chen, Huang, Li, arXiv:2604.01868 (2026)
- Hoang, Radosz, ARMA (2017)
- Joseph, Lundgren, ARMA 49 (1973)
- Sidi, Israeli, J. Sci. Comput. 3 (1988)
