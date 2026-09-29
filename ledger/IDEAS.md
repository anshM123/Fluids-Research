# IDEA REGISTRY — Fluids Parallel Discovery Program

Status codes: ACTIVE · PROMOTED · PAUSED · KILLED · MERGED · DISCOVERY · VERIFIED · CERTIFIED
Novelty verdicts (from literature death audit, `literature/AUDIT_*.md`): KNOWN · PARTIAL · OPEN

Compute envelope (measured 2026-09-28): 4 Xeon cores @2.8 GHz, 15 GB RAM, no GPU.
2D pseudo-spectral: 512² ≈ 20 RK4 steps/s, 1024² ≈ 5 steps/s; 3D 128³ ≈ 1–2 steps/s.

Each entry: **question** → *positive target* · cheapest decisive test.

---------------------------------------------------------------------------------------------------
## A. NEW FLUID PHYSICS

- **IDEA001** 2D plane-Poiseuille chaos below Re=5772: transient with super-exponential lifetimes? streamwise-localised "2D puffs" that split → DP-like transition without lift-up/streaks?
  *Target:* τ(Re) law + splitting/decay crossing Re_c in 2D. · Test: Fourier–Chebyshev 2D channel, 50 runs/Re at Lx=4π–32π.
- **IDEA003** Sharp nonlinear stability-threshold exponent γ for 2D Couette (and Kolmogorov/Poiseuille): ε_c ~ Re^{-γ}.
  *Target:* measured γ* with mechanism (echo cascade / critical layer). · Test: optimised perturbations at Re=10³–10⁵.
- **IDEA004** Re-independent energy dissipation in 2D no-slip wall collisions (dipole–wall).
  *Target:* ε(Re)→const law + mechanism. · Test: Chebyshev box, Re=10³–10⁵.
- **IDEA005** Critical slowing down at the vortex-merger threshold: T_merge ~ |d−d_c|^{-1/2} (saddle-node ghost) across gSQG-α family.
  *Target:* universal exponent + class. · Test: contour dynamics / pseudo-spectral pairs.
- **IDEA006** Law for number N of polar vortex-crystal members vs (Bu, β̂, E) (Jupiter 8/5, Saturn 1).
- **IDEA007** 2D elastic turbulence (Oldroyd-B/FENE-P Kolmogorov, Re→0) without artificial diffusion: universal spectrum/dissipation law vs Wi.
- **IDEA008** Horizontal convection ultimate regime in 2D (Nu–Ra exponent beyond 1/5).
- **IDEA009** Dipole↔bar condensate switching in stochastically forced 2D NS: Arrhenius/large-deviation rate law.
- **IDEA010** Axisymmetrisation vs long-lived tripoles of perturbed Gaussian vortices: sharp amplitude threshold ε_c(Re) law.
- **IDEA011** Kaplan–Yorke dimension of 2D Kolmogorov turbulence vs Re: D_KY ~ Re^α vs Constantin–Foias–Temam bound.
- **IDEA012** Isoperimetric problem for Taylor–Aris dispersion: is the disk the minimiser among cross-sections of fixed area?
  *Target:* theorem-like statement / new optimal shapes under extra constraints. · Test: FEM/spectral Poisson solves + shape gradient.
- **IDEA013** Caustic-formation rate of inertial particles in 2D turbulence ~ exp(−C/St).
- **IDEA014** Cyclone–anticyclone asymmetry in rotating shallow-water turbulence A(Ro,Bu) collapse.
- **IDEA015** Spontaneous stochasticity of Kelvin–Helmholtz roll-up: universal distributions.
- **IDEA016** 2D turbulence on negatively curved surfaces: inverse cascade arrested at curvature radius; new condensate law.
- **IDEA017** 2D Rayleigh–Taylor with surface tension: self-similar bubble-merger law.
- **IDEA018** Pilot-wave (walking-droplet) stroboscopic models: emergent statistical laws in lattices/corrals.
- **IDEA019** Autophoretic droplets: spontaneous-motion symmetry-breaking law.
- **IDEA057** 2D Kolmogorov flow with Rayleigh friction: first-order vs continuous condensate transition; tricritical point in (Re, α).
- **IDEA063** 2D cylinder flow at Re=10⁴–10⁶: drag-crisis analogue in purely 2D dynamics?
- **IDEA069** Onsager-vortex (negative-temperature) clustering transition in 2D GPE: critical exponents.
- **IDEA070** Odd-viscosity (chiral) 2D fluids with boundaries/compressibility: new instability or turbulence regime.
- **IDEA073** Active nematic / polar active turbulence: universality class of the transition to mesoscale turbulence.
- **IDEA074** Stokes waves near the limiting wave: law for spacing of successive super/subharmonic instability onsets.

## B. NEW COMPUTATIONAL METHODS

- **IDEA021** Deflation/deflated continuation to discover hidden steady/periodic branches in canonical 2D benchmarks (single-lid cavity, backward step, cylinder, sudden expansion).
  *Target:* previously unknown (possibly stable) branch in a canonical benchmark. · Test: Newton+deflation on 2D FD cavity, Re=500–8000.
- **IDEA022** Tensor-network (QTT/MPS) complexity of turbulence: derived & verified law χ(Re, tol) across Burgers, KS, 2D NS, 3D NS → cost O(χ³ log N) vs O(N^d).
- **IDEA023** Time-parallel (parareal-type) integration converging for *statistics* of chaotic flows.
- **IDEA024** Adjoint/AD-based UPO discovery at scale in JAX.
- **IDEA025** Casimir-exact Zeitlin (sine-bracket) truncation on T² vs pseudo-spectral+hyperviscosity for long-time 2D Euler.
- **IDEA026** Topology-aware AMR that provably preserves separation/vortex-birth events.
- **IDEA027** Certified a-posteriori error bounds for long chaotic simulations via numerical shadowing.
- **IDEA028** High-precision spectral Newton + dynamic rescaling for the full hierarchy of UNSTABLE self-similar blow-ups in 1D models (CCF, gCLM/De Gregorio, 1D Hou–Luo) to n≈20; asymptotic law for λ_n.
- **IDEA029** Walk-on-spheres (grid-free Monte Carlo) Stokes/potential-flow solver in complex geometry.
- **IDEA030** Rational-Krylov/exponential integrators removing artificial polymer-stress diffusion in viscoelastic DNS.
- **IDEA084** MPS/QTT time-stepper for 2D NS at N=2¹⁴ per direction beating DNS cost.
- **IDEA085** Interval-arithmetic enclosure of the Hopf bifurcation Reynolds number of the lid-driven cavity.

## C. SCIENTIFIC MACHINE LEARNING (only where scientifically necessary)

- **IDEA031** Data-assimilation synchronisation threshold (k_c η≈0.2 in 3D): universal cross-system law (2D enstrophy cascade, KS, shell models) with derivation.
- **IDEA032** Closure built on Mailybaev's hidden scale invariance: Re-independent LES/shell closure.
- **IDEA033** Equivariant ML search for non-Casimir conserved quantities in truncated 2D Euler.
- **IDEA034** Manifold-dimension scaling law d_M(Re) for Kolmogorov flow, d_M(L) for KS.
- **IDEA035** UPO-based prediction of extreme dissipation bursts (precursor variable).
- **IDEA049** Automated dimensionless-group discovery applied to transitional flow data.

## D. MATHEMATICAL FLUID DYNAMICS

- **IDEA037** Computer-assisted proof of chaos (horseshoe) in 2D Navier–Stokes (Kolmogorov flow).
- **IDEA038** "Speed limit" for mixing: sharp optimal exponential mix-norm decay rate under fixed enstrophy/W^{1,p}; structure of optimal stirring.
- **IDEA039** Super-exponential growth of ∇ω for 2D Euler on T² via optimised initial data.
- **IDEA040** Non-uniqueness of unforced Leray–Hopf solutions: continuation of self-similar solutions (Guillod–Šverák) + certified unstable eigenvalue.
- **IDEA041** Families of non-shear steady 2D Euler states near Kolmogorov/Poiseuille flow.
- **IDEA042** Self-similar collapse of N point vortices for large N: family counting law.
- **IDEA043** Sharp constants in enstrophy-growth inequalities.
- **IDEA044** 2D Boussinesq with fractional dissipation: critical exponent for blow-up.
- **IDEA045** Viscous Hou–Luo scenario: critical Re for singular-like growth.
- **IDEA055** Viscously-arrested singularity law: max amplification ~ ν^{−β(λ)} across 1D/2D models.
- **IDEA080** Smooth imploding self-similar profiles of compressible Euler (Merle–Raphaël–Rodnianski–Szeftel): full discrete family vs γ, stability indices.
- **IDEA081** Vorticity/entropy production law after first shock in 2D compressible Euler from smooth data.

## E. COMPUTATIONAL DISCOVERY ALGORITHMS

- **IDEA047** Quality-diversity (MAP-Elites) search over forcing/geometry space of 2D NS for novel attractors.
- **IDEA048** Ergodic optimisation of turbulence: is max time-averaged dissipation attained by a low-period UPO (Contreras-type principle for PDEs)?
- **IDEA050** Conserved-quantity-violation detectors as singularity diagnostics.

## F. CROSS-PDE PRINCIPLES

- **IDEA051** Unified law for maximal finite-time amplification exponents α in max Φ(T) ~ Φ₀^α from scaling symmetries (Burgers, KS, 2D NS, SQG, gSQG).
- **IDEA052** Super-exponential transient-chaos lifetimes as a cross-PDE universality (pipe, 2D channel, 1D models).
- **IDEA053** Minimal seeds localise into self-similar packets with universal envelope across shear flows.
- **IDEA054** Thermalisation-onset (tyger) time universality across truncated conservative PDEs.
- **IDEA056** Which invariants control long-time statistics of truncated vs untruncated 2D Euler.
- **IDEA058** Chaotic scattering of N sedimenting Stokeslets: escape-rate law vs N.
- **IDEA059** Minimum-dissipation selection principle for Saffman–Taylor finger width.
- **IDEA061** KdV–Burgers dispersive-dissipative shock envelope law.
- **IDEA064** Optimal wall-roughness shape for heat transport in 2D RB.
- **IDEA066** Maximal Lagrangian Lyapunov exponent / topological entropy of incompressible flows under energy or enstrophy constraint.
- **IDEA071** Time-dependent optimal wall-to-wall transport.
- **IDEA072** Sharp enhanced-dissipation laws for shear/cellular flows.

(Registry continues in literature/triage once audits return; second-generation ideas are appended as IDEA1xx.)

---------------------------------------------------------------------------------------------------
## STATUS TABLE (after audit A–G and first-generation experiments, 2026-09-28)

| Idea | Status | Note |
|---|---|---|
| IDEA028 (+T1–T3, IDEA044) | **PROMOTED → P02 (Tier A) → VERIFIED** | classical solver reproduces λ₀,λ₁,λ₂ of CCF to 8–12 digits; all lie on one connected branch; ladder = p(λ)=2 crossings; branch ends in a square-root cusp approached log-periodically (τ tanh(πτ/2)=½); p*=2.0058>2 ⇒ **no λ₃ on the branch (finite ladder)** — see dossier |
| IDEA005 | ACTIVE (P05, Tier B) | merger ratio vs Re measured at Re=10³, 4×10³; 1.6×10⁴ running |
| IDEA003 | ACTIVE-queued (P03) | shearing-frame solver not yet built |
| IDEA066 | ACTIVE-queued (P06) | cheap maps |
| IDEA027 | ACTIVE-queued (P07) | |
| IDEA109 | ACTIVE-queued (P08) | novelty unverified |
| IDEA002 | PAUSED (P04) | crowded (Modin–Viviani) |
| IDEA016, 101, 058, 108, 104, 021, 110 | Tier C | cheap probes pending |
| IDEA001, 008, 011, 012(as posed), 014, 017, 023, 034, 040, 048, 049, 053, 059, 061, 063, 069, 071, 072 | KILLED | see PORTFOLIO.md |
| IDEA038, 009, 007, 037, 073, 057, 004, 041, 026, 029, 103, 047, 031, 032, 111 | PAUSED | |

### Second-generation ideas spawned by P02
- **IDEA112** Spiral/ladder structure of self-similar solution families: log-periodic approach to a sonic (dead-zone) limit; universal ratio? → **VERIFIED/DERIVED (R015, R018)**: limit is a square-root cusp, not a dead zone; universal τ tanh(πτ/2)=½.
- **IDEA113** Two mechanisms for unstable-singularity ladders: degenerate vanishing order (Burgers, gCLM) vs self-induced sonic layers (CCF). Test gCLM(a), Burgers–Hilbert, fractional-velocity CCF family.
- **IDEA114** "Number of smooth self-similar profiles" as a function of a model parameter (fractional CCF / gSQG-1D): phase diagram of ladder length. → ACTIVE (R016: fractional solver validated; λ₀(s) fold near s≈0.50).
- **IDEA115** Classical adaptive spectral solvers vs PINNs for singular self-similar profiles: cost/accuracy benchmark (seconds on 1 core, 12 digits).
- **IDEA116** Computer-assisted proof (interval Newton–Kantorovich) of the λ₂ profile using the integral formulation + alternating-point quadrature error bounds. → PAUSED (next step for rigour).
- **IDEA117** Universal inner-layer problem f(X) = 1 + ½(H[F](X) − H[F](0)), F' = 1/f (fixes w/δ² ≈ 9.1 and the matching constants of the cusp spiral).
- **IDEA118** Cusp exponent family for fractional velocities HΛ^s: sonic-point singularity |x|^{(1+s)/2}, τ(s) from the Mellin symbol — test with the pinned method.

### Third-generation ideas (P03, 2D Boussinesq ladder)
- **IDEA119** The 2D Boussinesq (Hou–Luo scenario) ladder as level crossings m(λ) = 2 of ONE least-singular branch, computed with a classical log-polar Newton–Krylov solver (no PINN). → **VERIFIED** (R023–R028): λ₀…λ₇ = 1.9205593, 1.3990961, 1.2523487, 1.1842533, 1.1449857, 1.1194738, 1.1015817, 1.0883384. An independent global solver reproduces them to ≤ 1e-5. The instability index is n for all eight (two independent stability methods).
- **IDEA120** Hou–Luo 1D boundary model as a cheap proxy: same local smoothness condition (m = 2 ⇔ A = (3+λ)/2), same oscillating branch; 8+ smooth profiles accumulating at λ = 1 (spacing in 1/(λ−1) → 1.264). → ACTIVE (R025).
- **IDEA121** WKB quantisation of the quasi-stagnant boundary region as λ → 1 (D = V₁/x = O(ε), ε = 1+λ−A): complex local wavenumber (HL closed form iD̂κ² + κ − Ω/D̂ = 0; 2D vertical eigenproblem); smooth profiles ⇔ Re[K e^{iΦ}] = 0, Φ = ε⁻¹∫κ ds ⇒ 1/(λ_n−1) linear in n (explains the DeepMind empirical law) and exponentially decaying smoothness defect. → ACTIVE: spacing predicted to 0.2 % (HL). 2D: π/C between 1.476 and ≈ 1.51, depending on the correction exponent (R028). Formal; no proof claimed.
