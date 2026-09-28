# ACTIVE PORTFOLIO & COMPUTE ALLOCATION  (v2, 2026-09-28, final state; v1 after literature death audit A–G)

**Update (P03, 2D Boussinesq ladder):** new Tier-A program (folder `programs/P03_boussinesq_ladder`; the earlier
queued "P03 Couette threshold" is renamed P16). Result: 7 smooth profiles on one branch and a WKB-quantisation
theory of the infinite ladder (dossier/P03_DOSSIER.md; RLOG R023–R026).

**Final status:** P02 → VERIFIED main result (dossier/DOSSIER.md): the CCF self-similar branch carries exactly three smooth
profiles (λ₀, λ₁, λ₂) and terminates at a square-root cusp (λ* = 0.4535845, p* = 2.005772) approached with
log-periodic oscillations of universal frequency 2τ, τ tanh(πτ/2) = ½ — no λ₃ on the branch. P05 remains the
strongest backup (3 Re points). Other Tier B/C programs were not run (compute concentrated on P02 verification).


Budget accounting: 4 cores. "CPU-h" = core-hours. Reconnaissance phase target ≤ 15% of effort.

## TIER A
### P02 — Ladders of unstable self-similar singularities in 1D fluid models  (IDEA028 + T1–T3 + IDEA044)
- **Question.** For the CCF equation θ_t + (Hθ)θ_x = 0 (and gCLM(a), 1D Hou–Luo, Burgers–Hilbert family), what is the complete ladder {λ_n} of smooth self-similar blow-up profiles, and its n→∞ law?
- **Positive target.** New profiles n ≥ 3 (not found by PINN/Gauss–Newton searches) + a derived asymptotic law λ_n → λ* (+ correction) + consequence: blow-up thresholds α_c(n) = 1/(1+λ_n) for dissipative CCF; cross-PDE principle for ladders.
- **Nearest lit.** Wang et al. arXiv:2509.14185; Wang–Léger–Lai–Buckmaster arXiv:2511.22819; Eggers–Fontelos Nonlinearity 2020; Silvestre–Vicol 2016; Huang–Tong–Wang 2026.
- **Cheap decisive test.** Classical spectral method in log-variable (exact Hilbert multiplier) reproduces λ₁=0.6057, λ₂=0.4703; then continuation finds λ₃.
- **Structure.** Continuous family in λ with local exponent p(λ)=λ/(1+λ+h1); smooth profiles at p∈2ℕ (hypothesis) → ladder = level crossings of p(λ); limit = free-boundary profile.
- **Compute.** 1D dense Newton N≤8192: minutes per solve. Total ≲ 30 CPU-h.
- **Upside.** High (hot topic, explicit open question, positive by construction). **Generality.** gCLM(a) family, Hou–Luo 1D, IPM-like models.

## TIER B
- **P05 — Vortex merger as a dynamic saddle-node** (IDEA005). Target: (a/b)_onset(Re) − (a/b)_c(∞) ∝ Re^{−2/3}, merger delay ∝ Re^{1/3}, prefactor from normal form. Test: 2D pseudo-spectral pairs Re=10³–6×10⁴ (running). ~20 CPU-h. Upside: moderate (JFM Rapids).
- **P16 (formerly P03) — Sharp 2D Couette threshold γ** (IDEA003/T5). Target: ε_c ∼ Re^{−γ*} for H^s data; test Masmoudi–Zhao 1/3. Test: shearing-frame spectral solver, echo-seeded data, Re=10³–10⁶. ~60 CPU-h. Upside: moderate (math community).
- **P06 — Maximal Lagrangian chaos per unit enstrophy R*** (IDEA066). Target: sharp constant in [0.36, 0.5] and optimal protocol structure. Test: time-periodic shear maps, adjoint/gradient optimisation. ~10 CPU-h. Upside: moderate.
- **P07 — Certified trust horizon for chaotic PDE simulations** (IDEA027). Target: computable shadowing distance for KS/2D NS; identify shadowing-breakdown events. ~30 CPU-h. Upside: moderate–high (method).
- **P08 — Thermal-noise point rupture of thin films** (IDEA109/T8). Target: new noise-dominated self-similar rupture law & crossover. ~20 CPU-h. Upside: moderate–high (novelty unverified).
- **P04 — Integrability selection off the sphere** (IDEA002). Target: P(N | invariants) law on tori/disk. ~40 CPU-h. Upside: moderate; crowded (Modin–Viviani).

## TIER C (cheap probes only)
- P09 turbulence on hyperbolic plane (IDEA016) · P10 Couette-fraction threshold for sustained 2D chaos (IDEA101, solver ready) · P11 N-Stokeslet escape rate κ(N) (IDEA058) · P12 Taylor dispersion fixed-∇p maximiser (IDEA108) · P13 tensor-network volume law (IDEA104) · P14 cavity deflation (IDEA021) · P15 SQG discrete self-similarity (IDEA110).

## PAUSED
IDEA038 (crowded 2026), 009, 007 (scoop risk), 037 (feasibility), 073, 057, 004 (resolution cost), 041, 026, 029, 103, 047, 031, 032, 111.

## KILLED (reason)
001 known (Xiao–Tao–Zhang 2021; Zhang–Tao 2022; Huang 2024) · 008 known · 011 known (arXiv:2602.08960) · 012 false as posed · 014 incremental · 017 known · 023 trivial ensemble baseline · 034 known · 040 known (Hou–Wang–Yang 2025) · 048 target trivial (laminar maximises) · 049 known · 053 universality likely false · 059 known · 061 known · 063 budget · 069 known · 071 mostly known · 072 known.

## Allocation (planned)
Recon (done) ~5%; P02 40%; Tier B 35% (P05, P06, P03, P08 first); Tier C 5%; verification reserve 15%+.
