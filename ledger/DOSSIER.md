# DOSSIER — Fluids Parallel Discovery Program

**Note (October 2026).** This is the program dossier, written for the P02 result (the CCF ladder). The program's main
result is now the IPM paper prepared for *Nature Computational Science*, in `ncs-submission/`. The 2D Boussinesq
and Hou–Luo computations summarized in earlier versions of this file are its validation (Supplementary Note 1 of
that paper) and are documented in `programs/P03_boussinesq_ladder/`.

**P02 result (this document):** *The ladder of unstable self-similar singularities of the
Córdoba–Córdoba–Fontelos equation is finite. All known profiles lie on one connected branch, which ends in a
universal log-periodic square-root cusp.*

Status keys used below: VERIFIED = reproduced by ≥2 independent numerical routes and converged in all
discretisation parameters; DERIVED = analytic argument, numerically confirmed; PRELIMINARY = single route.
All numbers below come from files in this repository (paths given in N/O).

---------------------------------------------------------------------------------------------------------

## A. Initial portfolio
- 69 first-generation ideas across six categories (IDEA001–IDEA085 with gaps; `ledger/IDEAS.md`):
  - A — new fluid physics (24)
  - B — methods (12)
  - C — sci-ML (6)
  - D — mathematical fluids (12)
  - E — computational discovery (3)
  - F — cross-PDE principles (12)
- Further ideas were added from literature reconnaissance (IDEA101–111, T1–T10), plus 7 second-generation ideas spawned by P02 (IDEA112–118).
- Compute envelope: 4 Xeon cores, 15 GB RAM, no GPU. This steered the program towards 1D/2D problems where a single core can decide a question.

## B. Literature screen
- **First screen:** seven parallel death audits (A transition and stability, B 2D turbulence, C singularities, D mixing and transport, E methods and sci-ML, F complex and geophysical flows, G broad reconnaissance).
  - About 350 searches; arxiv.org and publisher pages were unreachable, so some verdicts rest on abstracts or memory. These are flagged in `literature/AUDIT_*.md`.
  - Outcome: 18 ideas KILLED as known, trivial or false; 15 PAUSED; 1 Tier A; 6 Tier B; 7 Tier C.
- **Decisive finding for P02** (AUDIT_C): Wang et al. (arXiv:2509.14185) found the CCF stable profile, the first unstable profile (λ₁ ≈ 0.6057; Eggers–Fontelos 2020) and a new second unstable profile (λ₂ ≈ 0.4703).
  - They gave empirical λ_n laws for IPM and Boussinesq, but called a CCF law "premature".
  - Wang–Léger–Lai–Buckmaster (arXiv:2511.22819) searched λ ∈ [0.455, 0.4713] for a third CCF profile and found "only a non-smooth signal at the origin".
  - Whether the CCF ladder continues was therefore explicitly open.
- **Second audit** on the final claims: `literature/AUDIT_H_second_death_audit_P02.md`. No kill.

## C. Active portfolio (final)
| Tier | Program | State |
|---|---|---|
| A | **P02 unstable-singularity ladders (CCF)** | DISCOVERY → VERIFIED (main result) |
| B | P05 vortex merger as slow passage through a saddle-node | ACTIVE: a_d/b₀ ≈ 0.209 + 5.8 Re^{−0.60} (3 Re values; prediction exponent 2/3) |
| B | P03 Couette threshold γ; P06 max Lagrangian chaos; P07 certified shadowing; P08 thermal-noise rupture; P04 integrability selection | queued (not run: compute focused on P02) |
| C | P09–P15 | cheap probes pending |
| P02 spin-offs | fractional-CCF ladder unfolding (R016); gCLM mechanism comparison (R013, R017) | PRELIMINARY / PAUSED |

## D. Parallel experimentation history (condensed; full log in `ledger/RLOG.md`)
1. **R001–R002.** Infrastructure. The 2D NS and channel solvers were verified to ≥ 8 digits: Orr–Sommerfeld c = 0.2375264888 + 0.0037396706i.
2. **R003.** Seven literature audits run in parallel, then triage.
3. **R004.** P01 (2D plane-Poiseuille chaos) launched and KILLED when audit A showed the core question was answered.
4. **R005, R010.** P05 merger law at three Reynolds numbers, running in parallel with P02.
5. **R006–R008.** P02 classical solver built.
   - A log-variable exact Hilbert multiplier, integral form, and uniform FFT and adaptive mapped grids reproduce λ₀, λ₁ and λ₂.
   - λ₂ = 0.471324227767, which corrects the earlier published 0.4703.
   - All three profiles lie on one connected branch, and the smooth profiles are its p = 2 crossings.
6. **R009.** External confirmation: arXiv:2511.22819 reports λ₂ = 0.47132422, agreeing with ours to all 8 published digits.
7. **R011.** Instability orders 0, 1, 2 for λ₀, λ₁, λ₂, from three independent diagnostics.
8. **R012.** Full branch map: the local exponent p performs a damped oscillation about 2.
9. **R013, R017.** gCLM solver validated on exact solutions; its q(c_l) scans are monotone (a different ladder mechanism).
10. **R014.** Diagnosis of why every λ-free continuation stalled near the sonic limit: a layer-translation null mode.
    - Fix 1: a geometric grid.
    - Fix 2: pinning the layer using the scaling symmetry.
    - With these, Newton converges quadratically.
11. **R015 and this dossier.** The sonic limit is a square-root cusp, and the branch approaches it with log-periodic oscillations.
    - The universal exponent was derived and confirmed.
    - p* > 2, so there is no third profile on the branch.
    - The large-λ end has no crossing either.
12. **R016.** Fractional-CCF family solver. First tracks show λ₀(s) folds near s ≈ 0.50 (PRELIMINARY).

## E. Eliminated programs (reasons)
- **KILLED as known:**
  - IDEA001: 2D PPF transient chaos; Xiao–Tao–Zhang 2021, Zhang–Tao 2022.
  - IDEA008, 011, 017, 034, 040 (Hou–Wang–Yang 2025 CAP non-uniqueness), 049, 059, 061, 069, 072.
- **KILLED as false or trivial as posed:**
  - IDEA012: Taylor-dispersion isoperimetry is false at fixed flux, since thin rectangles beat the disk.
  - IDEA048: the laminar state maximises.
  - IDEA023: a trivial ensemble baseline wins.
  - IDEA053: universality is unlikely.
- **KILLED for budget:** IDEA063.
- **PAUSED:** 15 ideas, for crowding (IDEA038, 007), feasibility (IDEA037) or resolution cost (IDEA004).
- **Within P02:**
  - The uniform-grid "λ₃ ≈ 0.57" candidate (R007) was an under-resolution artefact, eliminated by the adaptive grids.
  - The "fold at λ ≈ 0.4566" was a solver failure, not a fold.
  - gCLM scans were PAUSED because the solver fails at small c_l.

---------------------------------------------------------------------------------------------------------

## F. Main discovery
For the CCF equation θ_t + (Hθ)θ_x = 0, the self-similar profiles θ = (T−t)^λ Θ(x/(T−t)^{1+λ}) form one
connected branch parametrised by λ. The branch carries a local exponent p = λ/(1+λ+h₁), with Θ ~ ξ^p at the
origin. Smooth (analytic) profiles are exactly the points where p = 2.

1. **Exactly three smooth profiles on the branch** (VERIFIED):
   - λ₀ = 1.180777662899 (stable);
   - λ₁ = 0.6057337012 (one unstable mode, μ ≈ 0.366);
   - λ₂ = 0.471324227767 (two unstable modes, μ = 0.7229 and 0.3431).
   - The instability order equals the crossing order.
2. **The branch ends at both sides without further crossings:**
   - λ → ∞: p decreases monotonically from 2 to ≈ 1.13 at λ = 7.85 (p − 1 ≈ 1/λ). VERIFIED up to λ ≈ 8.
   - Small-λ end: after λ₂, p rises to 2.0134 and then approaches a **terminal cusp profile** with log-periodic oscillations of p at
     **λ* = 0.4535843 ± 5e−8, p* = 2.0057717 ± 5e−8** (VERIFIED).
3. **Terminal cusp** (DERIVED + VERIFIED). As the sonic depth δ = min(1+λ+HΘ/ξ) → 0:
   - the sonic factor develops a square-root cusp, den ≈ δ + k|η−η_s|^{1/2};
   - the profile develops an interior square-root singularity Θ − Θ_s ≈ B sgn(ξ−ξ_s)|ξ−ξ_s|^{1/2};
   - the amplitude is fixed locally by **B² = 2λΘ_s** (checked to 1e−4 at δ = 1.3e−6);
   - the layer width scales as **w ≈ 9.0 δ²**. The universal inner-layer problem predicts w/δ² = 4.7405/k² with no free parameter: 8.98930 against a computed 8.98925 (grid-extrapolated; 5e−6; R020/R022).
4. **Universal log-periodic approach** (DERIVED + VERIFIED). Near the cusp, the linearised profile equation is the universal nonlocal operator u' = −H[u]/(2|x|).
   - Its odd modes have purely imaginary exponents ±iτ, with **τ tanh(πτ/2) = 1/2, so τ = 0.649424**.
   - Hence p(δ) oscillates log-periodically about p*, with frequency 2τ = 1.2988 in ln δ and amplitude ∝ δ. λ(δ) carries a weaker oscillation on top of a monotone decrease to λ*, so there are no folds.
   - Fitted frequencies with the frequency left free: 1.2983–1.3009 for δ ≤ 1e−3 and δ ≤ 3e−3, from p and from λ, with or without a δ ln δ term. That is within 0.004–0.2 % of the prediction (Fig. 1c: (p − p*)/δ is a pure sinusoid in ln δ).
   - Successive extrema of p (dense sampling, 59 branch points): 2.013429 (max, δ = 1.51e−2), 2.004945 (min, δ = 1.68e−3), 2.005840 (max, δ = 1.46e−4), 2.005765 (min, δ = 1.34e−5).
   - Spacing ratios in δ are 11.52 and 10.88. The predicted alternation from the non-oscillatory term B is 11.53/10.95, and the mean prediction e^{π/2τ} is 11.23.
   - Successive amplitudes |p − p*| have ratios 0.108, 0.083 and 0.096. Their geometric mean over the last pair is 0.089, the predicted e^{−π/2τ} = 0.0890.
5. **Consequence: no third unstable CCF profile exists on the branch.**
   - p ∈ [2.004945, 2.013429] on the entire arc beyond λ₂, and p → p* > 2.
   - This explains the failed λ₃ search of arXiv:2511.22819: their window [0.455, 0.4713] lies exactly on this arc.
   - It also caps the ladder of blow-up thresholds for dissipative CCF (−Δ)^{α/2}, α_c,n = 1/(1+λ_n):
     0.4586, 0.6228, **0.6797**, and no more from this branch.

## G. Closest prior work
- **Wang et al., arXiv:2509.14185 (2025), "Discovery of unstable singularities".** They found CCF λ₀, λ₁ and λ₂ ≈ 0.4703 using PINNs plus Gauss–Newton, and gave empirical λ_n laws for IPM and Boussinesq.
  - Difference: our λ₂ = 0.471324227767 is 12-digit and classical.
  - We give the global branch structure, the termination and the finiteness.
- **Wang–Léger–Lai–Buckmaster, arXiv:2511.22819 (2025).** They obtained λ₂ = 0.47132422 to machine precision with gradient-normalised PINNs, and failed to find λ₃ in [0.455, 0.4713].
  - Difference: we explain why there is no λ₃ on the branch, and our method needs no machine learning.
- **Eggers–Fontelos, Nonlinearity 33 (2020).** They found λ₁ and the unstable-profile concept for CCF.
- **Continuous families of singular self-similar profiles in related models:**
  - Huang–Tong–Wang, arXiv:2603.25104 (gCLM with degenerate data);
  - Chen–Huang–Li, arXiv:2604.01868 (Hou–Luo, Boussinesq).
  - Difference: our CCF family is used to organise the smooth ladder, and its termination is new.
- **Hoang–Radosz, ARMA 2017.** Time-dependent cusp or needle formation in a CCF-inspired equation. This is a related phenomenon in a different object.
- **Classical spirals of branches approaching singular solutions (with infinitely many folds; here λ stays monotone and only p oscillates):** Joseph–Lundgren 1973 (the Gelfand problem) and Emden–Fowler-type problems.
  - Difference: here the operator is nonlocal (Hilbert transform), and the exponent equation τ tanh(πτ/2) = 1/2 is new.

## H. Novelty claim
1. First global solution-branch picture of CCF self-similar blow-up. All known smooth profiles are the p = 2 crossings of one connected branch, and the instability order equals the crossing order.
2. First identification of how the branch terminates: a self-similar profile with an interior square-root cusp, fixed by B² = 2λΘ_s. Its universal inner layer predicts the layer-width constant to 5e−6.
3. A derived universal law for the approach to the cusp: a log-periodic oscillation with τ tanh(πτ/2) = 1/2, confirmed numerically to 0.2 %.
4. A negative answer, for the principal branch, to the open question whether a third unstable CCF profile exists. The ladder is finite, with exactly three members.
5. A classical, ML-free numerical method that resolves sharp-layer self-similar profiles to 12 digits in seconds on a single CPU core. It follows the branch to sonic depth δ = 1.3e−6, where the layer width is 1.5e−11 in ln ξ.

## I. Method (details in `programs/P02_singularity_ladders/METHOD.md`)
- **Formulation.** Log variables η = ln ξ, with Θ = e^{cη}Ψ and β < c < 1.
  - The Hilbert transform is exact in these variables, with Mellin multiplier −tan(πa/2), a = c + ik.
  - The profile equation is solved in integral (shooting) form: φ = ln Ψ, φ(η) − φ(0) = ∫(λ/den − c).
  - 8th-order cumulative quadrature.
- **Grids.**
  - Uniform FFT grid (N up to 2¹⁸) for layer-free profiles.
  - Adaptive grids for thin layers: a tanh map, and a new geometric sinh-type map dη/ds = cosh u/(cosh u + k).
  - On both adaptive grids the Hilbert transform uses the Sidi–Israeli alternating-point rule, which is spectral for p.v. kernels.
  - Kernel differences are computed in local coordinates, to keep relative precision at centre spacings down to 1e−11.
- **Solvers.**
  - Newton–Krylov for layer-free profiles; dense LU for layers.
  - Bordered systems for three problems: smoothness (p = 2, λ unknown), arclength, and sonic depth.
- **Key device: layer pinning.** The scaling symmetry Θ → κΘ(ξ/κ) is a translation in η.
  - We use it to fix the minimum of the sonic factor at the centre of the geometric grid, with the equation den_{c+1} = den_{c−1} replacing the amplitude normalisation.
  - This removes the translation near-null mode (σ_min 2e−3) that stalled every λ-, δ- or arclength-parametrised method.
  - δ then becomes a regular parameter and can be reduced by 20–50 % per step.
- **Stability.**
  - Nonlinear eigen-condition 1 ∈ spec T_μ.
  - Linearised evolution in self-similar time.
  - Dense method-of-lines spectrum, with finite-domain loop modes removed.

## J. Derivation (sketch; full version for SI)
1. **Profile equation.** d ln Θ/dη = λ/den, where den = 1 + λ + HΘ/ξ.
   - At the origin, HΘ/ξ → h₁ = −(2/π)∫Θ/ξ², so Θ ~ ξ^p with p = λ/(1+λ+h₁).
   - Θ is analytic at the origin iff p ∈ 2ℕ. On the branch, p ∈ (1, 2.2), so p = 2 is the only possibility.
2. **Cusp.** Suppose den vanishes at ξ_s and Θ ≈ Θ_s + B sgn(ξ−ξ_s)|ξ−ξ_s|^{1/2}.
   - The identity H[sgn(y)|y|^{1/2}] = |y|^{1/2} gives den ≈ (B/ξ_s)|ξ−ξ_s|^{1/2}.
   - Inserting this in dΘ/dξ = λΘ/(ξ den) returns the same singularity with coefficient 2λΘ_s/B.
   - Self-consistency requires B² = 2λΘ_s.
   - For δ > 0 the cusp is regularised on the scale where k|x|^{1/2} ~ δ, so w ∝ δ².
3. **Linearisation near the cusp.** Write Θ = Θ₀(1+u).
   - Then u' = −λ den₁/den₀², with den₁ ≈ (Θ_s/ξ_s)H[u] and den₀² = k²|x| = (2λΘ_s/ξ_s)|x|.
   - All constants cancel, leaving **u' = −H[u]/(2|x|)**.
4. **Local exponents.** Take u = |x|^σ (even) or sgn(x)|x|^σ (odd), using H|x|^σ = −tan(πσ/2) sgn(x)|x|^σ and H[sgn|x|^σ] = cot(πσ/2)|x|^σ.
   - Even modes: σ = ½ tan(πσ/2), with roots 0 (scaling), ±½ (−½ is translation), ±2.891, …
   - Odd modes: σ = −½ cot(πσ/2), with roots **±iτ** where τ tanh(πτ/2) = ½, and ±1.830, ….
5. **Matching.** The inner layer (width w ∝ δ²) perturbs ln Θ by O(δ).
   - It excites the marginal odd modes with amplitude O(δ) and phase τ ln w = 2τ ln δ + const.
   - The inner problem itself is universal: f(X) − 1 = (1/π) p.v.∫₀^∞ F(y) X²/(y(X²−y²)) dy, with F' = 1/f. It has no free constants and gives f''(0) = 0.08900, hence w/δ² = 4.7405/k².
   - Therefore p(δ) − p* = δ[B₀ + C cos(2τ ln δ) + D sin(2τ ln δ)] + o(δ), and likewise for λ.
   - Consecutive extrema are a factor e^{π/2τ} = 11.23 apart in δ, and their amplitudes decay by the same factor.
6. **Finiteness.** Since p* ≠ 2 and the oscillation amplitude → 0, only finitely many crossings of p = 2 can occur.
   - Their number is set by the global part of the branch: here 3.

## K. Validation
| Check | Result |
|---|---|
| λ₀ vs literature γ = 0.5414465 ("6 digits") | β₀ = 0.5414479811 (lit. error 1.5e−6) |
| λ₁ | 0.6057337012 (lit. 0.6057); μ = 0.366 (lit. 0.36525) |
| λ₂ vs arXiv:2511.22819 (0.47132422) | 0.471324227767 (all 8 digits) |
| λ₂ grid-independence | 12 digits over h_s, ε, L₁ ∈ {30, 45}, L₂ ∈ {120, 170}, c ∈ [0.62, 0.78] |
| Two solvers (FFT uniform vs mapped) | λ₀ equal to 13 digits |
| Instability orders | three methods agree (0, 1, 2 unstable modes); trivial μ = 1 reproduced |
| Fixed-δ grid refinement (`verify_points.py`: 24/36/48 pts per width, h_s = 0.03/0.02) | δ = 5.4e−3: Δλ = 9e−13, Δp = 4.5e−12; δ = 1.47e−3 (the minimum of p): Δλ = 1.0e−12, Δp = 4.5e−14; δ = 1.3e−6: Δλ = 9e−10, Δp = 3.5e−11 |
| Independent runs (C vs D at 36 pts/width, h_s = 0.02; A vs C: old vs new grid construction) | agree to ≤ 1e−10 (δ ≥ 5e−3) and ≤ 1e−9 (δ ≥ 3e−4) |
| Re-gridding test at every step (run C, δ = 1e−2 … 1.3e−6) | Δλ ≤ 7e−9, Δp ≤ 3e−9 |
| w/δ² constant (cusp scaling) | 9.00–9.01 over δ ∈ [3e−4, 3e−3] (grids with ≥ 17 pts/width; coarser 6-pt estimates gave 9.10) |
| B² = 2λΘ_s: cusp coefficient, predicted k = √(2λΘ_s/ξ_s) = 0.725820 | measured 0.725749–0.725764 for |x| ∈ [3e−5, 1e−3] at δ = 1.3e−6 (1e−4); 1 % at δ = 1.7e−4 (finite-δ corrections) |
| ln Θ jump coefficient 2λ/k = 1.249855 | measured 1.249853–1.249967 (δ = 1.3e−6) |
| Oscillation frequency 2τ = 1.29885 | free fits 1.2983–1.3010 (p and λ, δ ≤ 1e−3 and ≤ 3e−3; 39 branch points) |
| Layer-width constant (universal inner problem, f''(0) = 0.08900) | predicted 4.7405/k² = 8.98930 (k² = 0.52735 on the same profile, δ = 1.47e−3) vs computed 8.99528, 8.99076 and 8.98963 at 24, 48 and 96 pts/width; Richardson limit 8.98925 (5e−6) |
| Operator accuracy on extreme grids (h_c down to 1e−11) | ‖Hm·Ψ − exact‖ ≤ 1.4e−14 |
| Large-λ end | p monotone decreasing, no crossing, up to λ = 7.85 (arclength, domain-converged) |
| Uniqueness probe (multi-start, 144 starts at 8 values of λ) | every converged start (23) reproduced the branch value of p to ≤ 1e−6 |

## L. Generalisation
- **Universality of the cusp law.** The constant ½ in u' = −H[u]/(2|x|) came from B² = 2λΘ_s. So τ is independent of the profile, of λ and of normalisation.
  - It applies to any transport equation θ_t + (Hθ)θ_x = 0 whose self-similar branch ends at a sonic point.
  - It also applies to any model whose sonic factor is regular plus a Hilbert transform of the transported field (the local structure only needs den ∝ HΘ near ξ_s).
- **Predicted cusp universality classes for fractional velocities** (DERIVED, not yet tested numerically; `cusp_exponents_s.py`). For θ_t + (HΛ^sθ)θ_x = 0:
  - the sonic-point singularity is Θ − Θ_s ∝ sgn(x)|x|^{(1+s)/2};
  - the linearisation is u' = −(a/c_s(a)) HΛ^s[u]/|x|^{1−s}, with a = (1+s)/2;
  - its odd roots are σ = s/2 ± iτ(s), with τ = 0.6494, 0.7046, 0.7594, 0.7872 at s = 0, 0.2, 0.5, 0.8;
  - so the oscillation has log-frequency ω(s) = 2τ(s)/(1−s) and amplitude ∝ δ^{1/(1−s)};
  - the s = 0 case is the CCF result above.
- **Two ladder mechanisms.**
  - Vanishing order: gCLM/CLM, where q = 1/c_l is monotone (R013). The ladder can be infinite.
  - Self-induced sonic layer: CCF, where p oscillates and terminates. The ladder is finite.
- **Fractional family θ_t + (HΛ^sθ)θ_x = 0** interpolates CCF (s = 0) and Burgers (s → 1, infinite explicit ladder λ_i = 1 + 1/i).
  - An exact Mellin-multiplier solver exists and is validated (R016).
  - First results: λ₀(s) folds near s ≈ 0.50. How the finite CCF ladder unfolds into the infinite Burgers ladder is the natural next study (PRELIMINARY).
- **The method transfers directly** to any 1D self-similar problem with nonlocal operators diagonal in Mellin space (IPM and Boussinesq 1D reductions, gSQG-type models).

## M. Limitations
1. **"Finite ladder" is proved numerically only for the connected branch containing all known profiles.** Disconnected branches cannot be excluded.
   - The multi-start probe found none, but only 23 of 144 starts converged.
2. **The arguments are asymptotic, not rigorous.** No computer-assisted proof yet (IDEA116). The derivation in J is formal matched asymptotics: log terms from resonances cannot be excluded, though fits prefer none.
3. **Large-λ end checked only to λ ≈ 8.** Beyond that, the left-boundary truncation contaminates (p → 1⁺ approaches the weight c). The trend p − 1 ≈ 1/λ is extrapolated.
4. **Instability counts come from discretised operators.** A 4th unstable mode near the cusp end is irrelevant to the claims, since no smooth profile exists there.
5. **The fractional-family and gCLM generalisations are preliminary.**
6. **The PINN comparison is only qualitative.** 1D problems favour classical methods, so no claim is made about 2D or 3D.

## N. Code (all in `programs/P02_singularity_ladders/`)
- `ccf_nk.py`: uniform FFT log-grid Newton–Krylov solver: fixed-λ, p-constrained and arclength.
- `ccf_mapped.py`: tanh-map adaptive grid with alternating-point Hilbert rule.
- `ccf_sinhgrid.py`: geometric grid in local coordinates.
- `pinned_delta.py`: layer-pinned sonic-depth continuation (main tool for the cusp end).
- `sinh_arclength.py`: translation-aware arclength.
- `dense_arclength.py`: dense Jacobian.
- `pin_refine.py`: refinement of p = 2 crossings on the pinned branch, plus instability spectrum.
- `stability_mapped.py`, `linear_evolution.py`, `block_evolution.py`, `dense_spectrum.py`: stability.
- `cusp_fit.py`: asymptotic fits. `cusp_check.py`: B² = 2λΘ_s check. `layer_shape.py`: cusp shape.
- `large_lambda_arc.py`: large-λ end. `multistart.py`: uniqueness probe.
- `fccf_nk.py`, `fccf_track.py`: fractional family.
- `gclm_nk.py` and related scripts: gCLM.
- Figures: `fig_ladder.py`, `fig_profiles.py`.

## O. Data
- **Profiles:**
  - `ladder_F16k_up_cross2.npy` (λ₀), `ladder_F16k_up_cross1.npy` (λ₁): FFT grid, [λ, φ].
  - `mapped_l2_hs0.02_eps0.025.npy` (λ₂): [η, φ].
- **Branch:**
  - `branch_map.npy`: λ ∈ [0.4625, 2], columns (λ, p, min den, layer, w, h₁).
  - `sac_S1_branch.npy`, `sac_S2_branch.npy`: arclength near λ = 0.458–0.456.
  - `pin_*_branch.npy` and logs `pin_A_final.log`, `pin_C.log`, `pin_D.log`: (δ, λ, p, layer, w, N, h_c). The post-regrid values are in the logs.
  - `large_lambda_arc.log`: λ ∈ [1.18, 11].
- **Tables:** `cusp_layer_profile.npy`, `multistart.json`, `fccf_track_l0.npy`, `fccf_track_l1.npy`.

## P. Figures
- **Fig. 1 (`fig_ladder.png`)**:
  - (a) the whole branch, p(λ), with the three smooth profiles and the cusp end;
  - (b) the oscillating end in the (λ, p) plane;
  - (c) p(δ) with the log-periodic asymptotic law;
  - (d) the square-root cusp of den.
- **Fig. 2 (`fig_profiles2.png`)**: the profiles Θ(ξ) of λ₀, λ₁ and λ₂, and the near-terminal profile at δ = 1.3e−6. It also shows the local slope λ/den and the sonic factor.
- **Fig. 3 (planned)**: the instability spectra; dense spectrum and T_μ crossings for λ₀–λ₂.
- **Fig. 4 (planned)**: method schematic. Log variables and the Mellin multiplier; geometric grid; layer pinning; cost vs accuracy.

## Q. Title
*"The unstable-singularity ladder of a nonlocal transport model is finite and ends in a log-periodic cusp"*
(alternative: *"Classical continuation reveals the finite ladder of unstable singularities in the CCF model"*)

## R. Abstract (≈190 words)
Unstable self-similar singularities are blow-up solutions reachable only from infinitely fine-tuned initial
data. They were recently computed with physics-informed neural networks for several models of fluid blow-up, and
empirical laws suggested infinite "ladders" of them. Here we show that for the Córdoba–Córdoba–Fontelos
(CCF) model of nonlocal transport the ladder is finite. An adaptive spectral continuation method in
logarithmic variables reproduces the known profiles to twelve digits in seconds on a single processor. It reveals
that the stable profile and both known unstable profiles are the three crossings of a single connected
branch of self-similar solutions through the smoothness condition p = 2 on the local exponent at the origin;
the instability order equals the crossing order. Beyond the second unstable profile, the branch develops a
self-induced sonic layer and terminates at a profile with an interior square-root cusp. We derive that the
linearisation about this cusp reduces to a universal Hilbert-transform operator, whose imaginary exponents
±iτ, τ tanh(πτ/2) = 1/2, predict log-periodic oscillations along the branch; this is confirmed numerically to 0.2%.
Because the oscillation centre has p* = 2.00577 > 2 (and p never falls below 2.00494 beyond λ₂), no third unstable profile exists on the branch.

## S. Outline (Nature Computational Science Article format)
1. **Introduction.** Unstable singularities, the PINN discovery pipeline, empirical ladders, and the open CCF question.
2. **Results:**
   - 2.1 One branch, three smooth profiles (Fig. 1a, profiles in Fig. 2);
   - 2.2 Instability order = crossing order (Fig. 3);
   - 2.3 The sonic layer and the terminal cusp (Fig. 1d, B² = 2λΘ_s, w ∝ δ²);
   - 2.4 Universal log-periodic approach (derivation, Fig. 1b,c);
   - 2.5 Consequences: no λ₃, capped dissipative thresholds, and why ML searches struggled (layers of width 1e−6 to 1e−11 at the cusp end).
3. **Discussion.** Two ladder mechanisms (vanishing order vs sonic layer); how to predict ladder length; classical vs ML discovery; outlook (fractional family, IPM and Boussinesq reductions, computer-assisted proof).
4. **Methods.** Log-variable Mellin formulation; integral form; geometric grids and the alternating-point rule; layer pinning; stability; error control.

## T. SI plan
- **S1:** full derivation of the cusp law, local exponents and matching, including the even and odd families and the resonance discussion.
- **S2:** convergence tables for λ₀–λ₂ (h_s, ε, L₁, L₂, c) and for the pinned branch (fixed-δ refinement at δ = 5.4e−3, 1.47e−3 and 1.3e−6; regrid differences; runs A/C/D/E).
- **S3:** stability methods and spectra, including identification of the loop-mode artefacts.
- **S4:** the SVD diagnosis of the translation near-null mode, and why λ-, δ- and arclength-parametrisations fail.
- **S5:** the large-λ end and the multi-start probe.
- **S6:** the fractional family (validated multiplier formula m_s(a), first tracks) and gCLM results.
- **S7:** reproducibility. Scripts, commands and runtimes: every number regenerates in under 1 h on 4 CPU cores.

## U. Three most dangerous reviewer objections
- **U1.** "Finite ladder" may be an artefact of following one branch: other smooth profiles may exist on disconnected branches.
- **U2.** The near-cusp numerics, with layers of width down to 1e−11 in η, could be contaminated by round-off or discretisation, and the log-periodic fit could be over-fitting.
- **U3.** Incremental: λ₂ was already known to 8 digits (arXiv:2511.22819), and spirals of solution branches near singular solutions are classical (Joseph–Lundgren).

## V. Responses
- **V1.**
  - We claim finiteness only for the connected branch containing all known profiles, and state it so.
  - Supporting evidence: every converged start in the multi-start probe returned the branch.
  - The large-λ end and the cusp end are both closed, so the branch is a complete arc.
  - The cusp analysis is local and applies to any branch that ends at a sonic point.
  - A deflation search and a computer-assisted proof (IDEA116) are proposed as follow-up.
- **V2.**
  - Fixed-δ grid refinement (24/36/48 points per width, two far-field spacings) agrees to 1e−12 (δ ≥ 1.5e−3) and to 1e−9 (δ = 1.3e−6).
  - Re-gridding at every step changes λ by ≤ 1e−8.
  - Kernel arguments are formed in local coordinates. We removed an earlier round-off effect (a 1.5e−6 shift in λ at h_c ~ 3e−9) and documented it.
  - The fitted frequency matches an a-priori prediction, with no free parameter, to 0.2 %.
  - The cusp coefficient k and the ln Θ jump agree with B² = 2λΘ_s to 1e−4.
  - The layer-width constant w/δ² agrees with the parameter-free inner-layer prediction to 5e−6.
  - The fit's conclusion (p* > 2) is also visible directly in the raw data: the minimum of p on the arc is 2.004945.
- **V3.**
  - The contribution is not λ₂. It is the global structure and the mechanism: why the ladder stops, where it ends, and a universal exponent with a new equation, for a nonlocal operator.
  - It is also a negative answer to an explicitly open question, plus a transferable classical method that reaches layers five orders thinner than those resolved by PINNs.
  - Unlike the Joseph–Lundgren spiral, λ stays monotone (no folds), only p oscillates. The exponent equation τ tanh(πτ/2) = ½ comes from a nonlocal (Hilbert) operator and is new.

## W. Realistic venue assessment
- **Nature Computational Science:** possible but uncertain.
  - Strengths: a topical AI-for-math context (unstable singularities); a clean negative result on an open question; a computational method with orders-of-magnitude cost advantage; derived universal law with numerical confirmation.
  - Weaknesses: a 1D model problem; the "finite ladder" is per-branch; no proof; the method is classical rather than new computational science per se. Editors may see it as specialist.
  - Estimated chance of being sent to review: moderate; of acceptance after review: low-to-moderate.
- **Better-matched venues, in order:**
  - Physical Review Letters (universal cusp law plus finite ladder);
  - Nonlinearity or SIAM J. Applied Dynamical Systems (full analysis);
  - J. Fluid Mech. (if tied to Euler/Boussinesq relevance);
  - J. Comput. Phys. (method paper).
- **To strengthen for NCS:**
  1. The same termination analysis for a 2D-relevant model: IPM or Boussinesq 1D reductions, or the fractional family.
  2. A computer-assisted proof of λ₂, or of p* > 2 bounds.
  3. A head-to-head cost/accuracy benchmark against the published PINN pipeline.
- **Nothing guarantees acceptance at any venue.** The claims above are exactly as strong as the evidence listed in K.
