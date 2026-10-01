# One phase, three models: the general form of the quantization

This note states the one-phase theory as a single object and checks it in each model. Evidence levels:
- [F] formal (matched asymptotics);
- [N] numerical;
- [P] proved (elementary).

§5 records the IPM endpoint as far as the deep-branch computations decide it (2026-10-01): its class is unresolved.

## 1. The object
Take a one-parameter branch of least-singular self-similar profiles, with parameter λ and a stagnation point on a
wall. Along the wall the self-similar radial speed is D(x) = V₁/x, and its stagnation-point value is D₀ = ε/m. Here
ε is the exponent of the transported scalar: λ − 1 for Boussinesq and Hou–Luo, λ for IPM.

**The phase.** In the stalled layer between the stagnation point and the front, linear waves have a local
wavenumber κ(s)/D₀ (s = ln x). κ is the root of a local eigenproblem across the layer:
- Hou–Luo: closed form;
- 2D Boussinesq: `bq_local_eig.py`;
- IPM: `ipm_local_eig.py`.

The profile phase is

    Φ₀(λ) = (1/D₀) ∫ κ ds          (stagnation point → front cut-off D/D₀ = D_cut).

**The master condition [F].** Put Φ(λ, μ) = Φ₀(λ) − π μ/ε(λ). Then:
- smooth profiles (μ = 0) sit where Re Φ(λ_n, 0) = nπ + δ;
- the unstable growth rates μ_{n,k} of the profile at λ_n satisfy Re Φ(λ_n, μ) = (k + ½)π − θ₀.

The shift −πμ/ε is the Mellin phase that the growth rate adds at the stagnation point, where the local exponent
changes from m to m(1 − μ/ε). In the layer, μ does not change the wave phase: κ(μ) = κ(0) + iμ/D̂, exactly.

**Consequences.**
- a lattice of growth rates with step ε;
- rigid spectral flow between profiles;
- the n-th profile has index n, because the two lattices in Re Φ interlace.

## 2. Measured in three models [N]

| | Hou–Luo | 2D Boussinesq | IPM |
|---|---|---|---|
| profiles computed | 11 | 8 | 7 |
| Re Φ step per profile | 3.02–3.09 → π | 3.07–3.12 (cut 3) | 3.131–3.151 for n = 1–5, 3.097 at n = 6 (cut 2) |
| δ_n = Re Φ(λ_n) − nπ | — | 1.81–2.27 (cut 3), ±0.03 at larger cut | 2.023–2.037 for n = 1–5, 1.984 at n = 6 |
| lattice step of μ | 1.07ε → 1.002ε | 1.10ε → 1.003ε | 1.12ε, 1.07ε, 1.05ε (U₂–U₄) |
| index | n (n ≤ 10) | n (n ≤ 7) | n (n ≤ 4, blind) |
| spacing constant | derived: π/C = 1.279 (limit problem) | measured 1.48–1.50 | — (see §5) |

**Status of the IPM checks.**
- The IPM phase was computed after the IPM profiles were known, so the δ_n row is a check, not a prediction.
- The IPM spectrum rows were predicted blind (PREDICTIONS_IPM.md, stages 1–2).

## 3. An exact scaling in IPM [F, checked N]
The IPM local problem is

    [iκ(D̂ + ĉY) + μ_v Y∂_Y] r̃ = −G ∂_Yψ̃ + iκ R_y ψ̃,     −∂_Y²ψ̃ + κ²ψ̃ = −iκ r̃,

with ĉ = G = ∂ₓR at the wall.

**The scaling.** It is invariant under Y → Y/D̂, κ → κD̂, ψ̃ → ψ̃/D̂. So κ = K(G, μ_v, R_y)/D̂ exactly, and

    Φ₀ = ∫ K ds / D,

the integral of a local wavenumber over the unnormalized wall speed. Numerically D̂κ is constant to 10⁻⁷ for
D̂ = 0.012 → 2 at fixed (G, μ, R_y) (`ipm_K_checks.out`).

**Why only IPM.** In Boussinesq the vorticity obeys its own transport equation, whose term [1 + …]ω̃ breaks the
scaling. In Hou–Luo the closed form carries the same 1/D̂, but K depends on the scaled vorticity.

### 3.1 The IPM phase is a travel time [N]
Along the wall, the transport equation gives d ln R/ds = m/D̂, so the density rises exponentially across the
stalled layer, and so does its gradient G = ∂ₓR.

The local wavenumber K = κD̂ follows G where G is small (K = 0.099 at G = 0.100), and saturates once G ≳ 2: K ≈ 0.85–1.0
over the outer layer, peaking at the dip. (Measured on the λ₆ profile, `ipm_dip_geometry.py`.)

Hence, with I_in ≈ 0.86 the λ-independent contribution of the inner layer (x < 0.5),

    Φ₀ ≈ (1/D₀) [ I_in + K∞ ∫_outer ds/D̂ ],          K∞ ≈ 0.95.

The phase is a travel time across the stalled layer, and the accumulation law of the hierarchy is fixed by the
geometry of the layer — how far the front recedes, and the mean wall speed behind it.

**The dip does not change this.** It narrows as it deepens, so its share of the phase does not grow with it.

**Measured with the repaired tracker** (`ipm_phase_anatomy3.out`, `ipm_travel_time.out`; z = 2.1–9.3).
- **Inner layer** (x < 0.5): 0.93 → 0.88 of I = D₀ Re Φ₀.
- **Dip core** (below half depth): it narrows from width 0.29 to 0.13 in s while D̂_min falls from 0.67 to 0.22. Its
  share goes 0.35 → 0.27, then 0.31–0.34 on the deepest states.
- **Front side**: 0.10 → 0.02.
- **Approach** (from x = 0.5 to the dip core): this is where the growth is, 0 → 0.58, as the dip and front recede
  (x_dip 0.61 → 0.99, x_cut 0.83 → 1.09).
- **Total.** I grows with slope 0.075–0.082 per unit z for z = 2–7.3 and 0.098 for z = 7.5–9.3, with a mild
  positive curvature. It tracks the travel time, I = 0.76 T + c (shallow) and 0.83 T + c (deep), with
  T = ∫ds/D̂.

## 4. Where a hierarchy accumulates, and whether it is infinite
**Where.** The rungs accumulate exactly where Re Φ₀ diverges along the branch. How it diverges fixes the law:

| divergence of Φ₀ | accumulation law | where it occurs |
|---|---|---|
| C/ε, layer regular as ε → 0 (stalled-layer ending) | 1/ε_n ≈ (π/C) n + b | Boussinesq, Hou–Luo (λ_c = 1) |
| C/(λ − λ_c) | 1/(λ_n − λ_c) ≈ an + b | IPM over λ₁–λ₆ (λ_c = 0.036–0.044, exponent 0.91–1.01) |
| (λ − λ_s)^{−1/2} (a dip of fixed shape closing linearly) | (λ_n − λ_s)^{−1/2} ≈ an + b | not IPM: the dip collapses, §3.1 |
| C/λ with a layer that lengthens like ln(1/λ) | z_n ln z_n ∝ n | IPM, if the front recession is logarithmic |
| C/λ² (layer length ∝ 1/λ, I linear in z) | z_n ∝ √n | IPM over z = 2–9.3 (I linear in z; slope 0.075–0.098) |
| logarithmic in the sonic depth δ (sonic cusp) | log-periodic | CCF |

**Finite versus infinite [P + F].**
- By the Lemma (ASYMPTOTICS §2), a defect F = R(cos Θ + η) with Θ → ∞ and |η| < 1 has infinitely many zeros.
- A divergent phase is therefore necessary, but not sufficient: the defect must also be centred (η → 0).
- **Stalled layer: centred.** Every order of the formal expansion is smooth, so m − 2 is beyond all orders and
  η → 0.
- **CCF sonic cusp: not centred.** The defect oscillates about p* − 2 = 0.0058 with an amplitude ∝ δ, so η → ∞ and
  only three profiles exist.

**The general quantization coordinate.** It is Re Φ₀ itself. Laws such as 1/|λ_n − λ_c| ≈ an + b are its special
cases, valid when Φ₀ ≈ C/|λ − λ_c|.

## 5. IPM: the endpoint — unresolved (recorded 2026-10-01)
*Earlier content of this section (kept for the record):* over λ₁–λ₆ the phase seemed to grow like 1/(λ − 0.037),
and the dip (D̂_min 0.82 → 0.45) extrapolated linearly to closure near λ ≈ 0.09. Neither inference survives the deep
branch.

**What the resolved deep branch shows** (h_s = 0.0125 continuation to z = 9.87; h_s = 0.00625 checks at z = 8.066
and 8.786; `ipm_endpoint_geometry.out`, `ipm_deepres_z*.out`).
- **The dip does not close linearly.** D̂_min falls at −0.107 per unit z over z = 7.5–8.5 and at −0.045 over the last
  unit. The value is converged to 0.3 % in h_s.
- **The dip and the front recede.**
- **The front steepens.** max ∂ₓR goes from 7.8 to 21. This is what limits the resolution: the indicator
  |Ω_b + ∂ₓR|/max reaches 13 % at h_s = 0.0125 and drops about 3.6× per halving.
- **The phase.** I grows linearly in z, so the fixed-layer class (I → const) is excluded.
  - Pole and log forms fitted to the whole range prefer z_c ≈ 10–13, with χ² ≈ 35 for 20 dof.
  - The deep range alone does not constrain z_c (≥ 9.5).
- **The termination indicators** (D̂_min, 1/max ∂ₓR) extrapolate linearly to z ≈ 14 and 12, but both decelerate.

**Open.** Three continuations remain:
- linear growth of I, with λ_c = 0 and λ_n ∝ n^{−1/2};
- divergence at a finite λ_c;
- termination at a singular front.

Deciding between them needs a mesh that follows the front. The candidate λ_c ≈ 0.037 from rung fits is not
supported by the phase data, and the endpoint is no longer inferred from m − 2.

**Phase-based λ₈–λ₁₀** are registered as a benchmark (`PREDICTIONS_IPM.md`, stage 4, 251bf1e).
