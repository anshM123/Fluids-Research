# Quantization and accumulation of unstable self-similar blow-up profiles in the 2D Boussinesq equations

*Draft for a letter-length submission. Every number comes from programs/P03_boussinesq_ladder, and each value has
a log file there. The claim hierarchy is stated in §0 and is binding for the text.*

## Abstract
We numerically discover four additional unstable self-similar blow-up profiles of the two-dimensional Boussinesq
equations with a boundary, the standard proxy for axisymmetric three-dimensional Euler blow-up. This extends the
known hierarchy to instability index seven.

**Independent checks.**
- Two independent nonlinear solvers, sharing no numerical ingredient, reproduce all eight profiles to 1e-5 or better.
  For five of them the agreement is better than 1e-6.
- Two independent linearizations recover the successive instability structure: the n-th profile has n resolved
  real unstable modes (n = 0,…,7).
- Right-half-plane eigenvalue counts find no further, oscillatory instabilities.

**Structure.**
- The discrete smooth profiles lie on a continuous branch of generically singular self-similar solutions. They are
  selected by the regularity condition m(λ) = 2 on the local exponent at the stagnation point.
- As λ → 1⁺ a stalled boundary layer and a front developing into a square-root cusp generate an oscillatory inner
  problem. Its phase quantizes the hierarchy, Re Φ(λ_n) = nπ + δ + o(1) with Φ ∼ C/(λ−1), and predicts an infinite
  sequence accumulating at λ = 1.
- The unstable eigenvalues form a second, arithmetic ladder with spacing λ_n − 1, carried by eigenmodes localized
  at the same front.
- Both ladders are governed by the same phase:
  - the growth rate enters the stalled-layer phase only through a real transport factor;
  - the stagnation point adds a phase πμ/(λ−1).
- The resulting eigen-condition predicts three things: unit spacing, a ladder offset locked to the profile phase
  along the continuous branch, and one new unstable mode per half-turn of the phase, which is the index n. We
  confirm each numerically.
- The infinite sequence is an asymptotic prediction, not a theorem.

**Context.** The one-dimensional Hou–Luo model shows the same mechanism over eleven rungs. In the
Córdoba–Córdoba–Fontelos model, by contrast, the cusp forms at a finite λ and the hierarchy terminates.

## 0. What is established, supported, predicted, and not proved
| level | statement |
|---|---|
| **Numerically established** | 8 Boussinesq rungs λ₀…λ₇, found by one solver and reproduced by an independent one. Rung n has n resolved real unstable modes, and these eigenvalues are reproduced by an independent linearization. Right-half-plane counts give exactly n unstable eigenvalues for n ≤ 7 (2D) and n ≤ 10 (Hou–Luo). |
| **Strongly supported** | A continuous branch of generically singular (least-singular) profiles; the smooth profiles are its regularity-selected points. |
| **Asymptotically predicted (formal)** | An unbounded phase Φ ∼ C/(λ−1) with C > 0, the quantization Re Φ(λ_n) = nπ + δ + o(1), and an infinite ladder accumulating at λ = 1. The eigen-condition Re Φ(λ) − πμ/(λ−1) + θ₀ = (k + ½)π for the unstable spectrum, with its consequences: unit spacing, phase-locked offset, and index = number of half-turns of Φ. |
| **Not proved** | Infinitely many smooth Boussinesq profiles. C > 0 for the 2D reduction (proved only for the Hou–Luo local root; checked pointwise in 2D). The value of the asymptotic spacing. The offsets c of the instability ladder, and how a new unstable mode enters near μ = 0. |

## Evidence chain
1. Method A (marching Newton–Krylov) finds 8 profiles.
2. Method B (global sparse Newton) finds the same 8 profiles.
3. Linearization A finds 0, 1, …, 7 real unstable modes.
4. Linearization B gives the same eigenvalues.
5. Continuation shows one underlying singular branch.
6. The regularity diagnostic m(λ) = 2 picks out discrete smooth intersections.
7. Inner asymptotics give an oscillatory, quantized mechanism.
8. As λ → 1⁺ the phase grows and a cusp/front limit forms, with the unstable modes localized at that front.
9. The same phase, together with a stagnation-point connection phase πμ/(λ−1), predicts the instability ladder.
   Between the rungs the ladder moves with the phase, and one mode is added per rung interval.

## 1. Introduction (outline)
**Background.**
- The Luo–Hou scenario and stable blow-up. Chen–Hou computed the stable Boussinesq/Euler profile with computer
  assistance.
- Unstable singularities and the Navier–Stokes question; the PINN discovery of unstable Boussinesq profiles and
  their empirical law (Wang et al. 2025).
- Rigorous existence of the unstable members is open.

**Question and viewpoint.** Is the hierarchy finite or infinite, and what organizes it? Smooth profiles are
regularity-selected members of a continuous family of singular profiles, and both the profile ladder and the
instability ladder are generated at a single front.

## 2. Two independent solvers find the same eight profiles (Table 1, Fig. 1a)
- **Method A.** Log-polar marching, the exact local structure factored out, Newton–Krylov on the stream function.
- **Method B.** A global sparse Newton solve on the unhatted fields: upwind finite differences, the elliptic
  Biot–Savart law, and λ treated as an eigenvalue fixed by smoothness.
- **Agreement.** ≤ 6e-7 for λ₁–λ₅, 1.7e-6 for λ₀, a few 1e-6 for λ₆ across truncations, and 2e-7 for λ₇.

## 3. Instability structure (Table 2, Fig. 2)
**Linearizations.**
- Linearization A: the eigen-condition ν(μ) = 1 of the linearized march + Biot–Savart map T_μ.
- Linearization B: a shift-invert eigen-solve of the global discretization.

**Real modes.** Rung n has n real unstable modes. Linearization B reproduces 27 of the 28 eigenvalues to ≤ 5e-5;
the 28th, λ₇'s lowest (0.0699), is confirmed by linearization A with two origin truncations.

**Right-half-plane counts** (argument principle for det(I − T_μ)).
- 2D: n + 1 zeros (the trivial time-translation mode included) for every n ≤ 7.
- Domain variants agree: moving the origin truncation from −20 to −30 and the contour's left edge from 0.06 to
  0.075/0.07 gives 5.954 (n = 5, against 5.952) and 6.945 (n = 6, against 6.942).
- Grid variant for n = 6 (hs = 0.0125): 7.003. Grid variant for n = 7: ⟨running⟩.

**Hou–Luo, as a control.**
- Eleven rungs, each with exactly n unstable eigenvalues.
- All real up to n = 7. From n = 8 on, weakly oscillatory complex pairs appear at the top of the unstable spectrum
  (for example 0.587 ± 0.022i for n = 8). The total index stays n.

**Truncation artifacts.** Spurious crossings near μ ≈ 0.6/|s_start| move with the truncation and are discarded.

## 4. A continuous singular branch and regularity selection (Fig. 1a)
- Local analysis at the stagnation point gives Θ ≈ −|y₁|^m with m = (λ−1)/(1+λ−A). Smooth ⇔ m = 2 ⇔
  λ = −3 − 2∂₁U₁(0).
- Continuation shows one branch on λ ∈ [1.069, 4.1]. m − 2 alternates in sign, with geometrically decaying extrema;
  its zeros are the eight profiles.

## 5. The λ → 1⁺ limit: stalled layer and limiting cusp (Fig. 1c)
- On the wall, the radial self-similar speed is O(λ−1) up to a front x_c ≈ 0.72. Outside the front,
  D ∝ (x − x_c)^{1/2}: the square-root cusp that terminates the CCF branch, reached here only as λ → 1.
- A dip just inside the front deepens slowly. In Hou–Luo down to ε = 0.012 (λ − 1 ≈ 0.025): depth ∝ ε^{0.31},
  width ∝ ε^{0.34}.

## 6. Phase quantization (formal; ASYMPTOTICS.md)
**Local problem.** Perturbations exp(i∫κ ds/ε) satisfy a boundary-layer eigenproblem, and every dropped term is
O(ε). The Hou–Luo root is explicit and has Re κ > 0 whenever Ω > 0.

**Representation.**
- m − 2 = |K| e^{−Im Φ}[cos(Re Φ + arg K) + O(ε)], with Re Φ = C/(λ−1) + E and C > 0.
- The measured dip scaling makes E non-analytic. In Hou–Luo the phase coefficient converges with corrections
  ∝ ε^{2/3} (data to z ≈ 40); in 2D the correction is estimated from pre-asymptotic data.

**Elementary lemma (proved).** This representation implies infinitely many resonances with
Re Φ(λ_n) = nπ + δ + o(1).

**Checks.**
- Phase gain per resonance: 2.98 → 3.12 (2D) and 3.02 → 3.09 (Hou–Luo), tending to π.
- Integrated through the front, Re Φ(λ_n) − nπ is constant to ±0.03.

**Asymptotic spacing.**
- 2D: π/C between 1.476 and ≈ 1.51; 1.478 with analytic corrections, 1.49–1.50 with the dip-type correction.
- Hou–Luo: ≈ 1.27, against measured spacings of 1.258–1.265.
- 3/2 is not claimed.

## 7. A second ladder, and one phase for both (Fig. 2; INSTABILITY_LADDER.md §6–7)
**The ladder.**
- In both models the lower unstable eigenvalues obey μ_{n,k} ≈ (λ_n − 1)(k + c).
- Extrapolated to λ → 1, the spacing is (λ_n − 1)(1.00 ± 0.03).
- The offset at the rungs is model-dependent: c ≈ 0.70–0.72 (2D) and 0.65–0.67 (Hou–Luo), in units of the local
  spacing.
- The eigenfunctions are localized at the dip/front of the stalled layer and decay like x² toward the stagnation
  point.

**Spectral flow along the continuous branch** (Hou–Luo, ten branch points between rungs 7 and 10).
- Between the rungs the whole lower ladder moves up rigidly by one spacing per rung interval, linearly in
  z = 1/(λ−1).
- The offset equals its rung value plus the fraction of the interval, to within 0.03 at every point.
- A new member enters at the bottom once per interval. So each rung has one more unstable mode than the previous one.

**Formal theory.** Two facts combine.
1. The stalled-layer WKB root depends on the growth rate only through κ(μ) = κ(0) + iμ/D̂. This is exact for the
   Hou–Luo root. For the 2D local eigenproblem it is verified to 1e-10. The profile-ladder phase Re Φ is therefore
   the same for every μ; μ enters only through the real transport factor |Θ̄|^{−μ/(λ−1)}, the WKB form of the exact
   identity L_μ(Θ̄φ) = Θ̄ L_{μ+λ−1}φ.
2. At the stagnation point the transport is the Euler operator ε[(μ/ε − 2) + X∂_X]. Its connection to the WKB wave
   is a Mellin (Γ-function) factor with phase π(μ/(λ−1) − 1).

The eigen-condition is therefore

    Re Φ(λ) − π μ/(λ−1) + θ₀ = (k + ½)π + o(1),

with the same Re Φ that quantizes the profiles, Re Φ(λ_n) = nπ + δ. It predicts:
- a spacing of exactly λ − 1;
- the phase-locked spectral flow;
- one new unstable mode per half-turn of the phase, hence index n at rung n (the constant fixed by any one rung);
- a uniform ladder that ends near C/π = 1/Δz_∞.

The highest eigenvalue lies 0.06–0.10 above the uniform ladder, and this displacement is not described by the
condition.

The offsets c and the entry of new modes near μ = 0, where the dilation mode sits, are not derived.

## 8. Finite versus infinite ladders
- **CCF (P02).** The sonic depth vanishes at a finite λ* = 0.4535843. The branch ends in a square-root cusp and
  carries three profiles. This is consistent with the unsuccessful search for a third profile in
  λ ∈ [0.455, 0.4713] by Wang, Léger, Lai and Buckmaster (arXiv:2511.22819).
- **Boussinesq and Hou–Luo.** The cusp forms only as λ → 1. The stalled layer in front of it supports the quantized
  oscillation, so an unbounded ladder is predicted.
- **Organizing principle (conjectural).** front formation → wave quantization → profile hierarchy →
  instability-index hierarchy.

## Methods (SI)
- Solvers A and B, their resolution tables and cost: one CPU core, 20–120 s per profile.
- Linearizations A and B, the argument-principle counter (stab_contour2.py, hl_contour.py) and its integrality
  diagnostics.
- Spectra at branch points between the rungs (hl_deep_spec.py, bq_flow_spec.py): parity flips of the count of real
  ν > 1 of T_μ, refined by bisection; the state is extended to a deeper origin truncation by its exact local
  structure.
- The 2D local eigenproblem with a growth rate (the shift κ(μ) − κ(0) = iμ/D̂).
- The truncation artifacts and how they are identified.
- The WKB phase, its cut-off study and the uniformity test.
- Hou–Luo: solver, crossings (N = 8192/32768), stability (hl_stability.py), dip scaling.

## Table 1 — the eight resonances
| n | λ_n (method A) | λ_n (method B) | 1/(λ_n − 1) |
|---|---|---|---|
| 0 | 1.9205593 | 1.9205610 | 1.08630 |
| 1 | 1.3990961 | 1.3990960 | 2.50566 |
| 2 | 1.2523487 | 1.2523481 | 3.96277 |
| 3 | 1.1842533 | 1.1842530 | 5.42733 |
| 4 | 1.1449857 | 1.1449853 | 6.89723 |
| 5 | 1.1194738 | 1.1194739 | 8.37004 |
| 6 | 1.1015817 | 1.1015771–1.1015821 | 9.84429 |
| 7 | 1.0883384 | 1.0883386 | 11.32010 |

## Table 2 — unstable eigenvalues μ_k (perturbations ∝ e^{μτ}, τ = −ln(1−t))
| n | 2D: μ_k (linearization B; λ₇'s lowest from A) | Hou–Luo: μ_k |
|---|---|---|
| 1 | 0.37379 | 0.37385 |
| 2 | 0.55419, 0.22048 | 0.57394, 0.22677 |
| 3 | 0.63437, 0.37541, 0.15492 | 0.66213, 0.40156, 0.15996 |
| 4 | 0.68004, 0.45973, 0.28629, 0.11908 | 0.71143, 0.49689, 0.30871, 0.12333 |
| 5 | 0.70985, 0.51143, 0.36869, 0.23121, 0.09654 | 0.74250, 0.55390, 0.40264, 0.25069, 0.10010 |
| 6 | 0.73139, 0.54481, 0.42577, 0.30749, 0.19415, 0.08154 | 0.76371, 0.58776, 0.46858, 0.33755, 0.21095, 0.08416 |
| 7 | 0.74632, 0.56676, 0.46713, 0.36192, 0.26371, 0.16678, 0.06985 | 0.77911, 0.60295, 0.52201, 0.39978, 0.29125, 0.18202, 0.07254 |
| 8 | — | 0.79095, 0.587 ± 0.022i, 0.44373, 0.35277, 0.25583, 0.16004, 0.06372 |
